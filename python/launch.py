#!/usr/bin/env python3
"""
Run the selector on samples, locally.

    ./launch.sh 'JetMET2024*' Mu2024C              # every era of JetMET2024, then Mu2024C
    ./launch.sh --dataset Stop --dataset Stau      # every sample of two datasets
    ./launch.sh --list '*'                         # all samples: type, dataset, binning
    ./launch.sh 'Gluino_Run3_MET_madgraph_*' --jobs 3 --out-dir output/gluino
    ./launch.sh TestMET2024 --dry-run

Samples are selected by fnmatch patterns against the names of cfg/datasets.json, or by
dataset (a group of samples, e.g. JetMET2024 for the eras C to I). Their ntuples and binning
come from the same file, under $HSCP_PROD_DIR (or --prod-dir, or the prod_dir of the file).
The selector is generated from the selection file and compiled once, then every sample runs
in its own root process on an explicit file list written next to its output:
<out-dir>/<sample>[_<label>].root, .log and .files.

launch.sh sources env.sh (ROOT) and calls this script.
"""

import argparse
import os
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hscp import config, driver, selector  # noqa: E402


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("patterns", nargs="*", help="sample names or fnmatch patterns")
    ap.add_argument("--dataset", action="append", default=[], metavar="NAME", help="all samples of a dataset (repeatable)")
    ap.add_argument("--list", action="store_true", help="print the matching samples and exit")
    ap.add_argument("--dry-run", action="store_true", help="print the root commands and exit")
    ap.add_argument("--out-dir", default=os.path.join(config.REPO, "output"),
                    help="where the outputs and logs go (default: output/)")
    ap.add_argument("--label", default="", help="suffix of the output files: <sample>_<label>.root")
    ap.add_argument("--jobs", type=int, default=1, help="samples run at the same time (each with 1 thread)")
    ap.add_argument("--threads", type=int, help="implicit-MT threads per sample (default: 1 if --jobs > 1, else ROOT's 4)")
    ap.add_argument("--keep-going", action="store_true", help="continue with the next sample after a failure")
    ap.add_argument("--selections", default=config.SELECTIONS_FILE, help="selection file (default: cfg/HSCPpreSelections.csv)")
    ap.add_argument("--prod-dir", help="production directory (default: $HSCP_PROD_DIR, else prod_dir of datasets.json)")
    ap.add_argument("--keep-build", action="store_true", help="keep the compiled selector (HSCPSelector_C*) afterwards")
    args = ap.parse_args(argv)
    if not args.patterns and not args.dataset:
        ap.error("give at least one sample name, pattern or --dataset (try: --list '*')")
    return args


def main(argv=None):
    args = parse_args(argv)
    try:
        cfg = config.load(prod_dir=args.prod_dir)
        names = cfg.select(args.patterns, args.dataset)
    except config.ConfigError as e:
        sys.exit(f"launch: {e}")

    if args.list:
        for n in names:
            print(f"{n:<34} {cfg.type(n):<6} {cfg.dataset(n):<16} "
                  f"{' '.join(driver._num(x) for x in cfg.binning(n))}")
        return 0

    files = {}
    for n in names:                                    # fail before compiling anything
        try:
            files[n] = cfg.files(n)
        except config.ConfigError as e:
            sys.exit(f"launch: {e}")

    os.makedirs(args.out_dir, exist_ok=True)
    suffix = f"_{args.label}" if args.label else ""
    base = {n: os.path.join(os.path.abspath(args.out_dir), f"{n}{suffix}") for n in names}
    plan = [(n, base[n] + ".files", base[n] + ".root", base[n] + ".log") for n in names]
    for n, flist, _, _ in plan:
        with open(flist, "w") as f:
            f.write("\n".join(files[n]) + "\n")

    if args.dry_run:
        print("cd macros &&", " ".join(driver.build_command()[:-1]), repr(driver.build_command()[-1]))
        for n, flist, out, _ in plan:
            print("cd macros &&", " ".join(driver.root_command(n, flist, out, True, cfg.binning(n))[:-1]),
                  repr(driver.macro_call(n, flist, out, True, cfg.binning(n))))
        return 0

    if shutil.which("root") is None:
        sys.exit("launch: root not found in PATH (source env.sh, or use ./launch.sh)")

    try:
        sels = selector.read_selections(args.selections)
    except selector.SelectionError as e:
        sys.exit(f"launch: {e}")
    selector.generate(sels)
    print(f"{len(sels)} selections; compiling the selector ...", flush=True)
    build_log = os.path.join(args.out_dir, "build.log")
    if driver.build_selector(log=build_log) != 0:
        sys.exit(f"launch: the selector does not compile, see {build_log}")

    threads = args.threads if args.threads is not None else (1 if args.jobs > 1 else None)
    failed = []

    def run(item):
        n, flist, out, log = item
        rc = driver.run_dataset(n, flist, out, True, cfg.binning(n), log, threads=threads)
        print(f"{'ok  ' if rc == 0 else 'FAIL'}  {n:<34} -> {out if rc == 0 else log}", flush=True)
        return n, rc

    if args.jobs > 1:
        with ThreadPoolExecutor(max_workers=args.jobs) as pool:
            for n, rc in pool.map(run, plan):
                if rc:
                    failed.append(n)
    else:
        for item in plan:
            n, rc = run(item)
            if rc:
                failed.append(n)
                if not args.keep_going:
                    break

    if not args.keep_build:
        driver.clean_build()
    print(f"\n{len(plan) - len(failed)} / {len(plan)} samples done" + (f", failed: {' '.join(failed)}" if failed else ""))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
