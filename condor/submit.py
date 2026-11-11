#!/usr/bin/env python3
"""
Prepare (and optionally submit) HTCondor jobs for the HSCP selector.

One job = one sample x one chunk of N ntuple files x one group of selections.
The job compiles libTools and the selector of its group in its scratch directory,
so nothing is shared between jobs and the worker only needs ROOT (see env.sh).

    python3 condor/submit.py --samples 'JetMET2024*' --tag jetmet                  # prepare only
    python3 condor/submit.py --dataset JetMET2024 --tag jetmet --submit             # same samples, and condor_submit
    python3 condor/submit.py --samples TestMET2024,TestMuon2024 --files-per-job 5 --selections-per-job 12

Samples, their ntuples and their binning come from cfg/datasets.json: names are matched
(fnmatch) with --samples, or every sample of a dataset taken with --dataset. The production
directory is --prod-dir, else $HSCP_PROD_DIR, else the prod_dir of the file.

Everything is written under condor/work/<tag>/:
    chunks/        one file list per chunk
    sel/g<k>/      selection file and generated selector of group k (macros/HSCPSelector.{h,C})
    jobs.txt       one line per job (the queue of hscp.sub)
    hscp.sub       the submit file
    out/           job outputs (after running)
    log/           Condor logs
Afterwards: python3 condor/merge.py condor/work/<tag>
"""

import argparse
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "python"))
from hscp import config, driver, selector  # noqa: E402

REPO = config.REPO


def _mk(d):
    os.makedirs(os.path.join(d, "macros"), exist_ok=True)
    return os.path.join(d, "HSCPpreSelections.csv")


# ----------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", default="",
                    help="comma-separated sample names or fnmatch patterns, e.g. 'JetMET2024*,Mu2024C'")
    ap.add_argument("--dataset", action="append", default=[], metavar="NAME", help="all samples of a dataset (repeatable)")
    ap.add_argument("--tag", required=True, help="name of the working directory condor/work/<tag>")
    ap.add_argument("--files-per-job", type=int, default=20)
    ap.add_argument("--selections-per-job", type=int, default=2,
                    help="selections per job; memory is about 1.7 GB per selection at the default binning")
    ap.add_argument("--selections", default=config.SELECTIONS_FILE)
    ap.add_argument("--prod-dir", help="production directory (default: $HSCP_PROD_DIR, else prod_dir of datasets.json)")
    ap.add_argument("--memory", default="4 GB", help="request_memory")
    ap.add_argument("--cpus", type=int, default=1, help="request_cpus, also HSCP_NTHREADS in the job")
    ap.add_argument("--runtime", type=int, default=7200, help="+MaxRuntime in seconds")
    ap.add_argument("--requirements", default="", help="extra Condor requirements expression")
    ap.add_argument("--submit", action="store_true", help="run condor_submit at the end")
    args = ap.parse_args()
    patterns = [p.strip() for p in args.samples.split(",") if p.strip()]
    if not patterns and not args.dataset:
        ap.error("give --samples and/or --dataset")

    work = os.path.join(REPO, "condor", "work", args.tag)
    if os.path.exists(os.path.join(work, "jobs.txt")):
        sys.exit(f"{work} already holds jobs.txt: choose another --tag or remove the directory")
    for d in ("chunks", "sel", "out", "log"):
        os.makedirs(os.path.join(work, d), exist_ok=True)

    # ---- samples ----
    try:
        cfg = config.load(prod_dir=args.prod_dir)
        selected = cfg.select(patterns, args.dataset)
    except config.ConfigError as e:
        sys.exit(f"submit: {e}")

    # ---- selection groups: one generated selector per group ----
    try:
        sels = selector.read_selections(args.selections)
    except selector.SelectionError as e:
        sys.exit(f"submit: {e}")
    groups = list(config.chunked(sels, args.selections_per_job))
    for g, group in enumerate(groups):
        gdir = os.path.join(work, "sel", f"g{g}")
        selector.write_selections(group, os.path.join(gdir, "HSCPpreSelections.csv") if os.path.isdir(gdir) else _mk(gdir))
        selector.generate(group, os.path.join(gdir, "macros"))
        shutil.copy(os.path.join(config.MACROS_DIR, "macro.cc"), os.path.join(gdir, "macros", "macro.cc"))
    print(f"{len(sels)} selections in {len(groups)} group(s) of up to {args.selections_per_job}")

    # ---- chunks and job list ----
    jobs = []
    for name in selected:
        binning = cfg.binning(name)
        try:
            files = cfg.files(name)
        except config.ConfigError as e:
            sys.exit(f"submit: {e}")
        chunks = list(config.chunked(files, args.files_per_job))
        for c, chunk in enumerate(chunks):
            chunk_path = os.path.join(work, "chunks", f"{name}_c{c}.txt")
            with open(chunk_path, "w") as f:
                f.write("\n".join(chunk) + "\n")
            for g in range(len(groups)):
                out = f"{name}_g{g}_c{c}.root"
                nosel = 1 if g == 0 else 0          # no-selection histograms from one group only
                b = " ".join(driver._num(x) for x in binning)
                jobs.append(f"{name} {chunk_path} {g} {out} {nosel} {b} {c}")
        print(f"{name:<32} {len(files):5d} files -> {len(chunks)} chunk(s) x {len(groups)} group(s)")

    with open(os.path.join(work, "jobs.txt"), "w") as f:
        f.write("\n".join(jobs) + "\n")

    # ---- submit file ----
    req = f"requirements = {args.requirements}\n" if args.requirements else ""
    sub = f"""# generated by condor/submit.py; queue read from jobs.txt
universe              = vanilla
executable            = {REPO}/condor/runJob.sh
initialdir            = {work}
arguments             = $(sample) $(chunk) $(group) $(out) $(nosel) $(ptcut) $(etabins) $(ihbins) $(pbins) $(massbins) $(fpixbins)
environment           = "HSCP_NTHREADS={args.cpus}"

should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files  = {REPO}/env.sh, {REPO}/compile.sh, {REPO}/inc, {REPO}/src, {work}/sel/g$(group)/macros, $(chunk)
transfer_output_files = $(out)
transfer_output_remaps = "$(out) = out/$(out)"

output                = log/$(sample)_g$(group)_c$(chunkid).out
error                 = log/$(sample)_g$(group)_c$(chunkid).err
log                   = log/condor.log

request_cpus          = {args.cpus}
request_memory        = {args.memory}
+MaxRuntime           = {args.runtime}
max_retries           = 2
on_exit_hold          = (ExitCode != 0)
{req}
queue sample, chunk, group, out, nosel, ptcut, etabins, ihbins, pbins, massbins, fpixbins, chunkid from jobs.txt
"""
    sub_path = os.path.join(work, "hscp.sub")
    with open(sub_path, "w") as f:
        f.write(sub)

    print(f"\n{len(jobs)} jobs prepared in {work}")
    print(f"submit file: {sub_path}")
    if args.submit:
        subprocess.run(["condor_submit", sub_path], check=True, cwd=work)
    else:
        print(f"to submit:   condor_submit {sub_path}")


if __name__ == "__main__":
    main()
