#!/usr/bin/env python3
"""
Make the validation plots of macros/CombineHistos.C from the input files named in cfg/plots.json.

    python3 python/plots.py                          # every plot of cfg/plots.json
    python3 python/plots.py --only 'nPV,nHSCP'       # some of them (fnmatch patterns on the plot names)
    python3 python/plots.py --list                   # the plots, their inputs, and which inputs are missing
    python3 python/plots.py --dry-run                # the root commands, nothing run
    python3 python/plots.py --config cfg/plots_V2.json --output-dir outputDisplay/V2

cfg/plots.json names the input files once ("inputs": data, qcd, ttbar, ...), then lists the macro
functions to call with them ("plots"); see python/hscp/plots.py for the format. Each plot runs in its
own root process, started in macros/: one failing plot does not stop the others. A plot whose input
files are missing is skipped. Outputs (pdf, root) go to output_dir (default outputDisplay/), and each
plot's log to <output_dir>/logs/<name>.log.

Needs root in the PATH: source env.sh first.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hscp import driver, plots  # noqa: E402


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", default=plots.PLOTS_FILE, help="plot file (default: %(default)s)")
    ap.add_argument("--only", default="", help="comma-separated plot names or fnmatch patterns")
    ap.add_argument("--output-dir", help="instead of output_dir of the plot file")
    ap.add_argument("--list", action="store_true", help="list the plots and their inputs")
    ap.add_argument("--dry-run", action="store_true", help="print the root commands only")
    ap.add_argument("--jobs", type=int, default=1, help="plots run at the same time (default: 1)")
    return ap.parse_args(argv)


def shell_quote(argv):
    import shlex
    return " ".join(shlex.quote(a) for a in argv)


def run_one(cfg, plot, labels, log):
    t0 = time.time()
    with open(log, "w") as out:
        rc = subprocess.run(plots.root_command(cfg, plot, labels), cwd=os.path.dirname(cfg.macro),
                            stdout=out, stderr=subprocess.STDOUT).returncode
    return rc, time.time() - t0


def describe(rc, log):
    """Why a plot failed: ROOT catches a crash (*** Break *** in the log) and exits with an unrelated code."""
    try:
        with open(log, errors="replace") as f:
            crash = next((line.split("***")[-1].strip() for line in f if "*** Break ***" in line), None)
    except OSError:
        crash = None
    if crash:
        return f"crashed ({crash}), often a histogram missing from an input file"
    return {2: "exit code 2: the macro does not load",
            3: "exit code 3: the call does not compile (function name or arguments)",
            4: "exit code 4: an input is missing or unusable, see the ERROR line of the log"}.get(rc) or driver.describe_exit(rc)


def missing_histograms(cfg, plot, names, ROOT):
    """The names read by the function that none of its input files holds (top-level keys)."""
    keys = set()
    for n in dict.fromkeys(cfg.input_names(plot)):
        path = cfg.inputs[n]
        if "%d" in path:
            continue
        f = ROOT.TFile.Open(path)
        if not f or f.IsZombie():
            continue
        keys |= {k.GetName() for k in f.GetListOfKeys()}
        f.Close()
    return [h for h in names if h not in keys]


def main(argv=None):
    args = parse_args(argv)
    try:
        cfg = plots.load(args.config)
        if args.output_dir:
            cfg.output_dir = os.path.abspath(args.output_dir)
        chosen = cfg.select([p.strip() for p in args.only.split(",") if p.strip()])
    except plots.PlotError as e:
        sys.exit(f"plots: {e}")
    try:
        with open(cfg.macro) as f:
            macro_text = f.read()
    except OSError as e:
        sys.exit(f"plots: cannot read {cfg.macro}: {e.strerror}")
    errors = plots.check_calls(cfg, chosen, macro_text)
    if errors:
        sys.exit("plots: " + "\n".join(errors))
    labels = {p.name: cfg.legend_labels(p, plots.signature(macro_text, p.function)) for p in chosen}

    if args.list:
        for p in chosen:
            missing = dict(cfg.missing_inputs(p))
            shown = ", ".join(a if isinstance(a, str) else "null" if a is None else cfg.cpp_arg(a) for a in p.args)
            print(f"{p.name}  ->  {p.function}({shown})")
            for n in dict.fromkeys(cfg.input_names(p)):
                print(f"    {n:<18} {cfg.inputs[n]}" + ("   MISSING" if n in missing else "")
                      + (f"   label: {cfg.labels[n]}" if n in cfg.labels else ""))
            if p.note:
                print(f"    note: {p.note}")
        return 0

    logs = os.path.join(cfg.output_dir, "logs")
    if args.dry_run:
        print(f"cd {os.path.dirname(cfg.macro)}")
        for p in chosen:
            print(shell_quote(plots.root_command(cfg, p, labels[p.name])) + f" > {os.path.join(logs, p.name + '.log')} 2>&1")
        return 0

    if shutil.which("root") is None:
        sys.exit("plots: root not found in PATH (source env.sh first)")
    os.makedirs(logs, exist_ok=True)

    runnable, skipped = [], []
    for p in chosen:
        missing = cfg.missing_inputs(p)
        if missing:
            skipped.append(p)
            print(f"skipped    {p.name:<28} missing input(s): " + ", ".join(f"{n} ({path})" for n, path in missing))
        else:
            runnable.append(p)

    # histograms the functions read by a fixed name, missing from all their inputs: warned before running
    try:
        import ROOT
        ROOT.gROOT.SetBatch(True)
    except ImportError:
        ROOT = None
    if ROOT is not None:
        for p in runnable:
            absent = missing_histograms(cfg, p, plots.histogram_names(macro_text, p.function), ROOT)
            if absent:
                print(f"warning    {p.name:<28} {len(absent)} histogram(s) not in its inputs: " + ", ".join(absent[:4])
                      + (", ..." if len(absent) > 4 else ""))

    failed = []

    def task(p):
        print(f"processing {p.name} ...", flush=True)
        rc, dt = run_one(cfg, p, labels[p.name], os.path.join(logs, p.name + ".log"))
        if rc == 0:
            print(f"completed  {p.name:<28} {dt:5.0f} s", flush=True)
        else:
            failed.append(p)
            log = os.path.join(logs, p.name + ".log")
            print(f"FAILED     {p.name:<28} {dt:5.0f} s  {describe(rc, log)}, see {log}", flush=True)

    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        list(pool.map(task, runnable))

    done = len(runnable) - len(failed)
    print(f"\n{done} / {len(chosen)} plot(s) made in {cfg.output_dir}"
          + (f", {len(failed)} failed" if failed else "") + (f", {len(skipped)} skipped" if skipped else ""))
    return 0 if not failed and not skipped else 1


if __name__ == "__main__":
    sys.exit(main())
