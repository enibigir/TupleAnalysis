#!/usr/bin/env python3
"""
Check the outputs of a condor/submit.py working directory and merge them with hadd.

    python3 condor/merge.py condor/work/<tag>            # check, then merge the complete datasets
    python3 condor/merge.py condor/work/<tag> --check    # check only
    python3 condor/merge.py condor/work/<tag> --force    # merge even if chunks are missing

For each sample, every chunk of every selection group must be present in out/.
Missing jobs are written to <work>/resubmit.txt; submit them again with
    condor_submit <work>/hscp.sub -queue 'dataset, ... from resubmit.txt'   (same column order)
or simply:  cd <work> && cp resubmit.txt jobs.txt && condor_submit hscp.sub

Complete samples are merged into <work>/merged/<sample>.root with "hadd -T"
(histograms only). Histograms present in several chunks are summed; the groups contribute
disjoint histogram names, and the no-selection histograms come from group 0 only.
Needs hadd in the PATH (source env.sh).
"""

import argparse
import os
import shutil
import subprocess
import sys


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("work")
    ap.add_argument("--check", action="store_true", help="only report, do not merge")
    ap.add_argument("--force", action="store_true", help="merge incomplete samples too")
    ap.add_argument("--min-size", type=int, default=1024, help="bytes below which an output is treated as missing")
    args = ap.parse_args()

    work = os.path.abspath(args.work)
    jobs_path = os.path.join(work, "jobs.txt")
    if not os.path.isfile(jobs_path):
        sys.exit(f"{jobs_path} not found")

    # sample -> list of (output file, job line); job line: sample chunk group out nosel binning... chunkid
    per_sample = {}
    with open(jobs_path) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 5:
                continue
            per_sample.setdefault(parts[0], []).append((parts[3], line.rstrip("\n")))

    missing_lines = []
    complete = []
    for key, entries in sorted(per_sample.items()):
        missing = [job for out, job in entries
                   if not os.path.isfile(os.path.join(work, "out", out))
                   or os.path.getsize(os.path.join(work, "out", out)) < args.min_size]
        status = "complete" if not missing else f"{len(missing)} / {len(entries)} missing"
        print(f"{key:<32} {len(entries):4d} jobs  {status}")
        missing_lines += missing
        if not missing or args.force:
            complete.append(key)

    resubmit = os.path.join(work, "resubmit.txt")
    if missing_lines:
        with open(resubmit, "w") as f:
            f.write("\n".join(missing_lines) + "\n")
        print(f"\n{len(missing_lines)} missing job(s) written to {resubmit}")
    elif os.path.exists(resubmit):
        os.remove(resubmit)

    if args.check:
        return 0 if not missing_lines else 1

    if shutil.which("hadd") is None:
        sys.exit("hadd not found in PATH: source env.sh first")

    merged_dir = os.path.join(work, "merged")
    os.makedirs(merged_dir, exist_ok=True)
    failed = 0
    for key in complete:
        inputs = [os.path.join(work, "out", out) for out, _ in per_sample[key]
                  if os.path.isfile(os.path.join(work, "out", out))]
        target = os.path.join(merged_dir, f"{key}.root")
        print(f"\nhadd -T -f {target}  ({len(inputs)} files)")
        r = subprocess.run(["hadd", "-T", "-f", target] + inputs, stdout=subprocess.DEVNULL)
        if r.returncode != 0:
            failed += 1
            print(f"  FAILED (exit {r.returncode})")
    print(f"\n{len(complete) - failed} sample(s) merged into {merged_dir}" + (f", {failed} failed" if failed else ""))
    return 1 if failed or missing_lines else 0


if __name__ == "__main__":
    sys.exit(main())
