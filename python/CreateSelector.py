#!/usr/bin/env python3
"""
Generate HSCPSelector.h and HSCPSelector.C from the templates and a selection file.

    python3 python/CreateSelector.py                                   # cfg/HSCPpreSelections.csv -> macros/
    python3 python/CreateSelector.py --selections sel.csv --out-dir d  # one selector per Condor group

launch.py and condor/submit.py do this themselves; this script is for a manual run.
The code lives in python/hscp/selector.py.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hscp import config, selector  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selections", default=config.SELECTIONS_FILE)
    ap.add_argument("--out-dir", default=config.MACROS_DIR)
    args = ap.parse_args()
    try:
        sels = selector.read_selections(args.selections)
    except selector.SelectionError as e:
        sys.exit(str(e))
    print("Selections are:")
    for label, _ in sels:
        print("  " + label)
    h, c = selector.generate(sels, args.out_dir)
    print(f"Written {h} and {c}")


if __name__ == "__main__":
    main()
