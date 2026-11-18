#!/usr/bin/env python3
"""
Weight the outputs of simulated samples by cross-section and add them into one file.

    python3 python/Weight.py --dataset QCD2024 --input-dir condor/work/qcd/merged --dry-run   # show the factors
    python3 python/Weight.py --dataset QCD2024 --input-dir condor/work/qcd/merged \\
            -o output/QCD2024_V16/QCD2024_mu_V16p2_weighted.root
    python3 python/Weight.py --samples 'JetMET2024*' --input-dir condor/work/jetmet/merged -o JetMET2024.root   # data: added only

Each selected sample is read from <input-dir>/<sample>.root: the merged/ directory of condor/merge.py,
or the --out-dir of launch.sh. Every histogram of a simulated sample is multiplied by

    xsec [pb] x luminosity [pb^-1] / sum of weights

and the samples are added. xsec comes from cfg/datasets.json, the luminosity from --lumi or the
'luminosity' of datasets.json. The sum of weights is 'sum_weights' of the sample in datasets.json if given,
else bin 1 of its EventCutflow histogram, the sum of weightPU over every event the selector read: give
sum_weights if the ntuples were skimmed, or if the generator weights are not all 1 (negative weights).
Profiles are added with these factors as weights; their mean values are not scaled. Data samples are
added with factor 1. Data and simulation cannot be mixed in one output.

The default output is <input-dir>/<dataset>_weighted.root when all the samples belong to one dataset.
Histograms are processed one name at a time, so memory stays small whatever the size of the files.
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hscp import config, weight  # noqa: E402


def parse_args(argv):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--samples", default="", help="comma-separated sample names or fnmatch patterns")
    ap.add_argument("--dataset", action="append", default=[], metavar="NAME", help="all samples of a dataset (repeatable)")
    ap.add_argument("--input-dir", required=True, help="directory holding <sample>.root for each sample")
    ap.add_argument("-o", "--output", help="output file (default: <input-dir>/<dataset>_weighted.root for a single dataset)")
    ap.add_argument("--lumi", type=float, help="integrated luminosity in pb^-1 (default: 'luminosity' of datasets.json)")
    ap.add_argument("--config", default=config.DATASETS_FILE, help="datasets file (default: %(default)s)")
    ap.add_argument("--dry-run", action="store_true", help="print the factors only")
    ap.add_argument("-f", "--force", action="store_true", help="overwrite an existing output")
    return ap.parse_args(argv)


# ---------------------------------------------------------------------- ROOT part
def _keys(directory):
    """{name: class name} of a directory, in file order, one entry per name (highest cycle)."""
    out = {}
    for k in directory.GetListOfKeys():
        out.setdefault(k.GetName(), k.GetClassName())
    return out


def _read(ROOT, directory, name):
    obj = directory.Get(name)
    if obj:
        ROOT.SetOwnership(obj, True)       # TH1.AddDirectory(False): nobody else owns it
    return obj


def _scaled(ROOT, obj, f):
    """A copy of obj with every fill weighted by f."""
    if f == 1.0:
        return obj
    if obj.InheritsFrom("TProfile") or obj.InheritsFrom("TProfile2D") or obj.InheritsFrom("TProfile3D"):
        # TProfile::Scale would scale the mean values; Add(p, f) weights the entries instead
        c = obj.Clone()
        c.Reset()
        c.Add(obj, f)
        return c
    obj.Scale(f)
    return obj


def merge_directory(ROOT, sources, out_dir, stats):
    """Add the objects of sources [(directory, factor)] into out_dir, recursing into sub-directories."""
    keys = [_keys(d) for d, _ in sources]
    names = list(dict.fromkeys(n for k in keys for n in k))
    for name in names:
        cls = next(k[name] for k in keys if name in k)
        tclass = ROOT.TClass.GetClass(cls)
        present = [(d, f) for (d, f), k in zip(sources, keys) if name in k]
        if len(present) < len(sources):
            stats["partial"].append(name)
        if tclass and tclass.InheritsFrom("TDirectory"):
            sub = out_dir.mkdir(name)
            merge_directory(ROOT, [(d.Get(name), f) for d, f in present], sub, stats)
            continue
        if tclass and tclass.InheritsFrom("TTree"):
            stats["skipped"].append(name)
            continue
        if not (tclass and tclass.InheritsFrom("TH1")):
            out_dir.WriteTObject(_read(ROOT, present[0][0], name), name)
            stats["copied"] += 1
            continue
        total = None
        for d, f in present:
            h = _scaled(ROOT, _read(ROOT, d, name), f)
            if total is None:
                total = h
            elif not total.Add(h):
                raise RuntimeError(f"cannot add {name}: different binning in the samples")
        out_dir.WriteTObject(total, name)
        stats["histograms"] += 1


def sum_weights_from(ROOT, path):
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise weight.WeightError(f"cannot open {path}")
    h = f.Get(weight.NORM_HISTO)
    if not h:
        raise weight.WeightError(f"{path} has no {weight.NORM_HISTO} histogram: it is written with the "
                                 "no-selection histograms (by selection group 0 on Condor; launch.sh always writes them)")
    value = h.GetBinContent(1)
    f.Close()
    return value


# ---------------------------------------------------------------------- main
def main(argv=None):
    args = parse_args(argv)
    patterns = [p.strip() for p in args.samples.split(",") if p.strip()]
    if not patterns and not args.dataset:
        sys.exit("Weight: give --samples and/or --dataset")
    try:
        cfg = config.load(args.config, prod_dir="/")      # the production directory is not used here
        samples = cfg.select(patterns, args.dataset)
    except config.ConfigError as e:
        sys.exit(f"Weight: {e}")

    paths = {s: os.path.join(args.input_dir, f"{s}.root") for s in samples}
    missing = [p for p in paths.values() if not os.path.isfile(p)]
    if missing:
        sys.exit("Weight: missing input file(s):\n  " + "\n  ".join(missing))
    output = args.output or weight.default_output(cfg, samples, args.input_dir)
    if not output and not args.dry_run:
        sys.exit("Weight: the samples belong to several datasets: give -o")

    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.TH1.AddDirectory(False)
    lumi = args.lumi if args.lumi is not None else cfg.luminosity
    try:
        entries = weight.plan(cfg, samples, lumi, lambda s: sum_weights_from(ROOT, paths[s]))
    except weight.WeightError as e:
        sys.exit(f"Weight: {e}")
    print(weight.table(entries, lumi))
    if args.dry_run:
        return

    if os.path.exists(output) and not args.force:
        sys.exit(f"Weight: {output} exists: use -f to overwrite it")
    if os.path.abspath(output) in {os.path.abspath(p) for p in paths.values()}:
        sys.exit(f"Weight: the output {output} is one of the inputs")
    os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)

    files = [ROOT.TFile.Open(paths[e.sample]) for e in entries]
    for e, f in zip(entries, files):
        if not f or f.IsZombie():
            sys.exit(f"Weight: cannot open {paths[e.sample]}")
    tmp = output + ".part"
    out = ROOT.TFile.Open(tmp, "RECREATE")
    stats = {"histograms": 0, "copied": 0, "skipped": [], "partial": []}
    try:
        merge_directory(ROOT, [(f, e.factor) for f, e in zip(files, entries)], out, stats)
    except RuntimeError as err:
        out.Close()
        os.remove(tmp)
        sys.exit(f"Weight: {err}")
    out.Close()
    for f in files:
        f.Close()
    os.replace(tmp, output)

    print(f"\n{stats['histograms']} histograms written to {output}"
          + (f", {stats['copied']} other objects copied" if stats["copied"] else ""))
    if stats["skipped"]:
        print(f"trees left out: {', '.join(stats['skipped'])}")
    if stats["partial"]:
        print(f"warning: {len(stats['partial'])} object(s) missing from some samples, e.g. {stats['partial'][0]}")


if __name__ == "__main__":
    main()
