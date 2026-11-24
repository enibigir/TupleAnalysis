"""
Cross-section weighting of the selector outputs (python/Weight.py). No ROOT here.

A simulated sample is scaled by

    factor = xsec [pb] x luminosity [pb^-1] / sum of weights

where the sum of weights is 'sum_weights' of the sample in cfg/datasets.json if given, else bin 1 of
the EventCutflow histogram of its output: the sum of weightPU over every event the selector read.
That default is right only if the ntuples hold every event produced (no skim) and the generator weights
are all 1; otherwise give sum_weights. Data samples get the factor 1 (they are only added).
"""

import os
from collections import namedtuple

NORM_HISTO = "EventCutflow"       # bin 1: sum of weightPU over all events read (filled with the no-selection histograms)

Entry = namedtuple("Entry", "sample type xsec sum_weights source factor")


class WeightError(ValueError):
    pass


def factor(xsec, luminosity, sum_weights):
    if not sum_weights > 0:
        raise WeightError(f"sum of weights {sum_weights} is not positive")
    return xsec * luminosity / sum_weights


def available(paths, skip_missing=False):
    """Split {sample: input file} into (samples whose file exists, missing files). A missing file is an
    error unless skip_missing; at least one sample must be left."""
    present = [s for s, p in paths.items() if os.path.isfile(p)]
    missing = [p for s, p in paths.items() if not os.path.isfile(p)]
    if missing and not skip_missing:
        raise WeightError("missing input file(s) (--skip-missing leaves those samples out):\n  " + "\n  ".join(missing))
    if not present:
        raise WeightError("none of the input files exists:\n  " + "\n  ".join(missing))
    return present, missing


def plan(cfg, samples, luminosity, sum_weights_of):
    """The weight of each sample, as a list of Entry.

    luminosity in pb^-1 (None if not given). sum_weights_of(sample) returns the sum of weights read from
    the output file of a sample; it is called only for simulated samples without sum_weights in datasets.json.
    """
    if not samples:
        raise WeightError("no sample selected")
    types = {cfg.type(s) for s in samples}
    if "data" in types and len(types) > 1:
        raise WeightError("data and simulation cannot be added into one file: select them separately")
    if types == {"data"}:
        return [Entry(s, "data", None, None, "", 1.0) for s in samples]
    missing = [s for s in samples if cfg.xsec(s) is None]
    if missing:
        raise WeightError(f"no 'xsec' in {cfg.path} for: {', '.join(missing)}")
    if luminosity is None:
        raise WeightError(f"no luminosity: give --lumi or 'luminosity' (pb^-1) in {cfg.path}")
    entries = []
    for s in samples:
        sw, source = cfg.sum_weights(s), "datasets.json"
        if sw is None:
            sw, source = sum_weights_of(s), NORM_HISTO
        try:
            f = factor(cfg.xsec(s), luminosity, sw)
        except WeightError as e:
            raise WeightError(f"{s}: {e}")
        entries.append(Entry(s, cfg.type(s), cfg.xsec(s), sw, source, f))
    return entries


def default_output(cfg, samples, input_dir):
    """<input_dir>/<dataset>_weighted.root if every sample belongs to the same dataset, else None."""
    datasets = {cfg.dataset(s) for s in samples}
    if len(datasets) == 1 and "" not in datasets:
        return os.path.join(input_dir, f"{datasets.pop()}_weighted.root")
    return None


def table(entries, luminosity):
    """The plan as text, one line per sample."""
    lines = [f"{'sample':<32} {'xsec [pb]':>12} {'sum of weights':>16} {'from':<14} {'factor':>12}"]
    for e in entries:
        if e.type == "data":
            lines.append(f"{e.sample:<32} {'-':>12} {'-':>16} {'data':<14} {e.factor:12.6g}")
        else:
            lines.append(f"{e.sample:<32} {e.xsec:12.6g} {e.sum_weights:16.6g} {e.source:<14} {e.factor:12.6g}")
    if luminosity is not None and any(e.type != "data" for e in entries):
        lines.append(f"luminosity: {luminosity:g} pb^-1")
    return "\n".join(lines)

