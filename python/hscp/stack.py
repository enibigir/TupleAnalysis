"""
Helpers of python/stack.py that need no ROOT: which samples to draw, which histograms, colours, file names.
"""

import fnmatch
import re
from collections import namedtuple

from .plots import PlotError

# stacked backgrounds, in stacking order, and signals (lines), when the input sets no "color"
BACKGROUND_COLORS = ["kGreen-4", "kRed", "kRed+2", "kBlue-7", "kOrange+1", "kMagenta-7", "kCyan+1", "kYellow+2"]
SIGNAL_COLORS = ["kOrange+8", "kViolet+1", "kAzure+2", "kPink+1", "kSpring-6"]

Samples = namedtuple("Samples", "data backgrounds signals")


def split_names(text):
    """'a, b' -> ['a', 'b']; None stays None (not given on the command line)."""
    if text is None:
        return None
    return [n.strip() for n in text.split(",") if n.strip()]


def samples(cfg, data=None, backgrounds=None, signals=None):
    """The data, backgrounds and signals to draw: the command line (data: a name, "" for none; lists of names)
    over the "stack" block of the configuration."""
    s = cfg.stack or {}
    data = s.get("data") if data is None else (data or None)
    backgrounds = s.get("backgrounds", []) if backgrounds is None else backgrounds
    signals = s.get("signals", []) if signals is None else signals
    names = ([data] if data else []) + backgrounds + signals
    for n in names:
        if n not in cfg.inputs:
            raise PlotError(f"'{n}' is not an input of {cfg.source} (inputs: {', '.join(sorted(cfg.inputs))})")
        if "%d" in cfg.inputs[n]:
            raise PlotError(f"input {n} is a pattern (%d): give one file per sample")
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        raise PlotError(f"sample(s) given twice: {', '.join(dup)}")
    if not backgrounds and not data and not signals:
        raise PlotError('no sample to draw: give --data/--bkg/--signal, or a "stack" block in ' + cfg.source)
    return Samples(data, backgrounds, signals)


def expand(patterns, keys):
    """The names of keys (in their order) matching any fnmatch pattern; each pattern must match one."""
    chosen = [k for k in keys if any(fnmatch.fnmatchcase(k, p) for p in patterns)]
    unmatched = [p for p in patterns if not any(fnmatch.fnmatchcase(k, p) for k in keys)]
    if unmatched:
        raise PlotError("no histogram named " + ", ".join(unmatched))
    return chosen


_COLOR = re.compile(r"^(k[A-Z][a-z]+)\s*([+-]\s*\d+)?$")


def parse_color(c):
    """A ROOT colour number, or ("kRed", 2) for "kRed+2"."""
    if type(c) is int:
        return c
    m = _COLOR.match(str(c).strip())
    if not m:
        raise PlotError(f"colour {c!r}: give a number or a ROOT name such as kRed+2")
    return m.group(1), int(m.group(2).replace(" ", "")) if m.group(2) else 0


def colors(cfg, names, palette):
    """{name: colour spec}: the "color" of the input, else the next one of the palette."""
    out, free = {}, iter(palette * (len(names) // len(palette) + 1))
    for n in names:
        out[n] = parse_color(cfg.colors[n]) if n in cfg.colors else parse_color(next(free))
    return out


def file_name(hist):
    """A file name for a histogram name."""
    return re.sub(r"[^\w.+-]", "_", hist)


def legend_label(cfg, name, scale=1.0):
    label = cfg.labels.get(name, name)
    return f"{label} #times {scale:g}" if scale != 1 else label
