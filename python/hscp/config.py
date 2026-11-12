"""
cfg/datasets.json: the samples of the analysis, where their ntuples are, and the binning.

    {
      "prod_dir": "/scratch/.../HSCP_prod/",
      "binning":  {"ptcut": 70, "etabins": 96, "ihbins": 400, "pbins": 1000, "massbins": 400, "fpixbins": 20},
      "samples": {
        "JetMET2024C": {"list": "JetMET2024/V12p310.txt", "type": "data", "dataset": "JetMET2024"},
        "Stop_Run3_MET_madgraph_700": {"file": "SIGNAL/V21p0/HSCP-Stop_Par-M-700_merged.root",
                                       "type": "signal", "dataset": "Stop", "mass": 700,
                                       "binning": {"etabins": 48}}
      }
    }

Only Python reads this file. macros/macro.cc receives an explicit file list and the binning as arguments.
"""

import fnmatch
import json
import os
from collections import namedtuple

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CFG_DIR = os.path.join(REPO, "cfg")
MACROS_DIR = os.path.join(REPO, "macros")
INC_DIR = os.path.join(REPO, "inc")
DATASETS_FILE = os.path.join(CFG_DIR, "datasets.json")
SELECTIONS_FILE = os.path.join(CFG_DIR, "selections.json")

Binning = namedtuple("Binning", "ptcut etabins ihbins pbins massbins fpixbins")
TYPES = ("data", "mc", "signal")
ENTRY_KEYS = {"list", "file", "type", "dataset", "mass", "binning", "_comment"}


class ConfigError(ValueError):
    pass


def full_path(path, prod):
    """A path of datasets.json: relative to the production directory, unless absolute or a URL."""
    if path.startswith("/") or "://" in path:
        return path
    return prod + path


def chunked(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def _binning(values, where):
    try:
        return Binning(float(values["ptcut"]), *(int(values[k]) for k in Binning._fields[1:]))
    except KeyError as e:
        raise ConfigError(f"{where}: binning key {e} missing")
    except (TypeError, ValueError) as e:
        raise ConfigError(f"{where}: bad binning value ({e})")


class Config:
    """The parsed datasets.json."""

    def __init__(self, doc, path=DATASETS_FILE, prod_dir=None):
        self.path = path
        if not isinstance(doc.get("samples"), dict) or not doc["samples"]:
            raise ConfigError(f"{path}: 'samples' must be a non-empty object")
        if not isinstance(doc.get("binning"), dict):
            raise ConfigError(f"{path}: 'binning' object missing")
        self.default_binning = _binning(doc["binning"], f"{path}: binning")
        d = prod_dir or os.environ.get("HSCP_PROD_DIR") or doc.get("prod_dir")
        if not d:
            raise ConfigError(f"{path}: 'prod_dir' missing and HSCP_PROD_DIR not set")
        self.prod_dir = d if d.endswith("/") else d + "/"
        self.samples = {}
        for name, entry in doc["samples"].items():
            self.samples[name] = self._check(name, entry)

    def _check(self, name, entry):
        where = f"{self.path}: sample {name}"
        if not isinstance(entry, dict):
            raise ConfigError(f"{where}: must be an object")
        unknown = set(entry) - ENTRY_KEYS
        if unknown:
            raise ConfigError(f"{where}: unknown key(s) {sorted(unknown)}")
        if ("list" in entry) == ("file" in entry):
            raise ConfigError(f"{where}: exactly one of 'list' or 'file' is required")
        if entry.get("type") not in TYPES:
            raise ConfigError(f"{where}: 'type' must be one of {TYPES}")
        if "binning" in entry:
            merged = dict(self.default_binning._asdict()); merged.update(entry["binning"])
            _binning(merged, where)
        return entry

    # ---- queries ----
    @property
    def names(self):
        return list(self.samples)

    def datasets(self):
        """The dataset names (groups of samples), sorted."""
        return sorted({e.get("dataset", "") for e in self.samples.values()} - {""})

    def dataset(self, name):
        return self.samples[name].get("dataset", "")

    def type(self, name):
        return self.samples[name]["type"]

    def binning(self, name):
        entry = self.samples[name]
        if "binning" not in entry:
            return self.default_binning
        merged = dict(self.default_binning._asdict()); merged.update(entry["binning"])
        return _binning(merged, name)

    def source(self, name):
        """(kind, full path) of a sample: kind is 'list' or 'file'."""
        entry = self.samples[name]
        kind = "list" if "list" in entry else "file"
        return kind, full_path(entry[kind], self.prod_dir)

    def files(self, name):
        """ROOT files of a sample. Raises ConfigError if the list is missing or empty."""
        kind, path = self.source(name)
        if kind == "file":
            return [path]
        try:
            with open(path) as f:
                files = [l.strip() for l in f if l.strip() and not l.lstrip().startswith("#")]
        except OSError as e:
            raise ConfigError(f"{name}: cannot read file list {path} ({e.strerror})")
        if not files:
            raise ConfigError(f"{name}: empty file list {path}")
        return files

    def select(self, patterns=(), datasets=()):
        """Samples matching any fnmatch pattern or belonging to any dataset, in file order. Raises ConfigError for a miss."""
        hit = set()
        for p in patterns:
            m = [n for n in self.samples if fnmatch.fnmatchcase(n, p)]
            if not m:
                raise ConfigError(f"no sample matches '{p}'")
            hit.update(m)
        for ds in datasets:
            m = [n for n, e in self.samples.items() if e.get("dataset") == ds]
            if not m:
                raise ConfigError(f"no sample in dataset '{ds}' (datasets: {', '.join(self.datasets())})")
            hit.update(m)
        return [n for n in self.samples if n in hit]


def load(path=DATASETS_FILE, prod_dir=None):
    """Read and validate datasets.json. prod_dir: argument, else $HSCP_PROD_DIR, else the file's value."""
    try:
        with open(path) as f:
            doc = json.load(f)
    except OSError as e:
        raise ConfigError(f"cannot read {path} ({e.strerror})")
    except json.JSONDecodeError as e:
        raise ConfigError(f"{path}: invalid JSON at line {e.lineno}: {e.msg}")
    return Config(doc, path, prod_dir)
