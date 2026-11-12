"""
Generation of macros/HSCPSelector.{h,C} from the templates and the selection file.

The templates are copied line by line. After each marker line, the code generated
from the selections is inserted:
    ADD-HSCP-SELECTION     (HSCPTemplate.h)  declarations of the PassHSCPpresel_<label> methods
    ADD-SELECTION-METHODS  (HSCPTemplate.C)  their implementations
    FILL-SELECTION-VECTOR  (HSCPTemplate.C)  registration of each method and of its label

A selection file is JSON (cfg/selections.json):

    {
      "common": ["Pt[i] > 50.0", "Ih_Strip[i] > C"],
      "selections": {
        "Eta2p4": ["std::abs(Eta[i]) < 2.4"],
        "Eta1":   {"cuts": ["std::abs(Eta[i]) < 1"], "enabled": false, "comment": "kept for reference"}
      }
    }

A cut is a C++ boolean expression written with the TTreeReader members of HSCPSelector and
the candidate index i. A selection is the AND of the common cuts and of its own cuts, each in
parentheses. A selection is either its list of cuts or an object with "cuts" and the optional
"enabled" (default true) and "comment". Keys starting with _ are ignored at both levels.
"""

import json
import os
import re
from collections import namedtuple

from .config import MACROS_DIR, SELECTIONS_FILE

IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")
SELECTION_KEYS = {"cuts", "enabled", "comment"}

Selection = namedtuple("Selection", "label cuts")


class SelectionError(ValueError):
    pass


def expression(cuts):
    """The C++ expression of a list of cuts."""
    return " && ".join(f"({c})" for c in cuts)


def _cut_list(value, where, allow_empty=False):
    if not isinstance(value, list) or not all(isinstance(c, str) and c.strip() for c in value):
        raise SelectionError(f"{where}: expected a list of cut strings")
    if not value and not allow_empty:
        raise SelectionError(f"{where}: no cut")
    return [c.strip() for c in value]


def read_selections(path=SELECTIONS_FILE):
    """Active selections of a selection file, each with its full list of cuts (common ones first)."""
    try:
        with open(path) as f:
            doc = json.load(f)
    except json.JSONDecodeError as e:
        raise SelectionError(f"{path}: not valid JSON ({e})")
    if not isinstance(doc, dict) or not isinstance(doc.get("selections"), dict):
        raise SelectionError(f"{path}: expected an object with a 'selections' object")
    unknown = {k for k in doc if not k.startswith("_")} - {"common", "selections"}
    if unknown:
        raise SelectionError(f"{path}: unknown key(s) {sorted(unknown)}")
    common = _cut_list(doc.get("common", []), f"{path}: common", allow_empty=True)

    sels = []
    for label, value in doc["selections"].items():
        if label.startswith("_"):
            continue
        where = f"{path}: selection '{label}'"
        if not IDENT.match(label):
            raise SelectionError(f"{where}: label is not a valid C++ identifier")
        if isinstance(value, dict):
            unknown = {k for k in value if not k.startswith("_")} - SELECTION_KEYS
            if unknown:
                raise SelectionError(f"{where}: unknown key(s) {sorted(unknown)}")
            if not value.get("enabled", True):
                continue
            if "cuts" not in value:
                raise SelectionError(f"{where}: 'cuts' missing")
            value = value["cuts"]
        sels.append(Selection(label, common + _cut_list(value, where)))
    if not sels:
        raise SelectionError(f"{path}: no active selection")
    return sels


def write_selections(sels, path):
    """Write a selection file (the inverse of read_selections), for one group of selections."""
    doc = {"common": [], "selections": {label: list(cuts) for label, cuts in sels}}
    with open(path, "w") as f:
        json.dump(doc, f, indent=2)
        f.write("\n")


# ---- code generators -------------------------------------------------------------
def _proto(label):
    return f"bool PassHSCPpresel_{label}(int hscpIndex);\n"


def _impl(label, cuts):
    code = f"bool HSCPSelector::PassHSCPpresel_{label}(int i){{\n"
    if label == "CalibPseudoMET":          # event-level selection: no candidate index check
        code += "   if (i < 0) {\n"
    else:
        code += "   if (i<0 || i>=(int)Pt.GetSize()) {\n"
    code += "      cout << i << endl;\n      return false;\n   }\n"
    code += f"   return {expression(cuts)};\n}}\n"
    return code


def _register(label):
    return (f'selections_.push_back(&HSCPSelector::PassHSCPpresel_{label});\n'
            f'selLabels_.push_back("{label}");\n')


def code_prototypes(sels):
    return "".join(_proto(l) + "\n" for l, _ in sels)


def code_implementations(sels):
    return "".join(_impl(l, c) + "\n" for l, c in sels)


def code_registration(sels):
    return "".join(_register(l) + "\n" for l, _ in sels)


def _copy_with_inserts(template, output, markers):
    with open(template) as src, open(output, "w") as dst:
        for line in src:
            dst.write(line)
            for marker, content in markers.items():
                if line.find(marker) > 0:
                    dst.write(content)


def generate(sels, out_dir=MACROS_DIR, templates_dir=MACROS_DIR):
    """Write HSCPSelector.h and HSCPSelector.C into out_dir. Returns their paths."""
    os.makedirs(out_dir, exist_ok=True)
    h = os.path.join(out_dir, "HSCPSelector.h")
    c = os.path.join(out_dir, "HSCPSelector.C")
    _copy_with_inserts(os.path.join(templates_dir, "HSCPTemplate.h"), h,
                       {"ADD-HSCP-SELECTION": code_prototypes(sels)})
    _copy_with_inserts(os.path.join(templates_dir, "HSCPTemplate.C"), c,
                       {"ADD-SELECTION-METHODS": code_implementations(sels),
                        "FILL-SELECTION-VECTOR": code_registration(sels)})
    return h, c
