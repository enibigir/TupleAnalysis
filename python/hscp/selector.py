"""
Generation of macros/HSCPSelector.{h,C} from the templates and the selection file.

The templates are copied line by line. After each marker line, the code generated
from the selections is inserted:
    ADD-HSCP-SELECTION     (HSCPTemplate.h)  declarations of the PassHSCPpresel_<label> methods
    ADD-SELECTION-METHODS  (HSCPTemplate.C)  their implementations
    FILL-SELECTION-VECTOR  (HSCPTemplate.C)  registration of each method and of its label

A selection file is a CSV: a header row, then "label, C++ expression" rows written with the
TTreeReader members of HSCPSelector and the candidate index i. A row whose label starts with #
is skipped. An expression must not contain a comma.
"""

import csv
import os
import re

from .config import MACROS_DIR, SELECTIONS_FILE

IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


class SelectionError(ValueError):
    pass


def read_selections(path=SELECTIONS_FILE):
    """Active (label, expression) pairs of a selection file."""
    with open(path, newline="") as f:
        rows = [row for row in csv.reader(f, skipinitialspace=True) if row]
    sels = []
    for n, row in enumerate(rows[1:], start=2):
        label = row[0].strip()
        if label.startswith("#"):
            continue
        if len(row) != 2:
            raise SelectionError(f"{path}:{n}: expected 'label, expression' (no comma in the expression), got {len(row)} fields")
        if not IDENT.match(label):
            raise SelectionError(f"{path}:{n}: label '{label}' is not a valid C++ identifier")
        sels.append((label, row[1].strip()))
    if not sels:
        raise SelectionError(f"{path}: no active selection")
    return sels


def write_selections(sels, path):
    """Write a selection file (the inverse of read_selections), for one group of selections."""
    with open(path, "w") as f:
        f.write("selections, true\n")
        for label, expr in sels:
            f.write(f"{label}, {expr}\n")


# ---- code generators -------------------------------------------------------------
def _proto(label):
    return f"bool PassHSCPpresel_{label}(int hscpIndex);\n"


def _impl(label, expression):
    code = f"bool HSCPSelector::PassHSCPpresel_{label}(int i){{\n"
    if label == "CalibPseudoMET":          # event-level selection: no candidate index check
        code += "   if (i < 0) {\n"
    else:
        code += "   if (i<0 || i>=(int)Pt.GetSize()) {\n"
    code += "      cout << i << endl;\n      return false;\n   }\n"
    code += f"   return {expression};\n}}\n"
    return code


def _register(label):
    return (f'selections_.push_back(&HSCPSelector::PassHSCPpresel_{label});\n'
            f'selLabels_.push_back("{label}");\n')


def code_prototypes(sels):
    return "".join(_proto(l) + "\n" for l, _ in sels)


def code_implementations(sels):
    return "".join(_impl(l, e) + "\n" for l, e in sels)


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
