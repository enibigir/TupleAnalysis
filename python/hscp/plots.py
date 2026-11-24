"""
Plot configuration for python/plots.py: named input files, and the macro functions to call with them.

cfg/plots.json:
    {
      "output_dir": "outputDisplay/V1",
      "lumi_label": "7.98 fb^{-1} (13.6 TeV)",
      "inputs": {"data": {"file": "output/JetMET2024_V1/JetMET2024.root", "label": "JetMET 2024C"},
                 "qcd":  "output/QCD2024_V1/QCD2024_weighted.root",
                 "gluino_mass": "output/Gluino_V1/Gluino_Run3_MET_madgraph_%d_weighted.root"},
      "plots": [
        {"function": "nPV", "args": ["qcd", "ttbar", "wjets", "data"]},
        {"name": "YieldAfterSF_MET", "function": "YieldAfterSF", "args": [true, false, true, false, "wjets", ...]}
      ]
    }

Paths (inputs, output_dir, macro) are relative to the repository unless absolute; an input may hold %d
(a pattern the function fills with a mass), and exists if one file matches. In "args", a string is the
name of an input, true/false and numbers are passed as they are, {"str": "text"} is a literal string, and
null an empty path: the sample is left out (backgrounds, gluino mass slots; see the functions of the macro).
A plot is named after its function unless "name" is given; names must be unique. "note" is free text.

Labels: an input is a path, or {"file": path, "label": text, "color": c} to set its legend entry (ROOT
TLatex, e.g. "#tilde{g} (m=1600 GeV)") and its colour in python/stack.py (a ROOT colour number, or a name
such as "kRed+2"); "lumi_label" sets the luminosity text. A label goes to the macro under the name
of the function parameter that receives the input (Label("qcd", ...) in macros/CombineHistos.C); legends
without a configured label keep the text written in the macro.

"stack" (optional) gives the samples of python/stack.py, by input name:
    "stack": {"data": "data", "backgrounds": ["qcd", "ttbar", "wjets"], "signals": ["gluino"]}
backgrounds in stacking order (the first at the bottom); "data": null for none. "plots" may be empty.
"""

import glob
import json
import os
import re
from collections import namedtuple

from . import config

PLOTS_FILE = os.path.join(config.CFG_DIR, "plots.json")
DEFAULT_MACRO = "macros/CombineHistos.C"
DEFAULT_OUTPUT_DIR = "outputDisplay"

Plot = namedtuple("Plot", "name function args note")
Param = namedtuple("Param", "text has_default")


class PlotError(Exception):
    pass


class PlotConfig:
    def __init__(self, doc, source="plots.json", base=config.REPO):
        if not isinstance(doc, dict):
            raise PlotError(f"{source}: expected a JSON object")
        unknown = set(doc) - {"output_dir", "inputs", "plots", "macro", "lumi_label", "stack", "_comment"}
        if unknown:
            raise PlotError(f"{source}: unknown key(s) {', '.join(sorted(unknown))}")
        self.source, self.base = source, base
        self.output_dir = self.path(doc.get("output_dir", DEFAULT_OUTPUT_DIR))
        self.macro = self.path(doc.get("macro", DEFAULT_MACRO))
        inputs = doc.get("inputs", {})
        if not isinstance(inputs, dict):
            raise PlotError(f"{source}: 'inputs' must map names to file paths")
        self.inputs, self.labels, self.colors = {}, {}, {}
        for k, v in inputs.items():
            if isinstance(v, dict) and set(v) <= {"file", "label", "color"} and isinstance(v.get("file"), str) \
                    and v["file"] and isinstance(v.get("label", ""), str) \
                    and (isinstance(v.get("color", ""), str) or type(v.get("color")) is int):
                file, label, color = v["file"], v.get("label"), v.get("color")
            elif isinstance(v, str) and v:
                file, label, color = v, None, None
            else:
                raise PlotError(f'{source}: input {k} must be a path or {{"file": path, "label": text, "color": c}}')
            self.inputs[k] = self.path(file)
            if label:
                self.labels[k] = label
            if color is not None and color != "":
                self.colors[k] = color
        self.lumi_label = doc.get("lumi_label")
        if self.lumi_label is not None and not isinstance(self.lumi_label, str):
            raise PlotError(f"{source}: 'lumi_label' must be a string")
        self.stack = self._read_stack(doc.get("stack"))
        plots = doc.get("plots", [])
        if not isinstance(plots, list):
            raise PlotError(f"{source}: 'plots' must be a list")
        self.plots = []
        for i, p in enumerate(plots):
            where = f"{source}: plots[{i}]"
            if not isinstance(p, dict) or not isinstance(p.get("function"), str):
                raise PlotError(f"{where}: needs a 'function'")
            unknown = set(p) - {"name", "function", "args", "note"}
            if unknown:
                raise PlotError(f"{where}: unknown key(s) {', '.join(sorted(unknown))}")
            args = p.get("args", [])
            if not isinstance(args, list):
                raise PlotError(f"{where}: 'args' must be a list")
            for a in args:
                self._check_arg(a, where)
            self.plots.append(Plot(p.get("name", p["function"]), p["function"], args, p.get("note", "")))
        names = [p.name for p in self.plots]
        dup = sorted({n for n in names if names.count(n) > 1})
        if dup:
            raise PlotError(f"{source}: plot name(s) used twice: {', '.join(dup)} (give a 'name' to each)")

    def _read_stack(self, s):
        if s is None:
            return None
        where = f"{self.source}: stack"
        if not isinstance(s, dict) or set(s) - {"data", "backgrounds", "signals"}:
            raise PlotError(f'{where}: expected {{"data": name, "backgrounds": [names], "signals": [names]}}')
        data, bkgs, sigs = s.get("data"), s.get("backgrounds", []), s.get("signals", [])
        if data is not None and not isinstance(data, str):
            raise PlotError(f"{where}: 'data' must be an input name or null")
        for key, value in (("backgrounds", bkgs), ("signals", sigs)):
            if not isinstance(value, list) or not all(isinstance(n, str) for n in value):
                raise PlotError(f"{where}: '{key}' must be a list of input names")
        for n in ([data] if data else []) + bkgs + sigs:
            if n not in self.inputs:
                raise PlotError(f"{where}: '{n}' is not an input (inputs: {', '.join(sorted(self.inputs))})")
        return {"data": data, "backgrounds": bkgs, "signals": sigs}

    def path(self, p):
        return p if os.path.isabs(p) else os.path.normpath(os.path.join(self.base, p))

    def _check_arg(self, a, where):
        if isinstance(a, str):
            if a not in self.inputs:
                raise PlotError(f"{where}: '{a}' is not an input (inputs: {', '.join(sorted(self.inputs))}); "
                                'write {"str": "..."} for a literal string')
        elif isinstance(a, dict):
            if set(a) != {"str"} or not isinstance(a["str"], str):
                raise PlotError(f'{where}: a literal string is written {{"str": "..."}}, not {a}')
        elif a is not None and not isinstance(a, (bool, int, float)):
            raise PlotError(f"{where}: unsupported argument {a!r}")

    def select(self, patterns):
        """Plots whose name matches one of the fnmatch patterns (all plots without patterns)."""
        if not patterns:
            return list(self.plots)
        import fnmatch
        chosen = [p for p in self.plots if any(fnmatch.fnmatchcase(p.name, pat) for pat in patterns)]
        unmatched = [pat for pat in patterns if not any(fnmatch.fnmatchcase(p.name, pat) for p in self.plots)]
        if unmatched:
            raise PlotError(f"no plot named {', '.join(unmatched)} (plots: {', '.join(p.name for p in self.plots)})")
        return chosen

    def input_names(self, plot):
        return [a for a in plot.args if isinstance(a, str)]

    def missing_inputs(self, plot):
        """[(name, path)] of the inputs of plot that do not exist."""
        return [(n, self.inputs[n]) for n in dict.fromkeys(self.input_names(plot)) if not input_exists(self.inputs[n])]

    def call(self, plot):
        """The C++ call, e.g. nPV("/abs/qcd.root", "/abs/data.root")."""
        return f"{plot.function}({', '.join(self.cpp_arg(a) for a in plot.args)})"

    def legend_labels(self, plot, params):
        """{key: label} for Label(key, ...) of the macro: the parameter name receiving each labelled input,
        and "lumi". params is the signature of the function."""
        out = {}
        for a, param in zip(plot.args, params):
            if isinstance(a, str) and a in self.labels:
                out[param_name(param)] = self.labels[a]
        if self.lumi_label:
            out["lumi"] = self.lumi_label
        return out

    def cpp_arg(self, a):
        if a is None:
            return '""'
        if isinstance(a, bool):
            return "true" if a else "false"
        if isinstance(a, (int, float)):
            return repr(a)
        if isinstance(a, dict):
            return cpp_string(a["str"])
        return cpp_string(self.inputs[a])


def load(path=PLOTS_FILE):
    try:
        with open(path) as f:
            doc = json.load(f)
    except OSError as e:
        raise PlotError(f"cannot read {path}: {e.strerror}")
    except json.JSONDecodeError as e:
        raise PlotError(f"{path}: invalid JSON ({e})")
    return PlotConfig(doc, path)


def input_exists(path):
    if "%d" in path:
        return bool(glob.glob(path.replace("%d", "*")))
    return os.path.isfile(path)


def cpp_string(s):
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"'


# ---------------------------------------------------------------------- the macro
_DEF = r"^[A-Za-z][\w:<>\*& ]*[\s\*&]{}\s*\("


def signature(macro_text, function):
    """[Param] of the top-level function of macro_text, or None if it is not defined (or defined twice)."""
    starts = [m.end() for m in re.finditer(_DEF.format(re.escape(function)), macro_text, re.M)]
    if len(starts) != 1:
        return None
    params, depth, cur, quote, i = [], 0, "", None, starts[0]
    while i < len(macro_text):
        c = macro_text[i]
        if quote:
            cur += c
            if c == "\\":
                cur += macro_text[i + 1]
                i += 1
            elif c == quote:
                quote = None
        elif c in "\"'":
            quote, cur = c, cur + c
        elif c in "(<{[":
            depth, cur = depth + 1, cur + c
        elif c in ")>}]" and depth:
            depth, cur = depth - 1, cur + c
        elif c == ")":
            break
        elif c == "," and not depth:
            params.append(cur)
            cur = ""
        else:
            cur += c
        i += 1
    if cur.strip():
        params.append(cur)
    return [Param(" ".join(p.split()), "=" in p) for p in params if p.strip() and p.strip() != "void"]


def param_name(param):
    """The name of a parameter: const char *qcd = "x.root" -> qcd."""
    names = re.findall(r"\w+", param.text.split("=", 1)[0])
    return names[-1] if names else ""


def function_body(macro_text, function):
    """The text of the top-level function, from its definition to the next one ("" if not defined)."""
    m = re.search(_DEF.format(re.escape(function)), macro_text, re.M)
    if not m:
        return ""
    nxt = re.compile(_DEF.format(r"\w+"), re.M).search(macro_text, m.end())
    return macro_text[m.start():nxt.start() if nxt else len(macro_text)]


def histogram_names(macro_text, function):
    """Names read with Get("literal") in the body of function, outside // comments. Names built at run time
    (Get(Form(...)), Get(variable)) are not seen."""
    body = function_body(macro_text, function)
    code = "\n".join(line.split("//", 1)[0] for line in body.splitlines())
    return list(dict.fromkeys(re.findall(r'->\s*Get\(\s*"([^"]+)"\s*\)', code)))


def check_calls(cfg, plots, macro_text):
    """Errors for plots whose function is not in the macro, or that pass a wrong number of arguments."""
    errors = []
    for p in plots:
        params = signature(macro_text, p.function)
        if params is None:
            errors.append(f"{p.name}: no function {p.function} in {cfg.macro}")
            continue
        required = sum(not q.has_default for q in params)
        if not required <= len(p.args) <= len(params):
            want = f"{required}" if required == len(params) else f"{required} to {len(params)}"
            errors.append(f"{p.name}: {p.function} takes {want} argument(s), {len(p.args)} given:\n    "
                          + "\n    ".join(q.text for q in params))
    return errors


def root_command(cfg, plot, labels=None):
    """argv running one plot in its own root process, started in the directory of the macro.

    The call is compiled with gInterpreter->Declare, which fails (exit 3) on an unknown function or a wrong
    argument: a plain call line would be dropped by cling with exit code 0. Exit 2: the macro did not load."""
    call = cfg.call(plot)
    if ")HSCP" in call:
        raise PlotError(f"{plot.name}: the arguments contain ')HSCP'")
    set_labels = []
    for key, text in (labels or {}).items():
        set_labels += ["-e", f"gLabels[{cpp_string(key)}] = {cpp_string(text)};"]
    return ["root", "-l", "-b", "-q",
            "-e", f"if (gROOT->LoadMacro({cpp_string(os.path.basename(cfg.macro))}) != 0) gSystem->Exit(2);",
            "-e", f"gOutDir = {cpp_string(cfg.output_dir)};",
            *set_labels,
            "-e", f'if (!gInterpreter->Declare(R"HSCP(void hscp_plot() {{ {call}; }})HSCP")) gSystem->Exit(3);',
            "-e", "hscp_plot();"]
