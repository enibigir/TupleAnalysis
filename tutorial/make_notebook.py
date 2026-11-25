#!/usr/bin/env python3
"""
Write tutorial/HSCP_tutorial.ipynb from the cells below (edit them here, then: python3 tutorial/make_notebook.py).
"""

import json
import os

CELLS = []


def md(text):
    CELLS.append(("markdown", text.strip("\n")))


def code(text):
    CELLS.append(("code", text.strip("\n")))


# ======================================================================
md(r"""
# HSCP search with TupleAnalysis: from ntuples to plots

This notebook runs the analysis chain of the repository on a **few files of a few samples**, on SWAN:

| Step | Tool | Output |
|---|---|---|
| 1. Select events, fill histograms | `launch.sh` (the `HSCPSelector` of `macros/`) | one ROOT file per sample, `output/tutorial/` |
| 2. Normalise the simulation | `python/Weight.py` | weighted files, `output/tutorial/weighted/` |
| 3. Data / simulation plots | `python/stack.py` | stacked histograms, `outputDisplay/tutorial/stack/` |
| 4. Analysis plots | `python/plots.py` (functions of `macros/CombineHistos.C`) | `outputDisplay/tutorial/` |

On the full samples, step 1 runs on HTCondor (`condor/submit.py`) and its outputs are added with
`condor/merge.py`: that part is described in section 6, not run here.

**Before you start**

* SWAN session: a recent LCG software stack, **4 cores and 16 GB** of memory.
* A grid certificate registered in the CMS VO: the ntuples are read with `davs://` from the Strasbourg
  storage, which needs a VOMS proxy (section 1).
* About 20 minutes for the whole notebook.
""")

md(r"""
## 0. Go to the repository

Open this notebook from `TupleAnalysis/tutorial/`. The cell below moves to the repository, where every
command of this notebook is run.
""")

code(r"""
import os, glob, json, subprocess, time

NOTEBOOK_DIR = os.getcwd()
if os.path.basename(NOTEBOOK_DIR) == "tutorial":
    os.chdir("..")
REPO = os.getcwd()
assert os.path.isfile("launch.sh"), f"not in the repository: {REPO}"
print("repository:", REPO)
""")

# ======================================================================
md(r"""
## 1. Grid proxy

The input files are on `davs://sbgdcache.in2p3.fr`, readable with a CMS VOMS proxy. A proxy is made
with your password, so **not in this notebook**: in a terminal.

**Option A: on LXPLUS, then copy it to CERNBox** (SWAN sees your CERNBox under `/eos/user/...`):
```bash
ssh lxplus.cern.ch
voms-proxy-init --voms cms --valid 168:00 --out /eos/user/${USER:0:1}/$USER/x509up
```

**Option B: in a SWAN terminal** (*New* → *Terminal*), with your certificate in `~/.globus/`:
```bash
export X509_VOMS_DIR=/cvmfs/grid.cern.ch/etc/grid-security/vomsdir
export VOMS_USERCONF=/cvmfs/grid.cern.ch/etc/grid-security/vomses
voms-proxy-init --voms cms --valid 168:00 --out $HOME/x509up
```

Then set `PROXY` below to the file you wrote. The cell copies it to `/tmp` with the permissions the
grid tools require, and points ROOT to it.
""")

code(r"""
PROXY = os.path.expanduser("~/x509up")            # or "/eos/user/<l>/<you>/x509up" (option A)

local = f"/tmp/x509up_u{os.getuid()}"
if os.path.realpath(PROXY) != os.path.realpath(local):      # copying a file onto itself would empty it
    with open(PROXY, "rb") as src:
        proxy = src.read()
    with open(local, "wb") as dst:
        dst.write(proxy)
os.chmod(local, 0o600)
os.environ["X509_USER_PROXY"] = local
for d in ("/cvmfs/grid.cern.ch/etc/grid-security/certificates", "/etc/grid-security/certificates"):
    if os.path.isdir(d):
        os.environ.setdefault("X509_CERT_DIR", d)
        break
print("X509_USER_PROXY =", os.environ["X509_USER_PROXY"])
print("X509_CERT_DIR   =", os.environ.get("X509_CERT_DIR"))
""")

md(r"""
Check that one input file can be read. If it fails: is the proxy still valid (`voms-proxy-info` in a
terminal)? An error *HTTP 429* means the storage is busy: wait a minute and run the cell again.
""")

code(r"""
import ROOT

first = open("input/TestJetMET2024C.txt").readline().strip()
f = ROOT.TFile.Open(first)
assert f and not f.IsZombie(), "cannot read the input: check the proxy"
print("completed  read test of", first.split("/")[-1], "-", [k.GetName() for k in f.GetListOfKeys()][:5])
f.Close()
""")

# ======================================================================
md(r"""
## 2. The samples and the selection

Samples are described in `cfg/datasets.json`: a file list in `input/`, a type (`data`, `mc`, `signal`),
a dataset, and for simulation a cross-section in pb. We take a few of them:
""")

code(r"""
SAMPLES = {
    "data":   ["JetMET2024C"],
    "qcd":    ["QCD2024_mu_pt170to300", "QCD2024_mu_pt300to470"],
    "ttbar":  ["TTbar2024"],
    "wjets":  ["WjetMuNu2024"],
    "signal": ["Gluino_Run3_MET_madgraph_1400", "Gluino_Run3_MET_madgraph_1600", "Gluino_Run3_MET_madgraph_2000"],
}
ALL = [s for names in SAMPLES.values() for s in names]
MAX_FILES = 2                       # files read per sample

datasets = json.load(open("cfg/datasets.json"))
for s in ALL:
    d = datasets["samples"][s]
    n = sum(1 for line in open(d["list"]) if line.strip())
    print(f"{s:<32} {d['type']:<7} xsec {str(d.get('xsec', '-')):>9} pb   {n:5d} files, {MAX_FILES} read")
""")

md(r"""
The selection is the AND of C++ cuts on the candidates (`[i]`) and the event (`[0]`). The full analysis
has twelve selections in `cfg/selections.json`; the tutorial keeps two in `tutorial/selections.json`,
so that a sample needs about 1.5 GB of memory instead of 8 GB:
""")

code(r"""
sel = json.load(open("tutorial/selections.json"))
print("cuts common to every selection:")
for c in sel["common"]:
    print("   ", c)
for name, cuts in sel["selections"].items():
    print(f"\n{name}: common cuts and", cuts)
""")

# ======================================================================
md(r"""
## 3. Run the selector: `launch.sh`

`compile.sh` builds the helper library `libTools.so`; `launch.sh` generates the selector from the
selection file, compiles it once, then runs every sample in its own `root` process. Two at a time
(`--jobs 2`) fits the memory of the session. Each sample writes `output/tutorial/<sample>.root`, the list
of files it read (`.files`) and a log (`.log`).

This is the long step: a few minutes per sample.
""")

code(r"""
!./compile.sh
""")

code(r"""
t0 = time.time()
!./launch.sh {" ".join(ALL)} --max-files {MAX_FILES} --selections tutorial/selections.json --out-dir output/tutorial --jobs 2 --keep-going
print(f"{(time.time() - t0) / 60:.1f} min")
""")

md(r"""
A sample marked `FAILED` has its reason in its log, e.g. `output/tutorial/TTbar2024.log`. An *HTTP 429*
error means the storage refused too many requests: run the cell again with only the failed samples
(replace `{" ".join(ALL)}` by their names).

Look into one output: it holds flat histograms, named `<selection>_<variable>`, and the no-selection ones
(`Nosel_...`). Bin 1 of `EventCutflow` is the sum of the event weights read: it normalises the simulation
in the next step.
""")

code(r"""
f = ROOT.TFile.Open("output/tutorial/Gluino_Run3_MET_madgraph_1600.root")
names = [k.GetName() for k in f.GetListOfKeys()]
print(len(names), "histograms, e.g.", [n for n in names if n.endswith("_Fpix")])
print("events read:", f.Get("EventCutflow").GetBinContent(1))

c = ROOT.TCanvas("c", "", 700, 450)
h = f.Get("METanalysis_TestPUppiMETCut_Eta2p4_Fpix")
h.SetLineWidth(2)
h.Draw("HIST")
c.Draw()
""")

# ======================================================================
md(r"""
## 4. Normalise the simulation: `Weight.py`

Every histogram of a simulated sample is scaled by

$$ w = \frac{\sigma \times L}{\sum \text{weights}} $$

with $\sigma$ from `cfg/datasets.json`, $L$ the luminosity of the data, and the sum of weights read
from `EventCutflow`. Several samples (the QCD $p_T$ bins) are added after scaling.

**Which luminosity?** We read only 2 of the files of `JetMET2024C`, so only a fraction of its luminosity.
Assuming the files hold the same number of events, $L \approx L_\text{list} \times 2 / N_\text{files}$:
""")

code(r"""
L_LIST = 7980.0                     # pb^-1, luminosity of the whole JetMET2024C list (7.98 fb^-1)
n_list = sum(1 for line in open(datasets["samples"]["JetMET2024C"]["list"]) if line.strip())
LUMI = L_LIST * MAX_FILES / n_list
print(f"luminosity of the {MAX_FILES} data files: {LUMI:.2f} pb^-1 ({MAX_FILES} of {n_list} files)")
""")

md(r"""
First a dry run: it prints the factor of each sample, and writes nothing. Check that they look sensible.
""")

code(r"""
!python3 python/Weight.py --samples {",".join(SAMPLES["qcd"])} --input-dir output/tutorial --lumi {LUMI} --dry-run
""")

md(r"""
Now the weighted files: the QCD bins added into one, one file per other background and per gluino mass.
The data are not scaled (factor 1): a single data sample is used as it is.
""")

code(r"""
W = "output/tutorial/weighted"
jobs = {"QCD2024_weighted": SAMPLES["qcd"], "TTbar2024_weighted": SAMPLES["ttbar"],
        "WjetMuNu2024_weighted": SAMPLES["wjets"]}
for s in SAMPLES["signal"]:
    jobs[f"Gluino_{s.rsplit('_', 1)[1]}_weighted"] = [s]

os.makedirs(W, exist_ok=True)
for out, samples in jobs.items():
    r = subprocess.run(["python3", "python/Weight.py", "--samples", ",".join(samples), "--input-dir", "output/tutorial",
                        "--lumi", str(LUMI), "-o", f"{W}/{out}.root", "-f"], capture_output=True, text=True)
    if r.returncode == 0:
        print(next(l for l in r.stdout.splitlines() if l.startswith("completed")))
    else:
        print(f"FAILED     {out:<24} {r.stderr.strip().splitlines()[-1] if r.stderr.strip() else 'see the output of Weight.py'}")
""")

# ======================================================================
md(r"""
## 5. Plots

Both plotting tools read `tutorial/plots.json`: the input files with their legend labels and colours,
which samples are data, backgrounds and signals (`stack`), and which analysis plots to make (`plots`).
We set the luminosity text and write the copy used below:
""")

code(r"""
cfg = json.load(open("tutorial/plots.json"))
cfg["lumi_label"] = f"{LUMI:.1f} pb^{{-1}} (13.6 TeV)"
PLOTS = "output/tutorial/plots.json"
json.dump(cfg, open(PLOTS, "w"), indent=1)
print(json.dumps(cfg["inputs"], indent=1)[:600], "...")
print("stack:", cfg["stack"])
""")

md(r"""
### 5a. Stacked histograms: `stack.py`

Data (points) on the stacked backgrounds, the signals as lines, and the ratio data / simulation. Give
histogram names or patterns; `--list` shows the names.
""")

code(r"""
!python3 python/stack.py --config {PLOTS} --list 'METanalysis_TestPUppiMETCut_Eta2p4_*' | head -20
""")

code(r"""
from IPython.display import Image, display

HISTS = ["EventCutflow", "METanalysis_TestPUppiMETCut_Eta2p4_Fpix", "METanalysis_TestPUppiMETCut_Eta2p4_PseudoCaloMET"]
!python3 python/stack.py --config {PLOTS} {" ".join(HISTS)} --logy --signal-scale 100 --format png,pdf

for h in HISTS:
    display(Image(f"outputDisplay/tutorial/stack/{h}_log.png", width=600))
""")

md(r"""
With 2 files per sample, the simulation has few events after the selection: expect large uncertainties,
and empty bins in the ratio.

### 5b. Analysis plots: `plots.py`

`plots.py` calls functions of `macros/CombineHistos.C` with the files of the configuration: event
cutflow, $F_\text{pixel}$, the $p_T$–$F_\text{pixel}$ plane of the background, kinematics. It writes PDF
files.
""")

code(r"""
!python3 python/plots.py --config {PLOTS} --list | head -12
!python3 python/plots.py --config {PLOTS} --jobs 2
""")

code(r"""
from IPython.display import IFrame

pdfs = sorted(glob.glob("outputDisplay/tutorial/*.pdf") + glob.glob("outputDisplay/tutorial/Nm1plots/*.pdf"))
print("\n".join(pdfs))
display(IFrame(os.path.relpath("outputDisplay/tutorial/Fpixel_plot.pdf", NOTEBOOK_DIR), width=700, height=500))
""")

# ======================================================================
md(r"""
## 6. On the full samples: HTCondor and `merge.py`

With all the files, step 3 takes hours: it runs on HTCondor, on LXPLUS or the FNAL LPC, split into
jobs of a few files and a few selections each:

```bash
python3 condor/submit.py --dataset QCD2024 --tag qcd_v1 --files-per-job 20 --selections-per-job 2 \
        --out-dir 'root://eosuser.cern.ch//eos/user/<l>/<you>/hscp/' --submit
condor_q                                             # follow the jobs
python3 condor/merge.py condor/work/qcd_v1 --check   # every job done?
python3 condor/merge.py condor/work/qcd_v1           # hadd: one file per sample, in condor/work/qcd_v1/merged/
python3 python/Weight.py --dataset QCD2024 --input-dir condor/work/qcd_v1/merged -o output/QCD2024_weighted.root
```

`merge.py` replaces step 3's single output per sample; the rest of this notebook is unchanged. See
`doc/condor_runbook.md`.

## 7. Exercises

1. **Another mass.** Add `Gluino_Run3_MET_madgraph_2600` to `SAMPLES["signal"]`, run sections 3 and 4 for
   it only, add it to `tutorial/plots.json` and to the `stack` signals.
2. **Statistics.** Rerun with `MAX_FILES = 5` for the backgrounds. How do the ratio panel and the
   uncertainty bands change? What happens to `LUMI`?
3. **A cut.** In `tutorial/selections.json`, tighten the $F_\text{pixel}$ or $I_h$ cut of one selection,
   rerun, and compare the signal and background yields printed by `stack.py`.
4. **Other variables.** Use `stack.py --list` to find the $I_h$, $p$ and mass histograms, and draw them
   with `--rebin` and `--xmin/--xmax`.
5. **Signal over background.** From the weighted files, compute $S/\sqrt{B}$ for
   `METanalysis_TestPUppiMETCut_Eta2p4_Fpix > 0.9`, for each gluino mass.
""")


# ----------------------------------------------------------------------
def notebook():
    cells = []
    for i, (kind, text) in enumerate(CELLS):
        lines = text.split("\n")
        source = [line + "\n" for line in lines[:-1]] + [lines[-1]]
        cell = {"cell_type": kind, "id": f"cell-{i:02d}", "metadata": {}, "source": source}
        if kind == "code":
            cell.update(execution_count=None, outputs=[])
        cells.append(cell)
    return {"cells": cells,
            "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
                         "language_info": {"name": "python", "pygments_lexer": "ipython3"}},
            "nbformat": 4, "nbformat_minor": 5}


if __name__ == "__main__":
    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "HSCP_tutorial.ipynb")
    with open(out, "w") as f:
        json.dump(notebook(), f, indent=1)
        f.write("\n")
    print(f"{out}: {len(CELLS)} cells")
