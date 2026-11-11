# TupleAnalysis

Histogram production for the CMS Run 3 search for Heavy Stable Charged Particles (HSCP).

The code runs a ROOT `TSelector` over the HSCP ntuples (tree `HSCPMiniAODAnalyzer/Events`) and writes, for each sample, one ROOT file holding

- control, cutflow, trigger-efficiency and signal-mass histograms, for each candidate selection;
- the histograms of the regions A, B, C, D used by the data-driven background estimate.

A second set of macros turns these files into plots and into the trigger scale factors applied to the simulation.


## Disclaimer

I wrote all the code in this project myself. However, this README and the comments in the scripts were generated using Claude.

## Requirements

- ROOT 6. `source env.sh` uses the ROOT found in the `PATH`; without one, it sources the LCG view
  from CVMFS matching the machine (el8 or el9, x86_64 or aarch64). No CMSSW is needed.
- Python 3 (standard library only), to generate the selector and to drive HTCondor.
- `numpy`, `uproot` and `matplotlib`, only for `python/DeriveSF.py`.
- HTCondor (`condor_submit`), only to run on a batch farm.

## Quick start

```bash
source env.sh                # ROOT in the PATH, or the LCG view from CVMFS
./compile.sh                 # builds libTools.so from src/ and inc/
./tests/run_smoke.sh         # optional: whole chain on a synthetic ntuple, a few minutes
./launch.sh --list '*'       # the samples known to cfg/datasets.json
./launch.sh TestMET2024      # generates and compiles the selector, runs it on that sample
```

The result is `output/<sample>.root`, with the log next to it.

To run on HTCondor, see [`doc/condor_runbook.md`](doc/condor_runbook.md): every command from a
fresh checkout to the merged histograms, in order.

Before the first run, adapt the site-specific paths listed [below](#site-specific-paths).

## Layout

| Path | Content |
|---|---|
| `cfg/datasets.json` | the samples, where their ntuples are, their type and dataset, and the binning |
| `env.sh` | provides ROOT: the one in the `PATH`, else an LCG view from CVMFS |
| `cfg/HSCPpreSelections.csv` | candidate selections, one per row |
| `src/`, `inc/` | helper library `libTools` (histogram containers, mass reconstruction, scale-factor readers); see [`doc/src_inc.md`](doc/src_inc.md) |
| `inc/SF_*.txt` | the three trigger scale-factor tables read by the selector |
| `compile.sh` | builds `libTools.so` |
| `python/hscp/` | shared code of the tools: configuration files (`config.py`), selector generation (`selector.py`), running the driver (`driver.py`) |
| `launch.sh`, `python/launch.py` | run the selector on samples chosen by name, pattern or dataset |
| `python/CreateSelector.py` | generate `macros/HSCPSelector.{h,C}` by hand (the launcher and the submitter do it themselves) |
| `macros/HSCPTemplate.{h,C}` | the analysis code: templates of the selector. **Edit these**, not `HSCPSelector.{h,C}`, which are regenerated at each launch |
| `macros/macro.cc` | driver: builds the chain of the ntuples of a file list and runs the selector; called by `launch.sh` and the Condor jobs |
| `condor/` | `submit.py`, `runJob.sh`, `merge.py`: run the samples as HTCondor jobs and merge the outputs |
| `tests/` | unit tests of `python/hscp` (`test_hscp.py`) and the smoke test on a synthetic ntuple (`run_smoke.sh`) |
| `macros/CombineHistos.C` | plotting macro, reads `output/` and writes `outputDisplay/` |
| `python/DeriveSF.py` | computes the trigger scale factors and writes the three tables |
| `output/` | where the plotting macros expect the selector outputs (content not tracked) |
| `outputDisplay/` | plots (content not tracked) |

## Processing chain

```
cfg/HSCPpreSelections.csv ─┐
macros/HSCPTemplate.{h,C} ─┴─ python/hscp/selector.py ─► macros/HSCPSelector.{h,C}

cfg/datasets.json ─► launch.sh ─► macros/macro.cc ─► TChain of ntuples ─► HSCPSelector ─► output/<sample>.root
                                                                  ▲
                                             libTools.so, inc/SF_*.txt
```

For each event the selector evaluates the trigger and, for simulation, the weight (pile-up weight times trigger scale factor). For each selection it keeps the most ionising candidate, the one with the largest `Ih`, and fills the histograms and the regions with it.

## Configuration

### `cfg/datasets.json`

One JSON file, read by `launch.sh` and `condor/submit.py` only (`macro.cc` receives a file list and the binning as arguments).

| Key | Meaning |
|---|---|
| `prod_dir` | directory of the ntuples and of their file lists; overridden by `--prod-dir` or the environment variable `HSCP_PROD_DIR` |
| `binning` | default pT cut (`ptcut`, separating regions A, B from C, D) and bin counts `etabins`, `ihbins`, `pbins`, `massbins`, `fpixbins` of the region histograms |
| `samples` | one entry per sample name, see below |

A **sample** is one entry, for example `JetMET2024C`; a **dataset** is the group of samples it belongs to, for example `JetMET2024` for the eras C to I, or `Stop` for the mass points. Each sample has exactly one of `list` (a text file with one ROOT file per line) or `file` (a single ROOT file), relative to `prod_dir` unless absolute or a `root://` URL; a `type` among `data`, `mc` and `signal`; its `dataset`, used by `--dataset`; optionally `mass` (GeV, informative) and a `binning` object overriding any of the default keys.

```json
"Stop_Run3_MET_madgraph_700": {"file": "SIGNAL/V21p0/HSCP-Stop_Par-M-700_merged.root", "type": "signal", "dataset": "Stop", "mass": 700}
```

The binning drives the memory used: with 96 η, 400 Ih, 1000 p, 400 mass and 20 Fpix bins, the region histograms take about 1.8 GB **per selection**.

### `cfg/HSCPpreSelections.csv`

The first row is a header. Each following row is `label, expression`:

- `label` names the selection and prefixes its histograms;
- `expression` is a C++ boolean expression written with the branch readers declared in `macros/HSCPTemplate.h` and the candidate index `i`, for instance `(Pt[i] > 50.0) && (std::abs(Eta[i]) < 2.4)`. It must not contain a comma.

A row whose label starts with `#` is skipped.

### Samples

The samples of `cfg/datasets.json`, by dataset:

| Sample name | Ntuples |
|---|---|
| `JetMET2024C` … `JetMET2024I` | file lists `JetMET2024/V12p310.txt` … `V12p316.txt` |
| `Mu2024C` … `Mu2024I` | file lists `Mu2024/V18p100.txt` … `V18p106.txt` |
| `MuonEG2024C` … `MuonEG2024I` | file lists `MuonEG2024/V17p40.txt` … `V17p46.txt` |
| `QCD2024_mu_pt15to20` … `QCD2024_mu_pt1000` | file lists `BKG/QCD2024/V16p21.txt` … `V16p212.txt` |
| `Wjets2024_1J_pt40to100` … `Wjets2024_2J_pt600` | file lists `BKG/Wjets2024/V14p601.txt` … `V14p6010.txt` |
| `WjetMuNu2024`, `TTbar2024`, `TTbarSemiLep2024` | one file list each |
| `Gluino_Run3_MET_pythia_<mass>` | one ROOT file per mass (`CodeV19p0` production) |
| `Gluino_Run3_MET_madgraph_<mass>` | one ROOT file per mass (`V19p6` production) |
| `Gluino_Run2_MET_madgraph_2000` | one ROOT file |
| `Stau_Run3_MET_<mass>`, `Stop_Run3_MET_madgraph_<mass>` | one ROOT file per mass |
| `TestMET2024`, `TestMuon2024`, `TestMuonEG`, `TestTTbar2024`, `TestWjets`, `TestWjetsMuNu` | short file lists, for tests |

A file list is a text file with one ROOT file per line.

## Running

```bash
./launch.sh 'JetMET2024*'                        # every era of the data stream
./launch.sh --dataset QCD2024 --dataset Wjets2024 --jobs 3   # three samples at a time, one thread each
./launch.sh 'Gluino_Run3_MET_madgraph_*' --out-dir output/Gluino_V20 --label V20
./launch.sh TestMET2024 --dry-run                # print the root commands only
./launch.sh --list '*'                           # names and binning
```

Samples are chosen by name or `fnmatch` pattern, or whole datasets with `--dataset`, in `cfg/datasets.json`. The selector is generated from `cfg/HSCPpreSelections.csv` (or `--selections`) and compiled once; each sample then runs in its own `root` process, writing `<out-dir>/<sample>[_<label>].root`, the file list it used (`.files`) and a `.log` next to it. The script stops at the first failure unless `--keep-going` is given, and exits non-zero if any sample failed. `launch.sh` sources `env.sh` and calls `python/launch.py`, which can be used directly once ROOT is in the `PATH`.

`macro.cc` can also be called directly, which is what the launcher and the Condor jobs do:

```bash
cd macros
root -l -b -q 'macro.cc("JetMET2024C","JetMET2024C.files","JetMET2024C.root",true,70,96,400,1000,400,20)'
```

sample name, file list, output file, whether to write the no-selection histograms, then the pT cut and the five binnings; `./launch.sh <sample> --dry-run` prints the exact command. `HSCP_NTHREADS` sets the number of implicit-MT threads (default 4).

## Output file

`<out-dir>/<sample>[_<label>].root` contains flat histograms, without directories. Use `--out-dir` and `--label` to keep several runs of the same sample apart; the plotting macros read `output/<family>_V<n>/`.

| Name | Content |
|---|---|
| `<selection>_<quantity>` | histograms of the candidates passing one selection, e.g. `<selection>_Ih`, `<selection>_9fp10_SignalMass_nominal` |
| names without a selection prefix | event-level and pre-selection histograms, e.g. `EventCutflow`, `CandidateCutflow`, `Nm1_<cut>`, `Nosel_<quantity>` |
| `<quantity>_region<A\|B\|C\|D>_<lo>fp<hi>_<selection>` | region histograms, e.g. `ih_eta_regionD_8fp9_<selection>` |

Regions A and B are filled by candidates with pT below the cut, regions C and D by candidates above it. `<lo>fp<hi>` is the interval of the pixel discriminator: `3fp8` means 0.3 < Fpix ≤ 0.8, `99fp10` means 0.99 < Fpix ≤ 1.

For data, several histograms are only filled by candidates with Fpix ≤ 0.9, to keep the signal region blind.

The mass is computed from the momentum and `Ih` as `m = p · sqrt((Ih − C) / K)`. `K` and `C` are chosen from the sample name at the beginning of `HSCPSelector::SlaveBegin`: one pair for 2024 data, one for simulation.

The list and binning of the region histograms are documented in [`doc/src_inc.md`](doc/src_inc.md).

## Plots and scale factors

- **`macros/CombineHistos.C`**, run from `macros/` with `root -l -b -q CombineHistos.C`. Its main function, at the end of the file, is a list of calls: enable a plot by uncommenting its call. Input files are read from `../output/<sample>_V<n>/`. Plots are written to `../outputDisplay/` and to its sub-directories `TriggEff/` and `Nm1plots/`, which must exist beforehand; a few functions take their output directory as an argument.
- **`python/DeriveSF.py`** computes the trigger scale factors, in bins of pseudo MET, of PUppi MET, and of both. Its paths are relative to the repository, so it can be started from anywhere: `python3 python/DeriveSF.py`.

Merging the selector outputs and weighting the simulated samples by cross-section (the `*_weighted.root` files read by `CombineHistos.C`) is done outside this repository.

### Trigger scale-factor tables

The selector reads three tables in `inc/` at start-up and stops if one is missing. All three are written by `python/DeriveSF.py`:

| File in `inc/` | Content |
|---|---|
| `SF_PseudoMET.txt` | scale factors in bins of pseudo MET |
| `SF_PUppiMET.txt` | scale factors in bins of PUppi MET |
| `SF_orMETtrg_PUppiMET_VS_PseudoMET__TriggerEffCalib_table_plain.txt` | scale factors in bins of (pseudo MET, PUppi MET), on a non-grid binning |

The script reads the data and simulation files named at the top of its `CONFIG` block, produced with a selection labelled `TriggerEffCalib`. It also writes the LaTeX version of the tables and the maps of the two-dimensional scale factors in `outputDisplay/TriggEff/`.

The format of the tables is described in [`doc/src_inc.md`](doc/src_inc.md).

## Site-specific paths

| Where | What |
|---|---|
| `cfg/datasets.json`, `prod_dir` | directory of the ntuples and of their file lists (overridden by the environment variable `HSCP_PROD_DIR` or `--prod-dir`) |
| `env.sh` | LCG release and platform of the CVMFS view used when no ROOT is installed |
| `macros/CombineHistos.C` | five absolute paths (search for `/safe/`) |

The selector reads the scale-factor tables from `HSCP_SF_DIR` (default `inc/`).

## Extending the analysis

- **New sample**: add an entry in `cfg/datasets.json`. The selector decides from the sample name which trigger to require, whether the sample is simulation or signal, and which `K` and `C` to use: check the `sample_.find(...)` conditions in `HSCPTemplate.C`.
- **New selection**: add a row to `cfg/HSCPpreSelections.csv`. No code change is needed.
- **New histogram**: book it with `AddHisto1F` / `AddHisto2F` in `SlaveBegin` and fill it with `FillHisto1F` / `FillHisto2F` in `Process`, both in `HSCPTemplate.C`. The two names must match exactly: a fill with an unknown name is silently ignored.
- **New branch**: declare a `TTreeReaderValue` or `TTreeReaderArray` in `HSCPTemplate.h`, then run `python3 tests/make_synthetic_source.py` so the synthetic ntuple of the smoke test gets the branch too.
- Any change in `src/` or `inc/` requires `./compile.sh`.