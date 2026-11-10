# TupleAnalysis

Histogram production for the CMS Run 3 search for Heavy Stable Charged Particles (HSCP).

The code runs a ROOT `TSelector` over the HSCP ntuples (tree `HSCPMiniAODAnalyzer/Events`) and writes, for each dataset, one ROOT file holding

- control, cutflow, trigger-efficiency and signal-mass histograms, for each candidate selection;
- the histograms of the regions A, B, C, D used by the data-driven background estimate.

A second set of macros turns these files into plots and into the trigger scale factors applied to the simulation.


## Disclaimer

I wrote all the code in this project myself. However, this README and the comments in the scripts were generated using Claude.

## Requirements

- ROOT 6, with `root-config` in the `PATH`. The code was developed inside `CMSSW_15_0_13_patch1`.
- Python 3 (standard library only), to generate the selector.
- `numpy`, `uproot` and `matplotlib`, only for `python/DeriveSF.py`.

## Quick start

```bash
./compile.sh                 # builds libTools.so from src/ and inc/
# activate one dataset: remove the leading # of one line of cfg/configFile.txt
./launchGeneral.sh           # generates the selector and runs it on that dataset
```

The result is `macros/<dataset>_<version>.root`.

Before the first run, adapt the site-specific paths listed [below](#site-specific-paths).

## Layout

| Path | Content |
|---|---|
| `cfg/configFile.txt` | dataset to process, binning and pT cut |
| `cfg/HSCPpreSelections.csv` | candidate selections, one per row |
| `src/`, `inc/` | helper library `libTools` (histogram containers, mass reconstruction, scale-factor readers); see [`doc/src_inc.md`](doc/src_inc.md) |
| `inc/SF_*.txt` | the three trigger scale-factor tables read by the selector |
| `compile.sh` | builds `libTools.so` |
| `python/CreateSelector.py`, `python/Functions.py` | generate `macros/HSCPSelector.{h,C}` from the templates and the selection file |
| `macros/HSCPTemplate.{h,C}` | the analysis code: templates of the selector. **Edit these**, not `HSCPSelector.{h,C}`, which are regenerated at each launch |
| `macros/macro.cc` | driver: reads the configuration, builds the chain of ntuples, runs the selector |
| `launch*.sh` | run one dataset or a whole family of datasets |
| `macros/CombineHistos.C` | plotting macro, reads `output/` and writes `outputDisplay/` |
| `python/DeriveSF.py` | computes the trigger scale factors and writes the three tables |
| `output/` | where the plotting macros expect the selector outputs (content not tracked) |
| `outputDisplay/` | plots (content not tracked) |

## Processing chain

```
cfg/HSCPpreSelections.csv ─┐
macros/HSCPTemplate.{h,C} ─┴─ python/CreateSelector.py ─► macros/HSCPSelector.{h,C}

cfg/configFile.txt ─► macros/macro.cc ─► TChain of ntuples ─► HSCPSelector ─► macros/<dataset>_<version>.root
                                                                  ▲
                                             libTools.so, inc/SF_*.txt
```

For each event the selector evaluates the trigger and, for simulation, the weight (pile-up weight times trigger scale factor). For each selection it keeps the most ionising candidate, the one with the largest `Ih`, and fills the histograms and the regions with it.

## Configuration

### `cfg/configFile.txt`

One line per dataset; lines starting with `#` are ignored. Exactly one line should be active (if several are, the last one is used).

| Column | Meaning |
|---|---|
| pT cut | pT threshold in GeV separating regions A, B (below) from C, D (above) |
| eta bins, ih bins, p bins, mass bins, FPIX bins | binning of the region histograms |
| type | dataset name, as known to `macros/macro.cc` |
| version | label of the run, used in the output file name. For the gluino samples it also selects the production |

The binning drives the memory used: with 96 η, 400 Ih, 1000 p, 400 mass and 20 Fpix bins, the region histograms take about 1.8 GB **per selection**.

### `cfg/HSCPpreSelections.csv`

The first row is a header. Each following row is `label, expression`:

- `label` names the selection and prefixes its histograms;
- `expression` is a C++ boolean expression written with the branch readers declared in `macros/HSCPTemplate.h` and the candidate index `i`, for instance `(Pt[i] > 50.0) && (std::abs(Eta[i]) < 2.4)`. It must not contain a comma.

A row whose label starts with `#` is skipped.

### Datasets

`knownInputs()` in `macros/macro.cc` is the list of valid dataset names and of where their ntuples are. Paths are relative to the production directory.

| Dataset name | Ntuples |
|---|---|
| `JetMET2024C` … `JetMET2024I` | file lists `JetMET2024/V12p310.txt` … `V12p316.txt` |
| `Mu2024C` … `Mu2024I` | file lists `Mu2024/V18p100.txt` … `V18p106.txt` |
| `MuonEG2024C` … `MuonEG2024I` | file lists `MuonEG2024/V17p40.txt` … `V17p46.txt` |
| `QCD2024_mu_pt15to20` … `QCD2024_mu_pt1000` | file lists `BKG/QCD2024/V16p21.txt` … `V16p212.txt` |
| `Wjets2024_1J_pt40to100` … `Wjets2024_2J_pt600` | file lists `BKG/Wjets2024/V14p601.txt` … `V14p6010.txt` |
| `WjetMuNu2024`, `TTbar2024`, `TTbarSemiLep2024` | one file list each |
| `Gluino_Run3_MET_pythia_<mass>` | one ROOT file per mass, version `V19p0` only |
| `Gluino_Run3_MET_madgraph_<mass>` | one ROOT file per mass, versions `V19p6` to `V19p12` |
| `Gluino_Run2_MET_madgraph_2000` | one ROOT file |
| `Stau_Run3_MET_<mass>`, `Stop_Run3_MET_madgraph_<mass>` | one ROOT file per mass |
| `TestMET2024`, `TestMuon2024`, `TestMuonEG`, `TestTTbar2024`, `TestWjets`, `TestWjetsMuNu` | short file lists, for tests |

A file list is a text file with one ROOT file per line.

## Running

| Script | What it runs |
|---|---|
| `launchGeneral.sh` | the active line of `cfg/configFile.txt` |
| `launchJetMET2024.sh`, `launchMu2024.sh`, `launchMuonEG.sh` | every era of the data stream, one after the other |
| `launchQCD.sh`, `launchWjets.sh` | every bin of the simulated sample |
| `launchGluinoMadgraph.sh`, `launchGluinoPythia.sh`, `launchStauRun3.sh`, `launchStopRun3.sh` | every mass point of the signal |

All scripts regenerate the selector first, and must be started from the root of the repository.

The family scripts rewrite `cfg/configFile.txt` in place: they look for the **commented** lines of their family, activate them one at a time and comment everything else. Leave the dataset lines commented before using them.

`TTbar2024`, `TTbarSemiLep2024`, `WjetMuNu2024`, `Gluino_Run2_MET_madgraph_2000` and the test samples have no family script: activate their line and use `launchGeneral.sh`.

## Output file

`macros/<dataset>_<version>.root` contains flat histograms, without directories.

| Name | Content |
|---|---|
| `<selection>_<quantity>` | histograms of the candidates passing one selection, e.g. `<selection>_Ih`, `<selection>_9fp10_SignalMass_nominal` |
| names without a selection prefix | event-level and pre-selection histograms, e.g. `EventCutflow`, `CandidateCutflow`, `Nm1_<cut>`, `Nosel_<quantity>` |
| `<quantity>_region<A\|B\|C\|D>_<lo>fp<hi>_<selection>` | region histograms, e.g. `ih_eta_regionD_8fp9_<selection>` |

Regions A and B are filled by candidates with pT below the cut, regions C and D by candidates above it. `<lo>fp<hi>` is the interval of the pixel discriminator: `3fp8` means 0.3 < Fpix ≤ 0.8, `99fp10` means 0.99 < Fpix ≤ 1.

For data, several histograms are only filled by candidates with Fpix ≤ 0.9, to keep the signal region blind.

The mass is computed from the momentum and `Ih` as `m = p · sqrt((Ih − C) / K)`. `K` and `C` are chosen from the dataset name at the beginning of `HSCPSelector::SlaveBegin`: one pair for 2024 data, one for simulation.

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
| `macros/macro.cc`, `kDefaultProdDir` | directory of the ntuples and of their file lists (overridden by the environment variable `HSCP_PROD_DIR`) |
| `launch<Family>.sh`, `CONFIG` | absolute path of `cfg/configFile.txt` |
| `macros/CombineHistos.C` | five absolute paths (search for `/safe/`) |

## Extending the analysis

- **New dataset**: add its name and location in `knownInputs()` (`macros/macro.cc`) and a line in `cfg/configFile.txt`. The selector decides from the dataset name which trigger to require, whether the sample is simulation or signal, and which `K` and `C` to use: check the `dataset_.find(...)` conditions in `HSCPTemplate.C`.
- **New selection**: add a row to `cfg/HSCPpreSelections.csv`. No code change is needed.
- **New histogram**: book it with `AddHisto1F` / `AddHisto2F` in `SlaveBegin` and fill it with `FillHisto1F` / `FillHisto2F` in `Process`, both in `HSCPTemplate.C`. The two names must match exactly: a fill with an unknown name is silently ignored.
- **New branch**: declare a `TTreeReaderValue` or `TTreeReaderArray` in `HSCPTemplate.h`.
- Any change in `src/` or `inc/` requires `./compile.sh`.