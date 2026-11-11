# `src/` and `inc/` — the `libTools` helper library

Reference for the seven source files in `src/` and `inc/`: what each class and function is for, and what its body does.

`inc/` also holds the three trigger scale-factor tables read by the selector (`SF_*.txt`), written by `python/DeriveSF.py`; their format is given in the `MassTools` section.

## Overview

| File | Content | In `libTools.so` |
|---|---|---|
| `inc/CPlots.h`, `src/CPlots.cc` | `CPlots`: named 1D/2D histograms booked and filled by string key | yes |
| `inc/RegionMassPlot.h`, `src/RegionMassPlot.cc` | `RegionMassPlot`: fixed set of control histograms for one ABCD region | yes |
| `inc/MassTools.h`, `src/MassTools.cc` | Mass from (p, Ih) with the K & C formula; scale-factor table readers | yes |
| `inc/ComputeATLASmass.h` | Mass from (p, Ih) with the ATLAS dE/dx parametrisation | no (header-only) |

### How the pieces are used

```
compile.sh ──► libTools.so  (CPlots, RegionMassPlot, MassTools)

macros/macro.cc
  ├─ gSystem->Load("../libTools.so")
  └─ chain->Process("HSCPSelector.C+")          # compiled on the fly by ACLiC
         └─ HSCPSelector.h
              ├─ #include "../inc/RegionMassPlot.h"   ─► MassTools.h
              ├─ #include "../inc/CPlots.h"
              └─ #include "../inc/ComputeATLASmass.h"
```

`HSCPSelector.{h,C}` are generated from `HSCPTemplate.{h,C}` by `python/CreateSelector.py`. The selector is the only consumer of the library. It does not include `MassTools.h` directly: `GetMass`, `loadSF1D` and `loadSF2D` reach it through `RegionMassPlot.h`, which is why that include must stay.

### Build

```
./compile.sh        # needs root-config in the PATH; writes libTools.so at the repository root
```

Any change in `src/` or `inc/` requires re-running `compile.sh`. `ComputeATLASmass.h` is recompiled with the selector, so it does not.

### Units and conventions

- Momenta in GeV, `Ih` in MeV/cm, masses in GeV.
- In `RegionMassPlot`, the argument and axes called `p` hold **10⁴/p** in GeV⁻¹, not the momentum.
- Histograms are detached from any `TDirectory` (`SetDirectory(nullptr)`) and booked with `Sumw2()`.
- Neither class owns its histograms: they are never deleted, and copying an object copies the pointers only. The selector hands them to `fOutput` at the end of the job.
- A value outside the range of an axis goes to the underflow or overflow bin; nothing is clamped.

---

## `CPlots` — `inc/CPlots.h`, `src/CPlots.cc`

A container of histograms addressed by name. The selector keeps one `CPlots` object per selection (`vcp`) plus one for the plots made before any selection (`vcp_nosel`).

Typical sequence: book with `AddHisto*` in `SlaveBegin`, fill with `FillHisto*` in `Process`, export with `AddToList(fOutput)` in `SlaveTerminate`.

### Data members

Three private maps, one per histogram type. The key is the histogram name.

| Member | Type |
|---|---|
| `mh1D_` | `std::map<std::string, TH1D*>` |
| `mh1F_` | `std::map<std::string, TH1F*>` |
| `mh2F_` | `std::map<std::string, TH2F*>` |

The three maps are independent: the same name can be booked once in each.

### `void AddHisto1F(name, nbins, xmin, xmax, title = "")`

1. If `name` is already a key of `mh1F_`, returns at once. The first booking wins; the new binning and title are ignored.
2. Creates `new TH1F(name, title, nbins, xmin, xmax)` and stores the pointer in `mh1F_[name]`. `title` is passed to ROOT unchanged, so the `"title;x title;y title"` syntax works.
3. Calls `SetDirectory(nullptr)`, then `Sumw2()`.
4. If `name` contains the text `_SignalMass`, calls `SetBinErrorOption(TH1::kPoisson)` and prints `Setting Poisson errors for histogram: <name>`.
5. Computes the bin width `(xmax − xmin) / nbins` and writes it with two decimals: in scientific notation if it is below 0.01 or above 1000, in fixed notation otherwise.
6. Sets the Y-axis title to `Normalized tracks / <bin width>`, for example `Normalized tracks / 25.00` or `Normalized tracks / 5.00e-03`. This replaces a Y title given through `title`.

### `void AddHisto1D(name, nbins, xmin, xmax, title = "")`

Same as `AddHisto1F` with a `TH1D` stored in `mh1D_`, without step 4: there is no `_SignalMass` case. The selector uses it for the cutflow histograms.

### `void AddHisto2F(name, nbinsx, xmin, xmax, nbinsy, ymin, ymax, title = "")`

1. If `name` is already a key of `mh2F_`, returns at once.
2. Creates `new TH2F(name, title, nbinsx, xmin, xmax, nbinsy, ymin, ymax)` and stores the pointer in `mh2F_[name]`.
3. Calls `SetDirectory(nullptr)`, then `Sumw2()`.

No axis title is set here; the axes only have the titles given through `title`.

### `bool FillHisto1F(name, value, weight = 1)`, `bool FillHisto1D(name, value, weight = 1)`, `bool FillHisto2F(name, xvalue, yvalue, weight = 1)`

Each one looks `name` up in the map of its own type.

- Found: calls `Fill(value, weight)` on the histogram (`Fill(xvalue, yvalue, weight)` for the 2D version) and returns `true`.
- Not found: does nothing and returns `false`.

The lookup is limited to one map: `FillHisto1D` on a name booked with `AddHisto1F` returns `false`.

### `bool AddToList(TList* list)`

1. Returns `false` if `list` is null.
2. Goes through `mh1D_`, then `mh1F_`, then `mh2F_`, each in alphabetical order of the names. A histogram is added to the list with `list->Add` unless the list already holds an object with the same name.
3. Returns `true`.

The list receives the pointers, not copies. Calling the function twice with the same list adds nothing the second time.

### Constructor and destructor

Both are empty. The histograms are created with `new` and never deleted by `CPlots`; a copy of a `CPlots` holds the same pointers as the original.

### Things to know

- **A misspelt name is silent.** `FillHisto*` returns `false` for an unknown name and the selector never checks the return value.
- **`kPoisson` and weights.** ROOT ignores the `kPoisson` option on a histogram whose sum of weights differs from its sum of squared weights. The `_SignalMass` histograms are filled with `weightPU × SF`, so they keep the usual `sqrt(sumw2)` errors; the option only acts on a histogram filled with unit weights.
- **Float interface.** All values go through `float` arguments, including the `double` branches of the tree.

---

## `RegionMassPlot` — `inc/RegionMassPlot.h`, `src/RegionMassPlot.cc`

A fixed set of 26 control histograms for one region of the background estimate. The selector builds 49 of them per selection, one for each of the regions A, B, C, D in each of their Fpix slices.

All histogram names end with the `suffix` given to the constructor, which the selector builds as `_region<A|B|C|D>_<lo>fp<hi>_<selection label>`, for example `ih_eta_regionD_8fp9_METanalysis_TestPUppiMETCut_Eta2p4`.

### Data members

All public.

| Members | Content |
|---|---|
| `suffix_` | the suffix given to the constructor |
| `np`, `plow`, `pup` | binning of the 10⁴/p axis |
| `npt`, `ptlow`, `ptup` | binning of the pT axis |
| `nih`, `ihlow`, `ihup` | binning of the Ih axis |
| `nias`, `iaslow`, `iasup` | binning of the Ias axis |
| `neta`, `etalow`, `etaup` | binning of the η axis |
| `nmass`, `masslow`, `massup` | binning of the mass axis |
| `nfpix`, `fpixlow`, `fpixup` | binning of the Fpix axis |
| 25 `TH2F*` and 1 `TH1F*` | the histograms, listed [below](#histograms) |

### `RegionMassPlot(suffix, etabins, ihbins, pbins, massbins, fpixbins, C_parameter)`

1. Sets the 26 histogram pointers to null.
2. Copies `suffix` into `suffix_`.
3. Calls `initHisto(etabins, ihbins, pbins, massbins, fpixbins, C_parameter)`.

The selector passes the bin counts of `cfg/datasets.json` and its `C` value as `C_parameter`.

### `void initHisto(etabins, ihbins, pbins, massbins, fpixbins, C_parameter)`

1. Calls `TH1::AddDirectory(false)` and `TH2::AddDirectory(false)`. This is a process-wide ROOT setting: from then on, no new histogram is attached to the current directory.
2. Fills the binning members from the arguments, with the ranges of the table below.
3. Books the 26 histograms with `new`, under the name `<name><suffix>`. Right after each booking, calls `Sumw2()` and `SetDirectory(nullptr)`.
4. Calls `SetBinErrorOption(TH1::kPoisson)` on `mass`.

It is public but meant for the constructor only: a second call replaces the 26 pointers with new histograms of the same names.

#### Binning

| Axis | Bins | Range |
|---|---|---|
| 10⁴/p | `np` = `pbins` | 0 – 200 GeV⁻¹ |
| pT (also used for p in `eta_p`) | `npt` = `pbins` | 0 – 10000 GeV |
| Ih | `nih` = `ihbins` | `C_parameter` – `C_parameter + 5` MeV/cm |
| Ias | `nias` = `ihbins` | 0 – 1 |
| η | `neta` = `etabins` | −2.4 – 2.4 |
| mass | `nmass` = `massbins` | 0 – 4000 GeV |
| Fpix | `nfpix` = `fpixbins` | 0 – 1 |
| nhits | 20 | 0 – 20 |
| npv | 100 | 0 – 100 |
| σ(pT)/pT | `np` in `eta_pterrOpt`, 100 in `pt_pterroverpt` | 0 – 1 |

The bin counts come from `cfg/datasets.json` (currently η 96, Ih 400, p 1000, mass 400, Fpix 20). With these values one `RegionMassPlot` takes about 38 MB: each cell holds a `float` content and a `double` sum of squared weights.

### `void fill(eta, nhits, p, pt, pterr, ih, ias, m, npv, fpix, w)`

Makes one `Fill` call per histogram, 26 in total, all with the weight `w`. The exact call for each histogram is in the last column of the [histogram table](#histograms). Two quantities are computed inside the function: the momentum `pt*cosh(eta)` for `eta_p`, and the relative uncertainty `pterr/pt` for `eta_pterrOpt` and `pt_pterroverpt`. There is no test on the arguments.

| Argument | Meaning | What the selector passes |
|---|---|---|
| `eta` | track η | `IsoTrack_eta` |
| `nhits` | number of dE/dx measurements | `DeDx_NoL1NOM` |
| `p` | **10⁴/p** [GeV⁻¹] | `10000 / (pT·cosh η)` |
| `pt`, `pterr` | pT and its uncertainty [GeV] | `IsoTrack_PseudoTrack_pt`, `IsoTrack_ptError` |
| `ih` | Ih [MeV/cm] | `DeDx_IhStrip` |
| `ias` | discriminator filled in the `ias` histograms | `DeDx_GiStrip` |
| `m` | reconstructed mass [GeV] | `GetMass(p, Ih, K, C)` |
| `npv` | number of good primary vertices | `PV_npvsGood` |
| `fpix` | pixel discriminator | `DeDx_FiPixelNoL1` |
| `w` | event weight | 1 for data; for simulation, pile-up weight times the PUppi MET scale factor |

All arguments are `float`. `GetMass` returns −1 when `Ih < C`: such a candidate lands in the underflow bin of the mass axis.

### Histograms

The object name is `<name><suffix>`. Three members have a name that differs from the object name.

| Member | Object name | X | Y | Call in `fill` |
|---|---|---|---|---|
| `ih_pt` | `ih_pt` | pT | Ih | `Fill(pt, ih, w)` |
| `ias_pt` | `ias_pt` | pT | Ias | `Fill(pt, ias, w)` |
| `ih_ias` | **`ias_ih`** | Ias | Ih | `Fill(ias, ih, w)` |
| `ih_fpix` | **`fpix_ih`** | Fpix | Ih | `Fill(fpix, ih, w)` |
| `eta_fpix` | **`fpix_eta`** | Fpix | η | `Fill(fpix, eta, w)` |
| `oP_fpix` | `oP_fpix` | 10⁴/p | Fpix | `Fill(p, fpix, w)` |
| `ih_nhits` | `ih_nhits` | nhits | Ih | `Fill(nhits, ih, w)` |
| `ias_nhits` | `ias_nhits` | nhits | Ias | `Fill(nhits, ias, w)` |
| `eta_pt` | `eta_pt` | pT | η | `Fill(pt, eta, w)` |
| `eta_1oP` | `eta_1oP` | 10⁴/p | η | `Fill(p, eta, w)` |
| `eta_p` | `eta_p` | p = pT·cosh η | η | `Fill(pt*cosh(eta), eta, w)` |
| `eta_pterrOpt` | `eta_pterrOpt` | η | σ(pT)/pT | `Fill(eta, pterr/pt, w)` |
| `nhits_pt` | `nhits_pt` | pT | nhits | `Fill(pt, nhits, w)` |
| `eta_nhits` | `eta_nhits` | nhits | η | `Fill(nhits, eta, w)` |
| `ih_eta` | `ih_eta` | η | Ih | `Fill(eta, ih, w)` |
| `ih_p` | `ih_p` | 10⁴/p | Ih | `Fill(p, ih, w)` |
| `ias_p` | `ias_p` | 10⁴/p | Ias | `Fill(p, ias, w)` |
| `pt_pterroverpt` | `pt_pterroverpt` | pT | σ(pT)/pT | `Fill(pt, pterr/pt, w)` |
| `ias_eta` | `ias_eta` | η | Ias | `Fill(eta, ias, w)` |
| `mass_eta` | `mass_eta` | η | mass | `Fill(eta, m, w)` |
| `eta_npv` | `eta_npv` | npv | η | `Fill(npv, eta, w)` |
| `p_npv` | `p_npv` | npv | 10⁴/p | `Fill(npv, p, w)` |
| `ih_npv` | `ih_npv` | npv | Ih | `Fill(npv, ih, w)` |
| `mass` | `mass` (TH1F, `kPoisson`) | mass | | `Fill(m, w)` |
| `mass_p` | `mass_p` | 10⁴/p | mass | `Fill(p, m, w)` |
| `mass_ih` | `mass_ih` | Ih | mass | `Fill(ih, m, w)` |

Axis titles that do not match the content: `mass_p` and `mass_ih` carry `Mass [GeV]` on X while the mass is on Y; `p_npv` is labelled `p [GeV]` but holds 10⁴/p; `oP_fpix` gives 10⁴/p in GeV instead of GeV⁻¹.

### `void addToList(TList* list)`

The function relies on a local helper, `safeAdd(object)`, which adds an object to the list unless

- the pointer is null, or
- the object name contains `CalibPseudoMET`, or
- the list already holds an object with the same name.

Steps:

1. Calls `safeAdd` on each of the 26 histograms.
2. Builds five derived objects from the current content of the histograms, detaches each with `SetDirectory(nullptr)` and passes it to `safeAdd`:

| Object | Built with | Content |
|---|---|---|
| `ih_pt<suffix>_py` (TH1D) | `ih_pt->ProjectionY` | Ih distribution, summed over all pT bins, underflow and overflow included |
| `ih_p<suffix>_pfx` (TProfile) | `ih_p->ProfileX` | mean Ih in each 10⁴/p bin |
| `ias_p<suffix>_pfx` (TProfile) | `ias_p->ProfileX` | mean Ias in each 10⁴/p bin |
| `ih_pt<suffix>_pfx` (TProfile) | `ih_pt->ProfileX` | mean Ih in each pT bin |
| `ias_pt<suffix>_pfx` (TProfile) | `ias_pt->ProfileX` | mean Ias in each pT bin |

The profiles use the Y bins inside the axis range only: entries in the Y underflow and overflow do not enter the mean.

One call therefore puts 31 objects in the list. The derived objects are a snapshot taken at the time of the call. There is no test on `list`, which must not be null. No selection of the current `cfg/HSCPpreSelections.csv` has `CalibPseudoMET` in its label, so nothing is skipped for that reason today.

### `~RegionMassPlot()`

Empty on purpose: it must not delete the histograms. The selector copies these objects by value (`push_back`, `for(auto obj : ...)`), and the copies share the histogram pointers, which also end up in `fOutput`. Deleting them here would destroy histograms that are still in use.

---

## `MassTools` — `inc/MassTools.h`, `src/MassTools.cc`

Three free functions.

### `float GetMass(float p, float ih, float k, float c)`

K & C mass. One line:

- if `ih − c < 0`, returns −1;
- otherwise returns `sqrt((ih − c) / k) · p`.

`ih = c` gives 0. The computation is done in `float`. The K and C values themselves are not in the library: they are globals of `HSCPTemplate.C`, set from the dataset name in `SlaveBegin`.

### `void loadSF1D(filepath, SF_MET, SF_Down, SF, SF_Up)`

Reads a 4-column scale-factor table into four `std::vector<float>`.

1. Clears the four vectors.
2. Opens `filepath`. If the file cannot be opened, throws `std::runtime_error("Cannot open file: <filepath>")`.
3. Prints `Loading SF from file: <filepath>`.
4. Reads four numbers at a time, separated by any white space, and appends them to `SF_MET`, `SF_Down`, `SF` and `SF_Up` in that order. Stops at the end of the file or at the first group that does not parse as four numbers.

### `void loadSF2D(filepath, SF_PseudoMETlo, SF_PseudoMEThi, SF_PUppiMETlo, SF_PUppiMEThi, SF_Down, SF, SF_Up)`

Same body for a 7-column table and seven vectors, filled in the order of the arguments. The message is `Loading 2D SF from file: <filepath>`.

### Things to know about the readers

- Line breaks have no special role: the numbers are read as one stream.
- A file with a header line, or with text in the middle, is read up to that text without any error. An incomplete last row is dropped.
- An empty result is not an error for the reader.
- The selector does not catch the exception: a missing table stops the job.

### File formats

- 4 columns: `MET  SF_down  SF  SF_up`. `MET` is the lower edge of the bin; the selector takes the next row as the upper edge, and the last row has no upper edge.
- 7 columns: `PseudoMET_lo  PseudoMET_hi  PUppiMET_lo  PUppiMET_hi  SF_down  SF  SF_up`. The selector uses the first row whose two intervals contain the event, and the last row if none does.

The three tables, loaded by the selector in `Begin`:

| File | Reader |
|---|---|
| `inc/SF_PseudoMET.txt` | `loadSF1D` |
| `inc/SF_PUppiMET.txt` | `loadSF1D` |
| `inc/SF_orMETtrg_PUppiMET_VS_PseudoMET__TriggerEffCalib_table_plain.txt` | `loadSF2D` |

---

## `ComputeATLASmass.h` — `inc/ComputeATLASmass.h`

Mass reconstruction from (p, Ih) with an ATLAS-style parametrisation of dE/dx against βγ. It is used for the alternative `*_SignalMass_ATLAS*` histograms, filled for simulation only.

The file is header-only and not part of `libTools.so`. It has no include guard, defines non-inline functions and global arrays, defines the macros `TOLERANCE`, `MIN` and `MAX`, and contains `using namespace std;`. It can therefore be included in one translation unit only, which today is the selector.

### Model

```
dEdx(βγ) = p1 · X^(p2/2) · ln(1 + (p3·βγ)^p4) − p5      with  X = ( sqrt(βγ⁴ + 4βγ²) − βγ² ) / 2
```

The curve falls steeply at low βγ, reaches a minimum, then rises slowly. A given `Ih` above the minimum is therefore reached at two values of βγ. The code keeps the lower one, on the falling branch, and returns `m = p / βγ`. An `Ih` below the minimum has no solution.

### Constants and parameter sets

| Macro | Value | Role |
|---|---|---|
| `MIN` | 0.3 | lower bound of βγ |
| `MAX` | 10000 | upper bound of βγ |
| `TOLERANCE` | 10⁻⁶ | tolerance of the root search, in βγ |

Each parameter set has the five parameters p1…p5 and their 5×5 covariance matrix, written as global arrays.

| Key passed to `findMass` | Arrays | Fit | Minimum of the curve | dEdx at βγ = 0.3 |
|---|---|---|---|---|
| `"2024data"` | `ATLASfit_data2024`, `covMatrixATLASfit_data2024` | 2024 data | 2.89 MeV/cm at βγ = 6.2 | 164 MeV/cm |
| `"2024bkg"` | `ATLASfit_bckg2024`, `covMatrixATLASfit_bckg2024` | 2024 background simulation | 2.97 MeV/cm at βγ = 6.9 | 316 MeV/cm |
| `"2024glupion"` | `ATLASfit_glupion2024`, `covMatrixATLASfit_glupion2024` | `glupion` variant | 3.02 MeV/cm at βγ = 5.4 | 19.8 MeV/cm |

The last two columns are computed with the functions of this file. The three covariance matrices are positive definite.

### `struct BetaGammaMinResult`

What `findBetaGammaWithCovariance` returns.

| Member | Content |
|---|---|
| `bg_nominal`, `bg_up`, `bg_down` | βγ at the minimum of the curve: nominal, largest, smallest |
| `params_nom[6]`, `params_up[6]`, `params_down[6]` | the parameters that gave each of them, in the layout `{Ih, p1, p2, p3, p4, p5}` |

### `TMatrixDSym InitializeCovMatrix(double covMatrix[5][5])`

Creates a 5×5 `TMatrixDSym`, copies the 25 elements of the array into it and returns it.

### `double AtlasFunction(double* x, double* par)`

The function given to `TF1`. With `x[0] = βγ`, `par[0] = Ih` and `par[1..5] = p1..p5`, it returns

```
p1 · X^(p2/2) · ln(1 + (p3·βγ)^p4) − p5 − Ih
```

that is `dEdx(βγ) − Ih`. Its root is the βγ at which the curve equals the measured `Ih`.

### `double findMinimumX(const double* par)`

1. Builds a `TF1` of `AtlasFunction` on [`MIN`, `MAX`] and fixes its six parameters to `par[0..5]`.
2. Asks ROOT for the position of the minimum with `TF1::GetMinimumX()`, and evaluates the function there.
3. If that value is positive, returns −1: the whole curve is above `Ih`, there is no solution.
4. Otherwise returns the position of the minimum, obtained with a second call to `GetMinimumX()`.

`Ih` only shifts the function vertically. The position of the minimum does not depend on it; only the test of step 3 does.

### `BetaGammaMinResult findBetaGammaWithCovariance(double Ih, const double* FitParam, const TMatrixDSym& Cov)`

Finds the position of the minimum for the nominal parameters and for parameters shifted along the eigen-directions of the covariance matrix.

1. Builds the nominal array `p_nom = {Ih, p1, …, p5}` from `FitParam`.
2. Diagonalises `Cov` with `TMatrixDSymEigen`: eigenvalues λ_k and eigenvectors V_k, for k = 0…4.
3. Computes `bg_nominal = findMinimumX(p_nom)`. Starts with `min_up = min_down = bg_nominal`, and with two independent copies of `p_nom` as the `up` and `down` parameters.
4. For each of the five eigen-directions:
   - builds the shift `δ = V_k · sqrt(λ_k)`;
   - builds two parameter arrays, `p_up = p + δ` and `p_down = p − δ`, with the same `Ih`;
   - computes `x_up = findMinimumX(p_up)` and `x_down = findMinimumX(p_down)`;
   - if `x_up > min_up`, keeps `x_up` and `p_up` as the new `up`;
   - if `x_down < min_down`, keeps `x_down` and `p_down` as the new `down`.
5. Returns `bg_nominal` with `p_nom`, `bg_up = min_up` with its parameters, `bg_down = min_down` with its parameters.

Consequences of step 4:

- `up` is the largest position among the nominal and the five `+δ` shifts; `down` is the smallest among the nominal and the five `−δ` shifts. A `−δ` shift is never a candidate for `up`, nor a `+δ` shift for `down`.
- `up` and `down` each come from a single eigen-direction. The shifts are not combined.
- If no shift moves the minimum in the wanted direction, the result keeps the nominal position and the nominal parameters.
- A `+δ` curve without solution gives `x_up = −1`, which is never kept. A `−δ` curve without solution gives `x_down = −1`, which is always kept: `bg_down` is then −1.

### `double ZeroBrentMethod(double a, double b, double t, const double* params)`

Root of `AtlasFunction` between `a` and `b` with Brent's method, which mixes bisection and interpolation. The code follows the `zero` routine of J. Burkardt's `brent` library.

1. Builds a `TF1` of `AtlasFunction` on [a, b] with the six parameters of `params`.
2. Evaluates the function at both ends. It keeps three points during the search: `sb`, the current best estimate; `sa`, the previous one; and `c`, a point where the function has the opposite sign to the one at `sb`, so that the root stays between `sb` and `c`.
3. Repeats:
   - if the function is smaller in absolute value at `c` than at `sb`, swaps the two, so that `sb` is always the best estimate;
   - computes the tolerance `tol = 2·ε·|sb| + t`, with ε = 2.22·10⁻¹⁶, and the half-width `m = (c − sb) / 2`;
   - stops if `|m| ≤ tol` or if the function is exactly zero at `sb`;
   - chooses the step. It tries an interpolation: a secant step when only two distinct points are known, an inverse quadratic one with three. The step is accepted if it stays within three quarters of the interval between `sb` and `c`, and is smaller than half of the step taken two iterations before. Otherwise, or when the previous steps were below `tol` or the last one did not reduce the function in absolute value, the step is a bisection, `m`;
   - moves `sb` by that step, or by `tol` towards `c` if the step is smaller than `tol`, and evaluates the function at the new point;
   - if the function now has the same sign at `sb` and at `c`, replaces `c` with the previous estimate.
4. Returns `sb`.

The function does not check that the two ends have opposite signs, and has no limit on the number of iterations. With `t = TOLERANCE`, βγ is found to about 10⁻⁶.

### `double findMass(const double p, const double Ih, const std::string year, bool up = false, bool down = false)`

Entry point. The selector only calls the nominal form, `findMass(p, Ih, key)`.

1. Selects the parameter array and builds the covariance matrix with `InitializeCovMatrix`, according to `year`. For any other key, prints an error message and returns −1.
2. Calls `findBetaGammaWithCovariance(Ih, params, CovMatrix)`.
3. Returns −1 if `bg_nominal`, `bg_up` or `bg_down` is negative.
4. Chooses the curve:
   - `up` is true: the upper end is `bg_up`, the parameters are `params_up`;
   - otherwise, `down` is true: `bg_down` and `params_down`;
   - otherwise: `bg_nominal` and `params_nom`.

   If both flags are true, `up` wins.
5. Finds βγ with `ZeroBrentMethod(MIN, <upper end>, TOLERANCE, <parameters>)`: the root between 0.3 and the minimum of the chosen curve.
6. Returns `p / βγ`.

When −1 is returned:

- unknown key;
- `Ih` below the minimum of the nominal curve;
- `Ih` below the minimum of one of the `−δ` curves. This raises the threshold slightly above the nominal minimum: about 2.90 MeV/cm for `2024data`, 2.97 for `2024bkg`, 3.02 for `2024glupion`. This also applies to the nominal form.

Things to know:

- **Very large `Ih`.** If `Ih` is above the value of the curve at βγ = 0.3 (last column of the parameter table), the root is not between the two ends. The search then converges to the lower end and the function returns `p / 0.3`. Only `2024glupion` has a reachable limit, 19.8 MeV/cm.
- **`up` and `down` refer to the position of the minimum in βγ**, not to the mass. Depending on `Ih`, the `up` mass is above or below the nominal one.
- **Cost.** Each call diagonalises the matrix and searches the minimum of eleven curves, although neither the eigen-directions nor the positions of the minima depend on `Ih` or `p`.
- **Error message.** For an unknown key it says `Use '2024' or '2018'`; the valid keys are those of the parameter table.
