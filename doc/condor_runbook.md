# Runbook: from a fresh checkout to merged histograms on HTCondor

Every command below is run from the root of the repository unless a `cd` says otherwise.
No CMSSW: the toolchain is ROOT (installed, or the LCG view that `env.sh` picks from CVMFS) and Python 3.

## 1. Environment and build

```bash
git clone <this repository> TupleAnalysis && cd TupleAnalysis
source env.sh                      # uses the ROOT in the PATH, else sources the LCG view for el8/el9
root-config --version              # must print a ROOT 6 version
./compile.sh                       # libTools.so at the root of the repository
```

Optional, once per machine: the unit tests and the smoke test. The latter builds a synthetic
ntuple and runs the whole chain, including the launcher and the Condor scripts executed locally.

```bash
python3 -m unittest discover -s tests -v
./tests/run_smoke.sh
```

## 2. Site configuration

```bash
ls input/                          # one file list per sample, see cfg/datasets.json ("prod_dir": "" = the checkout)
head -3 input/JetMET2024C.txt      # the ROOT files must be readable from the worker nodes
```

If the ntuples are on storage that the worker nodes cannot mount, the file lists must hold
`root://` URLs instead of local paths.

For remote files (`davs://`, `root://`), the jobs need a grid proxy that outlives them:

```bash
voms-proxy-init -voms cms -valid 192:00
voms-proxy-info --timeleft         # seconds left
```

`condor/submit.py` adds it to `hscp.sub` by itself (see `--proxy` in section 5).

## 3. Check the selections and the binning

```bash
cat cfg/selections.json            # the cuts and the selections; memory is about 0.7 GB per selection
./launch.sh --list 'JetMET2024*'   # the samples, their type, dataset and binning, from cfg/datasets.json
```

## 4. Interactive test on one small dataset

```bash
./launch.sh TestJetMET2024C --dry-run                     # shows the root commands
./launch.sh TestJetMET2024C --out-dir output/test         # output/test/TestJetMET2024C.root and .log
```

The same thing by hand, which is what the launcher and the Condor jobs run:

```bash
python3 python/CreateSelector.py
cd macros
root -l -b -q 'macro.cc("TestJetMET2024C","../output/test/TestJetMET2024C.files","TestJetMET2024C.root",true,70,96,400,1000,400,20)'
cd ..
```

## 5. Prepare the Condor jobs

```bash
python3 condor/submit.py --dataset JetMET2024 --tag jetmet_v1 \
        --files-per-job 20 --selections-per-job 2 --memory '4 GB' --runtime 7200
```

This writes `condor/work/jetmet_v1/` with the chunked file lists, one generated selector per
group of selections, `jobs.txt` and `hscp.sub`. Useful variants:

```bash
python3 condor/submit.py --dataset Gluino_madgraph --tag gluino_v1 --selections-per-job 12
python3 condor/submit.py --samples QCD2024_mu_pt80to120,QCD2024_mu_pt120to170 --tag qcd_80to170
python3 condor/submit.py --dataset QCD2024 --tag qcd --cpus 4 --memory '6 GB'      # 4 implicit-MT threads
python3 condor/submit.py --samples TestJetMET2024C --tag t1 --requirements '(OpSysAndVer == "AlmaLinux9")'
```

Grid proxy: by default (`--proxy auto`), if an input file is remote, `submit.py` writes
`x509userproxy = <file>` into `hscp.sub`, taking the file from `$X509_USER_PROXY`, else
`/tmp/x509up_u<uid>`. It stops if the proxy is missing or lives less than `--runtime` + 1 hour.
Condor copies the proxy to the worker node. `--proxy PATH` uses another file, `--proxy none` adds nothing.

### Where the job outputs go (`--out-dir`)

By default the outputs come back to `condor/work/<tag>/out/`. A full submission writes
100–150 GB of job outputs, which belongs on EOS rather than in the checkout:

```bash
python3 condor/submit.py --dataset QCD2024 --tag qcd_v1 \
        --out-dir 'root://cmseos.fnal.gov//store/user/$USER/hscp/'
# outputs:     root://cmseos.fnal.gov//store/user/<you>/hscp/qcd_v1  (copied by the jobs)
```

The outputs go to `<out-dir>/<tag>/`; `$USER` is expanded by `submit.py` if the shell did not.
With a URL (`root://`, `davs://`, ...), each job copies its own output at the end (`xrdcp` for
`root://`, `gfal-copy` otherwise, three attempts), and Condor brings back only the logs. A failed
copy makes the job exit with code 4, so it is held like any other failure. The jobs then need the
grid proxy, which `--proxy auto` gives them because the output is remote. With a local directory
instead of a URL, Condor transfers each output there (the submit machine must see it).
`submit.py` records the location in `condor/work/<tag>/output_dir.txt`, where `merge.py` reads it.

Check that you can write there before submitting many jobs:

```bash
xrdfs root://cmseos.fnal.gov ls /store/user/$USER          # needs the proxy
```

### Where the Condor logs go (`--log-dir`)

By default the `.out`, `.err` and `condor.log` files go to `condor/work/<tag>/log/`. A few thousand
jobs can fill a small home quota ("(errno 122) Disk quota exceeded", e.g. a checkout in
`/uscms/home` on the LPC). `--log-dir L` writes them to `L/<tag>/` instead:

```bash
python3 condor/submit.py --dataset QCD2024 --tag qcd_v1 --log-dir ~/nobackup/hscp_logs \
        --out-dir 'root://cmseos.fnal.gov//store/user/$USER/hscp/'
# logs:        /uscms_data/d3/<you>/hscp_logs/qcd_v1
```

`L` must be a local directory that the schedulers see, like the checkout (below); a URL is refused,
since Condor writes the logs itself. The location is recorded in `condor/work/<tag>/log_dir.txt`.
Jobs held for the quota can be released (`condor_release <ClusterId>`) once space is freed.

The paths in `hscp.sub` are real paths (symlinks resolved), because the Condor schedulers read
them, not the login node. On the FNAL LPC, `submit.py` refuses a checkout, a local `--out-dir` or a `--log-dir`
outside `/uscms/home` and `/uscms_data` (`~/nobackup`): the schedulers cannot see `/tmp`,
`/storage/local` or `/eos/uscms`, and every job would be held with "Cannot access initial working
directory" (hold code 14).

#### On CERN LXPLUS

The personal EOS area is `/eos/user/<first letter>/<username>` on `eosuser.cern.ch`. It is
reached with Kerberos (the LXPLUS login), not the proxy, and has a 1 TB quota by default.

```bash
xrdfs root://eosuser.cern.ch ls /eos/user/${USER:0:1}/$USER     # like the cmseos listing above
ls /eos/user/${USER:0:1}/$USER                                  # EOS is also mounted on LXPLUS
eos ls /eos/user/${USER:0:1}/$USER

python3 condor/submit.py --dataset QCD2024 --tag qcd_v1 \
        --out-dir 'root://eosuser.cern.ch//eos/user/e/$USER/hscp/'   # write your own first letter
```

Write the first letter of the username directly: `submit.py` expands `$USER`, but not
`${USER:0:1}`. CERN's HTCondor passes the Kerberos credentials to the jobs, so `xrdcp` can write
to `eosuser`. The inputs are still read over `davs://`, so `--proxy auto` still requires a proxy.
CMS group areas are under `root://eoscms.cern.ch//eos/cms/store/group/...`, and need write
permission for the group.

Keep the checkout on AFS (`/afs/cern.ch/user/<letter>/$USER` or an AFS workspace), not on `/eos`:
CERN's HTCondor cannot use the EOS mount for the job directory, the logs or a local `--out-dir`, and
`submit.py` refuses it on LXPLUS. The outputs go to EOS by URL, as above. The working directories
(chunk lists, selectors, logs) take a few hundred MB of AFS quota for a few thousand jobs.

Keep the proxy on AFS. `lxplus.cern.ch` sends each login to another node, and `/tmp` is local to
the node. `hscp.sub` records the proxy path, so a later `condor_submit hscp.sub` (resubmission)
from another node would not find a proxy left in `/tmp`. `--proxy auto` follows
`$X509_USER_PROXY`:

```bash
voms-proxy-init -voms cms -valid 192:00 --out ~/private/x509up
export X509_USER_PROXY=~/private/x509up        # also in ~/.bashrc
```

Memory: CERN's batch slots come with about 2 GB per core. A job asking for 3 GB on one core may wait
longer before it starts, and may count as more than one core in your share. If the measured
`MemoryUsage` is above 2 GB, ask for 2 cores with the memory: `--cpus 2 --memory '4 GB'` (the
selector uses the second core through `HSCP_NTHREADS`).

Measure the memory on one real job before submitting a whole dataset: submit a single sample with
`--files-per-job 5`, then read `MemoryUsage` in `log/condor.log` and set `--memory` from it (in the
`--log-dir` if you gave one).

## 6. Submit and follow

```bash
cd condor/work/jetmet_v1
condor_submit hscp.sub             # or add --submit to the submit.py command above
condor_q                           # running / idle / held
condor_q -hold -af HoldReason      # jobs held because runJob.sh exited with an error
tail -f log/condor.log             # with --log-dir: tail -f $(cat log_dir.txt)/condor.log
less log/JetMET2024C_g0_c0.err
cd ../../..
```

## 7. Check completeness and resubmit what is missing

```bash
python3 condor/merge.py condor/work/jetmet_v1 --check
# missing jobs, if any, are in condor/work/jetmet_v1/resubmit.txt
cd condor/work/jetmet_v1 && cp resubmit.txt jobs.txt && condor_submit hscp.sub && cd ../../..
```

Repeat until `--check` prints `complete` for every sample. With a remote `--out-dir`, `merge.py`
lists the directory with `xrdfs` (`root://`) or `gfal-ls` (other URLs); a job whose copy failed
has no file there and is listed as missing.

## 8. Merge

```bash
source env.sh                      # hadd must be in the PATH
python3 condor/merge.py condor/work/jetmet_v1
ls condor/work/jetmet_v1/merged/   # one <sample>.root per sample
```

The merged files are always written locally, in `condor/work/<tag>/merged/`. With a remote
`--out-dir`, `hadd` reads the job outputs directly from the storage (it needs the proxy).

`merge.py` finds the job outputs through `condor/work/<tag>/output_dir.txt`. To read them from
somewhere else (outputs moved or copied after the jobs), give the full directory, `<tag>` included:

```bash
python3 condor/merge.py condor/work/jetmet_v1 --check \
        --input-dir 'root://eosuser.cern.ch//eos/user/e/$USER/hscp/jetmet_v1'
```

`--input-dir` takes a local directory or a URL, like `--out-dir`; `$USER` and `~` are expanded.

## 9. Weight the simulation by cross-section

`python/Weight.py` scales every histogram of each simulated sample by

```
factor = xsec [pb] × luminosity [pb^-1] / sum of weights
```

then adds the samples into one file. It reads `<input-dir>/<sample>.root`, so it works on the
`merged/` directory of step 8 as well as on the `--out-dir` of `launch.sh`.

| Input | Where it comes from |
|---|---|
| `xsec` | the sample's `xsec` (pb) in `cfg/datasets.json`; for a filtered sample (QCD μ-enriched), σ × filter efficiency |
| luminosity | `--lumi` (pb^-1) if given, else the top-level `luminosity` of `cfg/datasets.json` |
| sum of weights | the sample's `sum_weights` in `cfg/datasets.json` if given, else bin 1 of its `EventCutflow` histogram (the sum of `weightPU` over every event read) |

Data samples are only added (factor 1). Data and simulation cannot go into the same output.

```bash
source env.sh                                      # PyROOT

# 1. check the factors: nothing is written
python3 python/Weight.py --dataset QCD2024 --input-dir condor/work/qcd/merged --dry-run

# 2. write the weighted file
python3 python/Weight.py --dataset QCD2024 --input-dir condor/work/qcd/merged \
        -o output/QCD2024_V16/QCD2024_mu_V16p2_weighted.root

# another luminosity than the one of datasets.json (era C only, for instance), in pb^-1
python3 python/Weight.py --dataset QCD2024 --input-dir condor/work/qcd/merged --lumi 7000 \
        -o output/QCD2024_V16/QCD2024_mu_eraC_weighted.root

# single samples or patterns instead of a whole dataset
python3 python/Weight.py --samples 'TTbar*' --input-dir condor/work/ttbar/merged -o output/TTbar2024/TTbar2024_weighted.root

# data: the eras are added, not scaled
python3 python/Weight.py --dataset JetMET2024 --input-dir condor/work/jetmet_v1/merged \
        -o output/JetMET2024_V12/JetMET2024_V12.root
```

| Option | Meaning |
|---|---|
| `--dataset NAME` / `--samples PATTERNS` | the samples to add, as in `submit.py` (`--dataset` can be repeated) |
| `--input-dir DIR` | the directory holding `<sample>.root` (required) |
| `-o FILE` | output; by default `<input-dir>/<dataset>_weighted.root` when all samples belong to one dataset |
| `--lumi L` | luminosity in pb^-1 for this run, instead of the one in `datasets.json` |
| `--dry-run` | print the table of factors only |
| `-f` | overwrite an existing output |
| `--skip-missing` | leave out the samples whose `<sample>.root` is missing (listed as a warning) instead of stopping |
| `--config FILE` | another datasets file |

The dry run prints, for each sample, `xsec`, the sum of weights and where it comes from
(`EventCutflow` or `datasets.json`), the factor, and the luminosity used. Check that:

- the factors of a binned sample (the 12 QCD bins) change smoothly from one bin to the next;
- the luminosity line is the one you meant (pb^-1: 7 fb^-1 is `--lumi 7000`).

The default sum of weights is right only if the ntuples hold every event of the dataset (no skim)
and the generator weights are all 1. Pythia LO samples (QCD) are fine. NLO samples (tt̄ Powheg,
W+jets) have negative weights that the selector does not read: set their `sum_weights` in
`cfg/datasets.json` to the sum of generator weights of the events processed. A sample with missing
jobs is still normalized correctly by the default, which counts only the events read.

Profiles (`*_pfx`) are added with the factors as weights, so their mean values are not scaled. The
file is written as `<output>.part` and renamed at the end, and histograms are processed one name
at a time: memory stays small (about 0.6 GB for 20,000 histograms).

### Several outputs in parallel

Each call takes about a minute whatever the size of its samples, because it writes every
histogram. Independent outputs can run at the same time, one process each (about 0.6 GB of memory
each). On the three test samples, the weighting took 74 s instead of 170 s, one after the other:

```bash
D=output/testchain
python3 python/Weight.py --samples TestQCD2024_mu_pt80to120 --input-dir $D -o $D/weighted/TestQCD2024_weighted.root -f > $D/weighted/TestQCD2024_weighted.log 2>&1 & p1=$!
python3 python/Weight.py --samples TestGluino_Run3_MET_madgraph_1600 --input-dir $D -o $D/weighted/TestGluino_1600_weighted.root -f > $D/weighted/TestGluino_1600_weighted.log 2>&1 & p2=$!
python3 python/Weight.py --samples TestJetMET2024C --input-dir $D -o $D/weighted/TestJetMET2024C.root -f > $D/weighted/TestJetMET2024C.log 2>&1 & p3=$!
wait $p1; e1=$?; wait $p2; e2=$?; wait $p3; e3=$?
echo "exit codes $e1 $e2 $e3"                    # 0 0 0
tail -qn1 $D/weighted/*.log                      # each: "20268 histograms written to ..."
```

Run the dry runs first, one after the other: they are quick and print the factors to check. Give
each process its own log, otherwise their outputs mix. Wait on the process IDs as above: `wait`
alone, or a loop over `jobs -p`, can return before the jobs end in some shells (zsh), and loses
their exit codes.

## 10. Hand the results to the plotting macros

The data/simulation validation plots of `macros/CombineHistos.C` are made by `python/plots.py` from
the files named in `cfg/plots.json`: name each input once (`data`, `qcd`, `ttbar`, `wjets`, `gluino`, ...),
the plots refer to these names.

```json
"lumi_label": "7.98 fb^{-1} (13.6 TeV)",
"inputs": {"data":   {"file": "output/JetMET2024_V1/JetMET2024.root", "label": "JetMET data 2024"},
           "qcd":    {"file": "output/QCD2024_V1/QCD2024_7p98fb_weighted.root", "label": "QCD (#mu enriched)"},
           "gluino": {"file": "output/gluino1600_v1/Gluino_1600_weighted.root", "label": "#tilde{g} (m=1600 GeV)"},
           "gluino_mass": "output/Gluino_V1/Gluino_Run3_MET_madgraph_%d_weighted.root", ...},
"plots":  [{"function": "nPV", "args": ["qcd", "ttbar", "wjets", "data"]}, ...]
```

An input is a path, or `{"file", "label"}` to set its legend entry (ROOT TLatex); `lumi_label` sets the
luminosity text of the plots that print one (`Cutflows`, `CompareKinematics`). The legends of `Cutflows`,
`nPV`, `nHSCP`, `FpixelPlot`, `Acceptance_EtaSlice_DataMC`, `CompareKinematics` (the three masses) and
`YieldAfterSF` (the signal) take these labels; an input without a label keeps the text of the macro.
`--list` shows the label of each input.

```bash
source env.sh
python3 python/plots.py --list                  # each plot, its inputs, MISSING for those not found
python3 python/plots.py                         # every plot; outputs in output_dir (outputDisplay/V1)
python3 python/plots.py --only 'nPV,Cutflows'   # some of them
python3 python/plots.py --jobs 4                # several at a time (about 5 s each to load the macro)
```

Each plot runs in its own `root` process (log in `<output_dir>/logs/<name>.log`), so one failure does
not stop the others; a plot whose input files are missing is skipped. Before running, `plots.py` warns
about histograms that a function reads by name and that none of its inputs holds: the function would
crash on them (`crashed (segmentation violation)`). That happens when the selections of
`cfg/selections.json` were renamed since the function was written: `nPV`, `nHSCP`, `FpixSlices` and
`ShowSignalEfficiency` read `METanalysis_Eta2p4_*` or `METanalysis_PseudoMETrescaled_*` histograms, and
`YieldAfterSF` the `CalibPseudoMET_*` ones of the calibration selections.

The converted functions are `Cutflows`, `nPV`, `nHSCP`, `FpixelPlot`, `FpixSlices`, `CompareMaping`,
`Acceptance_EtaSlice_DataMC`, `CompareKinematics`, `ShowSignalEfficiency` (one `%d` pattern for the
masses) and `YieldAfterSF`. Their old paths are the default values of their new parameters, so the calls
in `CombineHistos()` still work; the others keep their fixed paths:

```bash
cd macros && root -l -b -q CombineHistos.C && cd ..      # enable the wanted plots in its main function
```

## 11. Trigger scale factors (when the calibration selections were run)

```bash
python3 python/DeriveSF.py         # needs numpy, uproot, matplotlib; writes inc/SF_*.txt
```

A new set of tables changes the weights of the simulation: rerun steps 5 to 9 for the
simulated samples afterwards.

## Example: 2024 production (background, signal, JetMET2024C)

Steps 5 to 9 for the samples that are ready: the background samples with a cross-section, the
1600 GeV madgraph gluino, and JetMET2024C. One tag per group, so that each is checked, merged and
weighted on its own. The job outputs go to EOS. Memory and runtime are estimates (about 0.7 GB per
selection plus ROOT): submit the signal first, it is 6 small jobs, and read `MemoryUsage` and
`RemoteWallClockTime` before the large submissions.

```bash
source env.sh
voms-proxy-init -voms cms -valid 192:00               # if less than ~3 days are left
xrdfs root://cmseos.fnal.gov ls /store/user/$USER     # EOS reachable
OUT='root://cmseos.fnal.gov//store/user/$USER/hscp/'

# signal: Gluino 1600 (6 jobs)
python3 condor/submit.py --samples Gluino_Run3_MET_madgraph_1600 --tag gluino1600_v1 \
        --files-per-job 2 --selections-per-job 2 --memory '3 GB' --runtime 7200 --out-dir "$OUT" --submit

# background: QCD (12 bins), TTbar2024, TTbarSemiLep2024, WjetMuNu2024 (about 4,500 jobs)
python3 condor/submit.py \
        --dataset QCD2024 --dataset TTbar2024 --dataset TTbarSemiLep2024 --dataset WjetMuNu2024 \
        --tag bkg_v1 --files-per-job 20 --selections-per-job 2 --memory '3 GB' --runtime 7200 \
        --out-dir "$OUT" --submit

# data: JetMET2024C (about 300 jobs)
python3 condor/submit.py --samples JetMET2024C --tag jetmetC_v1 \
        --files-per-job 20 --selections-per-job 2 --memory '3 GB' --runtime 7200 --out-dir "$OUT" --submit
```

These are all the simulated samples of `cfg/datasets.json`: those without a cross-section were removed.
Without `--submit`, check `condor/work/<tag>/hscp.sub` first, then `condor_submit` it.

Follow, resubmit and merge:

```bash
condor_q
condor_q -hold -af HoldReason                         # exit 4: copy to EOS failed; 3: processing failed
condor_history -limit 5 -af MemoryUsage RemoteWallClockTime
for t in gluino1600_v1 bkg_v1 jetmetC_v1; do python3 condor/merge.py condor/work/$t --check; done
cd condor/work/bkg_v1 && cp resubmit.txt jobs.txt && condor_submit hscp.sub && cd ../../..   # per tag with missing jobs
for t in gluino1600_v1 bkg_v1 jetmetC_v1; do python3 condor/merge.py condor/work/$t; done
```

Weight the simulation with the luminosity of era C (`--lumi`, pb^-1), not the full-year
`luminosity` of `cfg/datasets.json`:

```bash
L=<luminosity of era C in pb^-1>
M=condor/work/bkg_v1/merged
python3 python/Weight.py --dataset QCD2024 --input-dir $M --lumi $L --dry-run       # factors smooth across the bins?
mkdir -p output/bkg output/signal
pids=() names=()
for d in QCD2024 TTbar2024 TTbarSemiLep2024 WjetMuNu2024; do       # the 5 outputs in parallel
    python3 python/Weight.py --dataset $d --input-dir $M --lumi $L -o output/bkg/${d}_weighted.root \
            > output/bkg/${d}_weighted.log 2>&1 &
    pids+=($!) names+=($d)
done
python3 python/Weight.py --samples Gluino_Run3_MET_madgraph_1600 --input-dir condor/work/gluino1600_v1/merged \
        --lumi $L -o output/signal/Gluino_1600_weighted.root > output/signal/Gluino_1600_weighted.log 2>&1 &
pids+=($!) names+=(Gluino_1600)
for i in "${!pids[@]}"; do
    wait "${pids[$i]}" && echo "ok    ${names[$i]}" || echo "FAIL  ${names[$i]} (see its .log)"
done
```

TTbar and WjetMuNu are NLO samples with negative weights: set their `sum_weights` in
`cfg/datasets.json` before weighting (step 9). JetMET2024C needs no weighting:
`condor/work/jetmetC_v1/merged/JetMET2024C.root` is used as it is.

If `MemoryUsage` stays well below 3 GB, use `--selections-per-job 4` or 6 (with more `--memory`)
for the next round: fewer jobs.

## Clean-up

```bash
rm -rf condor/work/jetmet_v1       # chunks, generated selectors, logs and job outputs (if in out/)
# job outputs on EOS (--out-dir), once merged:
eos root://cmseos.fnal.gov rm -r /store/user/$USER/hscp/jetmet_v1     # the eos client of the LPC nodes
eos rm -r /eos/user/${USER:0:1}/$USER/hscp/jetmet_v1                  # on LXPLUS
rm -f macros/HSCPSelector_C*       # ACLiC products of an interactive run
```
