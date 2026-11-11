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
export HSCP_PROD_DIR=/scratch/ui3_1/gcoulon/HSCP_prod/    # directory holding the ntuples and their file lists
ls $HSCP_PROD_DIR/JetMET2024/V12p310.txt                  # one file list per dataset, see cfg/datasets.json
head -3 $HSCP_PROD_DIR/JetMET2024/V12p310.txt             # the ROOT files must be readable from the worker nodes
```

If the ntuples are on storage that the worker nodes cannot mount, the file lists must hold
`root://` URLs instead of local paths.

## 3. Check the selections and the binning

```bash
cat cfg/HSCPpreSelections.csv      # one selection per row; memory is about 1.7 GB per selection
./launch.sh --list 'JetMET2024*'   # the samples, their type, dataset and binning, from cfg/datasets.json
```

## 4. Interactive test on one small dataset

```bash
./launch.sh TestMET2024 --dry-run                     # shows the root commands
./launch.sh TestMET2024 --out-dir output/test         # output/test/TestMET2024.root and .log
```

The same thing by hand, which is what the launcher and the Condor jobs run:

```bash
python3 python/CreateSelector.py
cd macros
root -l -b -q 'macro.cc("TestMET2024","../output/test/TestMET2024.files","TestMET2024.root",true,70,96,400,1000,400,20)'
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
python3 condor/submit.py --samples Mu2024C,Mu2024D --tag mu_CD
python3 condor/submit.py --dataset QCD2024 --tag qcd --cpus 4 --memory '6 GB'      # 4 implicit-MT threads
python3 condor/submit.py --samples TestMET2024 --tag t1 --requirements '(OpSysAndVer == "AlmaLinux9")'
```

Measure the memory on one real job before submitting a whole dataset: submit a single sample with
`--files-per-job 5`, then read `MemoryUsage` in `log/condor.log` and set `--memory` from it.

## 6. Submit and follow

```bash
cd condor/work/jetmet_v1
condor_submit hscp.sub             # or add --submit to the submit.py command above
condor_q                           # running / idle / held
condor_q -hold -af HoldReason      # jobs held because runJob.sh exited with an error
tail -f log/condor.log
less log/JetMET2024C_g0_c0.err
cd ../../..
```

## 7. Check completeness and resubmit what is missing

```bash
python3 condor/merge.py condor/work/jetmet_v1 --check
# missing jobs, if any, are in condor/work/jetmet_v1/resubmit.txt
cd condor/work/jetmet_v1 && cp resubmit.txt jobs.txt && condor_submit hscp.sub && cd ../../..
```

Repeat until `--check` prints `complete` for every sample.

## 8. Merge

```bash
source env.sh                      # hadd must be in the PATH
python3 condor/merge.py condor/work/jetmet_v1
ls condor/work/jetmet_v1/merged/   # one <sample>.root per sample
```

## 9. Hand the results to the plotting macros

```bash
mkdir -p output/JetMET2024_V12
cp condor/work/jetmet_v1/merged/JetMET2024*.root output/JetMET2024_V12/   # the directory name carries the run version
cd macros && root -l -b -q CombineHistos.C && cd ..      # enable the wanted plots in its main function
ls outputDisplay/
```

## 10. Trigger scale factors (when the calibration selections were run)

```bash
python3 python/DeriveSF.py         # needs numpy, uproot, matplotlib; writes inc/SF_*.txt
```

A new set of tables changes the weights of the simulation: rerun steps 5 to 8 for the
simulated samples afterwards.

## Clean-up

```bash
rm -rf condor/work/jetmet_v1       # chunks, generated selectors, logs and job outputs
rm -f macros/HSCPSelector_C*       # ACLiC products of an interactive run
```
