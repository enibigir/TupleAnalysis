#!/bin/bash
# End-to-end smoke test on a synthetic ntuple, without any real data:
#   unit tests of python/hscp, build libTools, generate the selector, write two synthetic
#   files, run macro.cc in argument mode with a small binning, run launch.sh, check the
#   histogram entries, then do the same through the Condor tooling (submit.py + runJob.sh
#   locally + merge.py). Needs ROOT only. Everything is written under tests/work/.
set -euo pipefail
cd "$(dirname "$0")/.."
REPO=$PWD
source ./env.sh

# step "<name>": "completed" for the previous step, "processing" for this one; a failing command prints FAILED
STEP=""
step() {
    [ -n "$STEP" ] && echo "completed  $STEP"
    STEP=$1
    echo "processing $STEP ..."
}
trap 'rc=$?; [ $rc -ne 0 ] && echo "FAILED     $STEP (exit code $rc)"' EXIT

W=$REPO/tests/work
D=$REPO/tests/datasets.json                # the synthetic sample TestMET2024, under $W/prod
rm -rf "$W" && mkdir -p "$W/prod/JetMET2024"

step "0. unit tests"
python3 -m unittest discover -s tests

step "1. libTools"
./compile.sh

step "2. selector"                      # generated from cfg/selections.json
python3 python/CreateSelector.py > /dev/null

step "3. synthetic ntuples"
for i in 0 1; do root -l -b -q "tests/MakeSyntheticNtuple.C(\"$W/prod/synthetic_$i.root\", $i)" > /dev/null; done
ls "$W"/prod/synthetic_*.root > "$W/prod/JetMET2024/testJetMET.txt"   # TestMET2024 in tests/datasets.json

step "4. macro.cc"                      # argument mode, both files, small binning
( cd macros && root -l -b -q "macro.cc(\"TestMET2024\",\"$W/prod/JetMET2024/testJetMET.txt\",\"$W/direct.root\",true,70,12,20,20,20,10)" > "$W/direct.log" 2>&1 ) || { tail -30 "$W/direct.log"; exit 1; }
root -l -b -q "tests/CheckSmoke.C(\"$W/direct.root\", 2, true)"

step "4b. launch.sh"                    # same files, default binning of tests/datasets.json (takes longer)
HSCP_PROD_DIR=$W/prod python3 python/launch.py TestMET2024 --config "$D" --list
HSCP_PROD_DIR=$W/prod python3 python/launch.py TestMET2024 --config "$D" --dry-run > /dev/null
HSCP_PROD_DIR=$W/prod ./launch.sh TestMET2024 --config "$D" --out-dir "$W/launch" --label smoke > "$W/launch.log" 2>&1 || { tail -30 "$W/launch.log"; exit 1; }
root -l -b -q "tests/CheckSmoke.C(\"$W/launch/TestMET2024_smoke.root\", 2, true)"

step "5. submit.py"                     # prepare the Condor jobs: 1 file per job, 6 selections per group
rm -rf condor/work/smoke
HSCP_PROD_DIR=$W/prod python3 condor/submit.py --config "$D" --samples TestMET2024 --tag smoke --files-per-job 1 --selections-per-job 6
WORK=$REPO/condor/work/smoke

step "6. Condor jobs, run locally"      # as Condor would: fresh scratch directory, transferred inputs only
while read -r sample chunk group out nosel ptcut eta ih p mass fpix chunkid; do
    S=$W/scratch_g${group}_c${chunkid}; mkdir -p "$S"
    cp -r "$REPO/env.sh" "$REPO/compile.sh" "$REPO/inc" "$REPO/src" "$WORK/sel/g$group/macros" "$chunk" "$S/"
    ( cd "$S" && "$REPO/condor/runJob.sh" "$sample" "$chunk" "$group" "$out" "$nosel" 70 12 20 20 20 10 > job.log 2>&1 ) || { tail -30 "$S/job.log"; exit 1; }
    mv "$S/$out" "$WORK/out/"
done < "$WORK/jobs.txt"

step "7. merge.py"
python3 condor/merge.py "$WORK"
root -l -b -q "tests/CheckSmoke.C(\"$WORK/merged/TestMET2024.root\", 2, true)"
# group 1 holds the SigmaPtoverPt selections and must not hold the no-selection histograms
root -l -b -q "tests/CheckSmoke.C(\"$WORK/out/TestMET2024_g1_c0.root\", 1, false, \"METanalysis_TestPUppiMETCut_SigmaPtoverPt_0p5\")"
echo "completed  $STEP"
echo "ALL SMOKE STEPS PASSED"
