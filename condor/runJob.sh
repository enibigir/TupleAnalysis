#!/bin/bash
# HTCondor job: build libTools and the selector in the scratch directory, run macro.cc on one
# chunk of files for one group of selections, leave the output file in the scratch directory.
#
# Arguments (written by condor/submit.py in jobs.txt):
#   sample chunkList group outFile fillNosel ptcut etabins ihbins pbins massbins fpixbins
#
# Transferred inputs: env.sh compile.sh inc/ src/ macros/ (macro.cc + the generated selector
# of the group) and the chunk list. Nothing is read from the repository checkout.
set -euo pipefail

SAMPLE=$1; CHUNK=$2; GROUP=$3; OUT=$4; NOSEL=$5
PTCUT=$6; ETABINS=$7; IHBINS=$8; PBINS=$9; MASSBINS=${10}; FPIXBINS=${11}

SCRATCH=$PWD
echo "[job] host $(hostname), scratch $SCRATCH, $(date)"
echo "[job] sample=$SAMPLE chunk=$CHUNK group=$GROUP nosel=$NOSEL"

source ./env.sh                              # ROOT: installed, or LCG view matching the OS
./compile.sh                                 # libTools.so, built against that ROOT

export HSCP_SF_DIR="$SCRATCH/inc"
export HSCP_NTHREADS="${HSCP_NTHREADS:-1}"   # set from request_cpus by submit.py

CHUNK_ABS="$SCRATCH/$(basename "$CHUNK")"
OUT_ABS="$SCRATCH/$OUT"

cd macros
root -l -b -q "macro.cc(\"$SAMPLE\",\"$CHUNK_ABS\",\"$OUT_ABS\",$NOSEL,$PTCUT,$ETABINS,$IHBINS,$PBINS,$MASSBINS,$FPIXBINS)"
cd "$SCRATCH"

# Fail loudly if nothing was written: Condor then keeps the job in the log as failed.
test -s "$OUT_ABS"
echo "[job] wrote $OUT ($(du -h "$OUT_ABS" | cut -f1)), $(date)"
