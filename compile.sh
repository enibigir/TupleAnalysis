#!/bin/bash
# Build libTools.so from src/ and inc/. Run from the root of the repository.
# The C++ standard comes from root-config, so the library matches the ROOT in use.
set -euo pipefail
cd "$(dirname "$0")"
source ./env.sh

CXXFLAGS="$(root-config --cflags) -O2 -fPIC -Wall -I inc"
for f in RegionMassPlot MassTools CPlots; do
    g++ $CXXFLAGS -c "src/$f.cc" -o "$f.o"
done
g++ -shared -fPIC CPlots.o MassTools.o RegionMassPlot.o $(root-config --glibs) -o libTools.so
rm -f CPlots.o MassTools.o RegionMassPlot.o
echo "libTools.so built with ROOT $(root-config --version)"
