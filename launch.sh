#!/bin/bash
# Run the selector on datasets: ./launch.sh 'JetMET2024*' Mu2024C   (see python/launch.py --help)
# Provides ROOT through env.sh, then hands over to python/launch.py.
set -euo pipefail
cd "$(dirname "$0")"
source ./env.sh
exec python3 python/launch.py "$@"
