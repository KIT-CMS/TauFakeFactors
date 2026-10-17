#!/bin/bash
# Preselection of the 2017 b-tag efficiency probe jets (production smrun2_btag_eff_v1), six channels.
# preselection.py starts one worker per process and, inside it, one per sample (about 70 at once for
# the full list); NCORES (default 1) sets the ROOT threads of each, so the total stays below ~128.
cd "$(dirname "$0")/.." || exit 1
set -e

for channel in et mt tt ee em mm; do
    python preselection.py --config-file configs/btag_efficiency/2017/preselection_${channel}.yaml --ncores "${NCORES:-1}"
done
