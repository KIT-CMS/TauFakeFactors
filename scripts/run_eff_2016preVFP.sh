#!/bin/bash
# b-tag efficiency payloads of 2016preVFP, six channels, then the strict-consumer check against the era's BTV payload.
cd "$(dirname "$0")/.." || exit 1
set -e

for channel in et mt tt ee em mm; do
    python btag_efficiency.py --config-file configs/btag_efficiency/2016preVFP/btag_efficiency_${channel}.yaml --workers "${WORKERS:-4}" --threads "${THREADS:-4}"
done
python3 /work/sdaigler/bbtautau/check_btag_payload.py --era 2016preVFP --sf-check workdir/btag_eff_2016preVFP/2016preVFP/*/btag_efficiency.json.gz
