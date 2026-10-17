#!/bin/bash
# b-tag efficiency payloads of 2018, six channels, then the strict-consumer check against the era's BTV payload.
cd "$(dirname "$0")/.." || exit 1
set -e

for channel in et mt tt ee em mm; do
    python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_${channel}.yaml --workers "${WORKERS:-4}" --threads "${THREADS:-4}"
done
python3 scripts/check_btag_payload.py --era 2018 --sf-check workdir/btag_eff_2018/2018/*/btag_efficiency.json.gz
