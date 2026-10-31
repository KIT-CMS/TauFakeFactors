#!/bin/bash
# b-tag efficiency payloads of 2016postVFP, six channels, then the strict-consumer check against the era's BTV payload and,
# if it passes, the installation into CROWN bbtautau (scripts/install_btag_payloads.sh; INSTALL=0: check only).
# Order (v7 spec 7.2): scripts/run_presel_2016postVFP.sh, then python3 scripts/btag_binning_scan.py --config-dir
# configs/btag_efficiency/2016postVFP until it reports no cell below 5 jets (flavor_bins updated), then this script.
cd "$(dirname "$0")/.." || exit 1
set -e

for channel in et mt tt ee em mm; do
    python btag_efficiency.py --config-file configs/btag_efficiency/2016postVFP/btag_efficiency_${channel}.yaml --workers "${WORKERS:-4}" --threads "${THREADS:-4}"
done
if [[ "${INSTALL:-1}" == 1 ]]; then
    scripts/install_btag_payloads.sh 2016postVFP
else
    python3 scripts/check_btag_payload.py --era 2016postVFP --sf-check workdir/btag_eff_2016postVFP/2016postVFP/*/btag_efficiency.json.gz
fi
