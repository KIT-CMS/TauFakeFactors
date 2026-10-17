#!/bin/bash
# Install the b-tag efficiency payloads of one era into CROWN bbtautau (v7 spec 5.3):
#   workdir/btag_eff_<era>/<era>/<ch>/btag_efficiency.json.gz
#     -> $CROWN_BBTAUTAU/payloads/btagging_efficiencies/upart/<era>/btag_efficiency_<ch>.json.gz
# (the file name common_config.py builds per scope; 2018 replaces the v2 payloads). Nothing is copied unless
# all six payloads exist and check_btag_payload.py --era <era> --sf-check passes on them.
#   CROWN_BBTAUTAU  CROWN analysis config (default: the sibling KingMaker checkout)
#   WORKDIR         TauFakeFactors workdir (default: workdir)
#   SF_PAYLOAD      BTV payload for the check (default: the era's pin on cvmfs)
# Afterwards drop <era> from PENDING["btag_efficiency_payloads"] in CROWN tests/helpers.py and commit the
# payloads in CROWN.
cd "$(dirname "$0")/.." || exit 1
set -euo pipefail

era="${1:?usage: scripts/install_btag_payloads.sh <era>}"
crown="${CROWN_BBTAUTAU:-../KingMaker/CROWN/analysis_configurations/bbtautau}"
workdir="${WORKDIR:-workdir}"
channels=(et mt tt ee em mm)

payloads=()
for channel in "${channels[@]}"; do
    payload="${workdir}/btag_eff_${era}/${era}/${channel}/btag_efficiency.json.gz"
    if [[ ! -f "${payload}" ]]; then
        echo "install_btag_payloads: ${payload} is missing, nothing installed" >&2
        exit 1
    fi
    payloads+=("${payload}")
done
if [[ ! -d "${crown}/payloads" ]]; then
    echo "install_btag_payloads: no CROWN bbtautau config at ${crown} (set CROWN_BBTAUTAU), nothing installed" >&2
    exit 1
fi

check=(python3 scripts/check_btag_payload.py --era "${era}" --sf-check)
if [[ -n "${SF_PAYLOAD:-}" ]]; then
    check+=(--sf-payload "${SF_PAYLOAD}")
fi
if ! "${check[@]}" "${payloads[@]}"; then
    echo "install_btag_payloads: the check failed, nothing installed" >&2
    exit 1
fi

target="${crown}/payloads/btagging_efficiencies/upart/${era}"
mkdir -p "${target}"
for channel in "${channels[@]}"; do
    cp -v "${workdir}/btag_eff_${era}/${era}/${channel}/btag_efficiency.json.gz" \
        "${target}/btag_efficiency_${channel}.json.gz"
done
echo "installed into ${target}; drop ${era} from PENDING[\"btag_efficiency_payloads\"] in ${crown}/tests/helpers.py"
