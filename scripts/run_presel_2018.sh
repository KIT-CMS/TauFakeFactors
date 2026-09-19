#!/bin/bash
cd "$(dirname "$0")/.." || exit 1

python preselection.py --config-file configs/btag_efficiency/2018/preselection_et.yaml
python preselection.py --config-file configs/btag_efficiency/2018/preselection_mt.yaml
python preselection.py --config-file configs/btag_efficiency/2018/preselection_tt.yaml
python preselection.py --config-file configs/btag_efficiency/2018/preselection_ee.yaml
python preselection.py --config-file configs/btag_efficiency/2018/preselection_em.yaml
python preselection.py --config-file configs/btag_efficiency/2018/preselection_mm.yaml
