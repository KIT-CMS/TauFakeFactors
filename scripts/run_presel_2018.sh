#!/bin/bash
cd "$(dirname "$0")/.." || exit 1

python preselection.py --config-file configs/btag_efficiency/2018/preselection_et.yaml --ncores 16
python preselection.py --config-file configs/btag_efficiency/2018/preselection_mt.yaml --ncores 16
python preselection.py --config-file configs/btag_efficiency/2018/preselection_tt.yaml --ncores 16
python preselection.py --config-file configs/btag_efficiency/2018/preselection_ee.yaml --ncores 16
python preselection.py --config-file configs/btag_efficiency/2018/preselection_em.yaml --ncores 16
python preselection.py --config-file configs/btag_efficiency/2018/preselection_mm.yaml --ncores 16
