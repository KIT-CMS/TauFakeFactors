#!/bin/bash
cd "$(dirname "$0")/.." || exit 1

python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_et.yaml --workers 4 --threads 4
python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_mt.yaml --workers 4 --threads 4
python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_tt.yaml --workers 4 --threads 4
python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_ee.yaml --workers 4 --threads 4
python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_em.yaml --workers 4 --threads 4
python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_mm.yaml --workers 4 --threads 4
