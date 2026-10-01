#!/bin/bash
cd "$(dirname "$0")/.." || exit 1

python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_et.yaml --workers 8 --threads 8
python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_mt.yaml --workers 8 --threads 8
python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_tt.yaml --workers 8 --threads 8
python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_ee.yaml --workers 8 --threads 8
python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_em.yaml --workers 8 --threads 8
python btag_efficiency.py --config-file configs/btag_efficiency/2018/btag_efficiency_mm.yaml --workers 8 --threads 8
