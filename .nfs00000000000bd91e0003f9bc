#!/bin/bash

eras=("2022preEE" "2022postEE" "2023preBPix" "2023postBPix" "2024" "2025")
# eras=("2023postBPix")
# eras=("2024" "2025")
channels=("et" "mt" "tt" "em")
#channels=("em")

for era in "${eras[@]}"; do
    for channel in "${channels[@]}"; do
        echo "Running btag efficiency calculation for era=${era}, channel=${channel}"
        python btag_efficiency.py \
            --config-file "configs/btag_efficiency/${era}/btag_efficiency.yaml"
    done
done
