#!/bin/bash

eras=("2022preEE" "2022postEE" "2023preBPix" "2023postBPix")
channels=("et" "mt" "tt")

for era in "${eras[@]}"; do
    for channel in "${channels[@]}"; do
        echo "Running btag efficiency calculation for era=${era}, channel=${channel}"
        python btag_efficiency.py \
            --config-file "configs/btag_efficiency/${era}/btag_efficiency_${channel}.yaml"
    done
done
