#!/bin/bash

eras=("2022preEE" "2022postEE" "2023preBPix" "2023postBPix")
channels=("et" "mt" "tt")

for era in "${eras[@]}"; do
    for channel in "${channels[@]}"; do
        echo "Running preselection for era=${era}, channel=${channel}"
        python preselection.py \
            --config-file "configs/btag_efficiency/${era}/preselection_${channel}.yaml"
    done
done
