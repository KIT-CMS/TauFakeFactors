#!/bin/bash

eras=("2018" "2022preEE" "2022postEE" "2023preBPix" "2023postBPix" "2024" "2025")
channels=("et" "mt" "tt" "em")

for era in "${eras[@]}"; do
    for channel in "${channels[@]}"; do
        config="configs/btag_efficiency/${era}/preselection_${channel}.yaml"
        if [ ! -f "${config}" ]; then
            echo "Skipping era=${era}, channel=${channel}: no ${config}"
            continue
        fi
        echo "Running preselection for era=${era}, channel=${channel}"
        python preselection.py --config-file "${config}"
    done
done
