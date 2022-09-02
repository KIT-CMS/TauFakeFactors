# TauFakeFactors
Reimplementation of the FakeFactor calculation with python.

## Setup
For this calculations ROOT is needed (tested version is 6.26). The environment is set up with conda.
```bash
conda create -c conda-forge --name root_env root python=3.9
```

## Event preselection
The inputs are ntuples which are produced with CROWN. All information for the preselection are defined in a configuration in the `configs` folder. The preselection is then performed with a python script:
```bash
python preselection.py --config CONFIG_NAME 
```