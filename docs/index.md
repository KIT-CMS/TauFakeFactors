# TauFakeFactors
FakeFactor framework for the estimation of jets misidentified taus with pyROOT.

## Setup
Clone the repository via
```bash
git clone --recurse-submodules https://github.com/KIT-CMS/TauFakeFactors.git
```
The environment can be set up with conda via
```bash
conda env create --file environment.yml
```

## General configuration
General definitions like paths for all steps of the fake factor measurements should be defined in a `configs/ANALYSIS/ERA/common_settings.yaml` file (this file name should always stay the same).

The expected ntuple folder structure is NTUPLE_PATH/ERA/SAMPLE_TAG/CHANNEL/*.root

  parameter | type | description
  ---|---|---
  `ntuple_path` | `string` | absolute path to the folder with the n-tuples on the dcache, a remote path is expected like "root://cmsdcache-kit-disk.gridka.de//store/user/USER/..."
  `friends` | `list` | (optional) list of friend names that exist for the n-tuples in `ntuple_path`
  `tree` | `string` | name of the tree in the n-tuple files ("ntuple" in CROWN)
  `era` | `string` | data taking era (e.g. "2018, "2017", "2016preVFP", "2016postVFP")
  `nanoAOD_version` | `string` | definition of the nanoAOD version is relevant to calculate the correct generator weights for the later specified simulated sample as well as getting the correct cross section
  `tau_vs_jet_wps` | `list` | list of tau ID vsJet working points to be written out in the preselection step (e.g. ["Medium", "VVVLoose"])
  `tau_vs_jet_wgt_wps` | `list` | list of tau ID vsJet working point scale factors to be written out in the preselection step (e.g. ["Medium"])

The output folder structure is OUTPUT_PATH/preselection/ERA/CHANNEL/*.root

  parameter | type | description
  ---|---|---
  `output_path` | `string` | absolute path where the files with the preselected events will be stored, a local path is expected like "/ceph/USER/..."
  `file_path` | `string` | absolute path to the folder with the preselected files (should be the same as `output_path`) to be used for the fake factor calculation
  `workdir_name` | `string` | relative path where the fake factor measurement output files will be stored; folder is produced in `workdir/`

The parameters in `common_settings.yaml` will be overwritten if they are also specified in the config for the individual steps of the fake factor determination, like event preselection, fake factor calculation and correction calculation.  

## Hints
* check out `configs/general_definitions.py`, this file has many relevant definition for plotting (dictionaries for names) and correctionlib output information
* check `ntuple_path` and `output_path` (preselection) and `file_path` and `workdir_name` (fake factors, corrections) in the used config files to avoid wrong inputs or outputs