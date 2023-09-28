# TauFakeFactors
FakeFactor framework for the estimation of jets misidentified taus with pyROOT.

## Setup
The environment can be set up with conda via
```bash
conda env create --file environment.yaml
```

<details>
<summary>

## Event preselection

</summary>

This framework is designed for n-tuples produced with CROWN as input. 

All information for the preselection step is defined in a configuration file in the `configs/` folder. 
The preselection config has the following parameters:

* The expected input folder structure is NTUPLE_PATH/ERA/SAMPLE_TAG/CHANNEL/*.root
    parameter | type | description
    ---|---|---
    `ntuple_path` | `string` | absolute path to the folder with the n-tuples on the dcache, a remote path is expected like "root://cmsxrootd-kit.gridka.de//store/user/USER/..."
    `era` | `string` | data taking era ("2018, "2017", "2016preVFP", "2016postVFP")
    `channel` | `string` | tau pair decay channels ("et", "mt", "tt")
    `tree` | `string` | name of the tree in the n-tuple files ("ntuple" in CROWN)
    `analysis` | `string` | analysis name, needed to get the output features which are saved/needed for the later steps e.g. `"smhtt_ul"`

* The output folder structure is OUTPUT_PATH/preselection/ERA/CHANNEL/*.root
    parameter | type | description
    ---|---|---
    `output_path` | `string` | absolute path where the files with the preselected events will be stored, a local path is expected like "/ceph/USER/..."

* In `processes` all the processes are defined that should be preprocessed. \
  The names are also used for the output file naming after the processing. \
  Each process needs two specifications:
    parameter | type | description
    ---|---|---
    `tau_gen_modes` | `list` | split of the events corresponding to the origin of the hadronic tau
    `samples` | `list` | list of all sample tags corresponding to the specific process
  
  The `tau_gen_modes` have following modes:
    parameter | type | description
    ---|---|---
    `T` | `string` | genuine tau
    `J` | `string` | jet misidentified as a tau
    `L` | `string` | lepton misidentified as a tau
    `all` | `string` | if no split should be performed

* In `event_selection`, parameter for all selections that should be applied are defined. \
  This is basically a dictionary of cuts where the key is the name of a cut and the value is the cut itself as a string e.g. `had_tau_pt: "pt_2 > 30"`. The name of a cut is not really important, it is only used as an output information in the terminal. A cut can only use variables which are in the ntuples.

* In `mc_weights` all weights that should be applied for simulated samples are defined. \
  There are two types of weights.
  1. Like for `event_selection` a weight can directly be specified and is then applied to all samples the same way e.g. `lep_id: "id_wgt_mu_1"`
  2. Some weights are either sample specific or need additional information. Currently implemented options are:
      parameter | type | description
      ---|---|---
      `generator` | `string` | `""` if a normal generator weight should be applied to all samples, if `"stitching"` for DY+jets and W+jets a special stitching weights is applied
      `lumi` | `string` | luminosity scaling, this depends on the era and uses the `era` parameter of the config to get the correct weight, so basically it's not relevant what is in the string
      `Z_pt_reweight` | `string` | reweighting of the Z boson pt, the weight in the ntuple is used and only applied to DY+jets
      `Top_pt_reweight` | `string` | reweighting of the top quark pt, the weight in the ntuple is used and only applied to ttbar

* In `emb_weights` all weights that should be applied for embedded samples are defined. \
  Like for `event_selection` a weight can directly be specified and is then applied to all samples the same way e.g. `single_trigger: "trg_wgt_single_mu24ormu27"`

Scale factors for b-tagging and tau ID vs jet are applied on the fly during the FF calculation step. 

To run the preselection step, execute the python script and specify the config file (relative path possible):
```bash
python preselection.py --config-file PATH/CONFIG.yaml
```

</details>

<details>
<summary>

## Fake Factor calculation

</summary>

In this step the fake factors are calculated. This should be run after the preselection step.

All information for the FF calculation step is defined in a configuration file in the `configs/` folder. \
The FF calculation config has the following parameters:

* The expected input folder structure is FILE_PATH/preselection/ERA/CHANNEL/*.root
    parameter | type | description
    ---|---|---
    `file_path` | `string` | absolute path to the folder with the preselected files
    `era` | `string` | data taking era ("2018, "2017", "2016preVFP", "2016postVFP")
    `channel` | `string` | tau pair decay channels ("et", "mt", "tt")
    `tree` | `string` | name of the tree in the preselected files (same as in preselection e.g. "ntuple")

* The output folder structure is workdir/WORKDIR_NAME/ERA/fake_factors/CHANNEL/*outputfiles*
    parameter | type | description
    ---|---|---
    `workdir_name` | `string` | relative path where the output files will be stored

* General options for the calculation:
    parameter | type | description
    ---|---|---
    `use_embedding` | `bool` | True if embedded sample should be used, False if only MC sample should be used

* In `target_processes` the processes for which FFs should be calculated (normally for QCD, Wjets, ttbar) are defined. \
  Each target process needs some specifications:
    parameter | type | description
    ---|---|---
    `split_categories` | `dict` | names of variables for the fake factor measurement in different phase space regions <ul><li>the FF measurement can be split based on variables in 1D or 2D (1 or 2 variables)</li><li>each category/variable has a `list` of orthogonal cuts (e.g. "njets" with "==1", ">=2")</li><li>implemented split variables are "njets", "nbtag" or "deltaR_ditaupair"</li><li>at least one inclusive category needs to be specified</li></ul>
    `split_categories_binedges` | `dict` | bin edge values for each `split_categories` variable <ul><li>number of bin edges should always be N(variable cuts)+1</li></ul>
    `SRlike_cuts` | `dict` | event selections for the signal-like region of the target process
    `ARlike_cuts` | `dict` | event selections for the application-like region of the target process
    `SR_cuts` | `dict` | event selections for the signal region (normally only needed for ttbar)
    `AR_cuts` | `dict` | event selections for the application region (normally only needed for ttbar)
    `var_dependence` | `string` | variable the FF measurement should depend on (normally pt of the hadronic tau e.g. `"pt_2"`)
    `var_bins` | `list` | bin edges for the variable specified in `var_dependence`
    
    Event selections can be defined the same way as in the preselection step `event_selection`. Only the tau vs jet ID cut is special because the name should always be `had_tau_id_vs_jet` (or `had_tau_id_vs_jet_*` in tt channel), this is needed to read out the working points from the cut string and apply the correct tau vs jet ID weights.

* In `process_fractions` specifications for the calculation of the process fractions are defined.
    parameter | type | description
    ---|---|---
    `processes` | `list` | sample names (from the preprocessing step) of the processes for which the fractions should be stored in the correctionlib json, the sum of fractions of the specified samples is 1.
    `split_categories` | `dict` | see `target_processes` (only in 1D)
    `AR_cuts` | `list` | see `target_processes`
    `SR_cuts` | `list` | see `target_processes`, (optional) not needed for the fraction calculation
  
To run the FF calculation step, execute the python script and specify the config file (relative path possible):
```bash
python ff_calculation.py --config-file PATH/CONFIG.yaml
```

</details>

<details>
<summary>

## Fake Factor corrections

</summary>

In this step the corrections for the fake factors are calculated. This should be run after the FF calculation step.

Currently two different correction types are implemented: 
1. non closure correction depending on a specific variable
2. DR to SR interpolation correction depending on a specific variable

All information for the FF correction calculation step is defined in a configuration file in the `configs/` folder. Additional information is loaded from the used config in the previous FF calculation step (this is done automatically). \
The FF correction config has the following parameters:

* The expected input folder structure is workdir/WORKDIR_NAME/ERA/fake_factors/CHANNEL/*
    parameter | type | description
    ---|---|---
    `workdir_name` | `string` | the name of the work directory for which the corrections should be calculated (normally the same as in the FF calculation step)
    `era` | `string` | data taking era ("2018, "2017", "2016preVFP", "2016postVFP")
    `channel` | `string` | tau pair decay channels ("et", "mt", "tt")

* In `target_processes` the processes for which FF corrections should be calculated (normally for QCD, Wjets, ttbar) are defined. \
  Each target process needs some specifications:
    parameter | type | description
    ---|---|---
    `non_closure` | `dict` | one or two non closure corrections can be specified indicated by the variable the correction should be calculated for (e.g. `leading_lep_pt`), if more than one correction is specified, `leading_lep_pt` should come first (due to code specifics) because the second corrections is calculated with the first already applied
    `DR_SR` | `dict` | this correction should be specified only once per process in `target_processes`

  Each correction has following specifications:
    parameter | type | description
    ---|---|---
    `var_dependence` | `string` | variable the FF correction measurement should depend on (e.g. `"pt_1"` for "leading_lep_pt")
    `var_bins` | `list` | bin edges for the variable specified in `var_dependence`
    `SRlike_cuts` | `dict` | event selections for the signal-like region of the target process that should be replaced compared to the selection used in the previous FF calculation step
    `ARlike_cuts` | `dict` | event selections for the application-like region of the target process that should be replaced compared to the selection used in the previous FF calculation step
    `AR_SR_cuts` | `dict` | event selections for a switch from the determination region to the signal/application region, this is only relevant for `DR_SR` corrections
    `non_closure` | `dict` | this is only relevant for `DR_SR` corrections, since for this corrections additional fake factors are calculated it's possible to calculated and apply non closure corrections to these fake factors before calculating the actual DR to SR correction
    
To run the FF correction step, execute the python script and specify the config file (relative path possible):
```bash
python ff_corrections.py --config-file PATH/CONFIG.yaml 
```
An optional parameter is `--only-main-corrections`. By using this parameter the precalculation step for the DR to SR corrections is skipped. This is helpful is the precalculations step is already done.
</details>

## Hints

* check out `configs/general_definitions.py`, this file has many relevant definition for preselection (which variables to save), plotting (dictionaries for names) or correctionlib output information
* check `ntuple_path` and `output_path` (preselection) or `file_path` and `workdir_name` (fake factors, corrections) in the used config files to avoid wrong inputs or outputs 
