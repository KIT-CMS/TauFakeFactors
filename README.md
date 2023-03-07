# TauFakeFactors
FakeFactor framework for estimation of misidentified tau's with ROOT.

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
The preselection config has the following options:

* The expected input folder structure is NTUPLE_PATH/ERA/SAMPLE_TAG/CHANNEL/*.root
    parameter | type | description
    ---|---|---
    `ntuple_path` | `string` | absolute path to the folder with the n-tuples or remote via xrootd
    `era` | `string` | data taking era ("2018, "2017", "2016preVFP", "2016postVFP")
    `channel` | `string` | tau pair decay channels ("et", "mt", "tt")
    `tree` | `string` | name of the tree in the n-tuple files (in CROWN "ntuple")

* The output folder structure is OUTPUT_PATH/preselection/ERA/CHANNEL/*.root
    parameter | type | description
    ---|---|---
    `output_path` | `string` | absolute path where the files with the preselected events will be stored

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
  Currently implemented options are:
    parameter | type | description
    ---|---|---
    `had_tau_pt` | `string` | threshold for the transverse momentum of the hadronic tau in GeV (e.g. ">30")
    `had_tau_eta` | `string` | threshold for the pseudo rapidity of the hadronic tau as absolute value (e.g. "<2.3")
    `had_tau_decay_mode` | `list` | of all tau decay modes to consider (e.g. ["0","1"])
    `had_tau_id_vs_ele` | `string` | working point for the tau ID vs electron (e.g. "Tight")
    `had_tau_id_vs_mu` | `string` | working point for the tau ID vs muon (e.g. "VLoose")
    `lep_iso` | `string` | threshold for the lepton (e/mu) isolation (e.g. "<0.15")
    `trigger` | `bool` | True if a trigger selection should be applied, False otherwise

* In `mc_weights` all weights that should be applied for simulated samples are defined. \
  Currently implemented options are:
    parameter | type | description
    ---|---|---
    `generator` | `string` | generator weight from MC production
    `lumi` | `string` | luminosity scaling
    `pileup` | `string` | pileup weight
    `lep_iso` | `string` | lepton (e/mu) isolation scale factor
    `lep_id` | `string` | lepton (e/mu) identification scale factor
    `had_tau_id_vs_ele` | `string` | tau ID vs electron scale factor for the working point chosen in the `event_selection`
    `had_tau_id_vs_mu` | `string` | tau ID vs muon scale factor for the working point chosen in the `event_selection`
    `trigger` | `string` | trigger scale factor
    `Z_pt_reweight` | `string` | reweighting of the Z boson pt
    `Top_pt_reweight` | `string` | reweighting of the top quark pt

* In `emb_weights` all weights that should be applied for embedded samples are defined. \
  Currently implemented options are:
    parameter | type | description
    ---|---|---
    `generator` | `string` | generator weight from MC production
    `lep_iso` | `string` | lepton (e/mu) isolation scale factor
    `lep_id` | `string` | lepton (e/mu) identification scale factor
    `trigger` | `string` | trigger scale factor

Scale factors for b-tagging and tau ID vs jet are applied on the fly during the FF calculation step. 

To run the preselection step execute the python script and specify the config file name:
```bash
python preselection.py --config CONFIG_NAME 
```

</details>

<details>
<summary>

## Fake Factor calculation

</summary>

In this step the fake factors are calculated. This should be run after the preselection step.

All information for the FF calculation step is defined in a configuration file in the `configs/` folder. \
The FF calculation config has the following options:

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
    `generate_json` | `bool` | True if a correctionlib json file with the FFs should be produced, False otherwise

* In `target_process` the processes for which FFs should be calculated (normally for QCD, Wjets, ttbar) are defined. \
  Each target process needs some specifications:
    parameter | type | description
    ---|---|---
    `split_categories` | `list` | names of variables for the fake factor measurement in different phase space regions <ul><li>the FF measurement can be split based on variables in 1D or 2D (1 or 2 variables)</li><li>each category/variable has a `list` of orthogonal cuts (e.g. "njets" with "==1", ">=2")</li><li>implemented split variables are "njets", "nbtag" or "deltaR_ditaupair"</li><li>at least one inclusive category needs to be specified</li></ul>
    `split_categories_binedges` | `list` | bin edge values for each `split_categories` variable <ul><li>number of bin edges should always be N(variable cuts)+1</li><li>is only used if `generate_json` is True</li></ul>
    `SRlike_cuts` | `list` | event selections for the signal-like region of the target process
    `ARlike_cuts` | `list` | event selections for the application-like region of the target process
    `SR_cuts` | `list` | event selections for the signal region (normally only needed for ttbar)
    `AR_cuts` | `list` | event selections for the application region (normally only needed for ttbar)
    `var_dependence` | `string` | variable the FF measurement should depend on (normally pt of the hadronic tau e.g. "pt_2")
    `var_bins` | `list` | bin edges for the variable specified in `var_dependence`
    
    Implemented event selection cuts are (besides the already mentioned cuts in the preselection step):
    parameter | type | description
    ---|---|---
    `tau_pair_sign` | `string` | two options "same" or "opposite"
    `nbtag` | `string` | number of b-tagged jets (e.g. ">=1")
    `lep_mt` | `string` | threshold for the transverse mass of the lepton (e/mu) + MET pair in GeV (e.g. "<50")
    `no_extra_lep` | `bool` | True if no other leptons than the tau pair are allowed, False if other leptons should be present 
    `had_tau_id_vs_jet` | `string`/`list` | select events above a working point (e.g. "Tight") / select events between two working points (e.g. ["VVVLoose","Tight"])

* In `process_fractions` specifications for the calculation of the process fractions are defined.
    parameter | type | description
    ---|---|---
    `processes` | `list` | sample names (from the preprocessing step) of the processes for which the fractions should be stored in the correctionlib json
    `split_categories` | `list` | see `target_process` (only in 1D)
    `AR_cuts` | `list` | see `target_process`
    `SR_cuts` | `list` | see `target_process`, (optional) not needed for the fraction calculation
  
To run the FF calculation step execute the python script and specify the config file name:
```bash
python ff_calculation.py --config CONFIG_NAME 
```

</details>

<details>
<summary>

## Fake Factor corrections

</summary>

In this step the corrections for the fake factors are calculated. This should be run after the FF calculation step.

Currently two different correction types are implemented: 
1. non closure correction dependent of a specific variable
2. DR to SR interpolation correction dependent of a specific variable

All information for the FF correction calculation step is defined in a configuration file in the `configs/` folder. Additional information is loaded from the used config in the previous FF calculation step. \
The FF correction config has the following options:

* The expected input folder structure is workdir/WORKDIR_NAME/ERA/fake_factors/CHANNEL/*
    parameter | type | description
    ---|---|---
    `workdir_name` | `string` | the name of the work directory for which the corrections should be calculated (normally the same as in the FF calculation step)
    `era` | `string` | data taking era ("2018, "2017", "2016preVFP", "2016postVFP")
    `channel` | `string` | tau pair decay channels ("et", "mt", "tt")

* General options for the calculation:
    parameter | type | description
    ---|---|---
    `generate_json` | `bool` | True if a correctionlib json file with the FF corrections should be produced, False otherwise

* In `target_process` the processes for which FF corrections should be calculated (normally for QCD, Wjets, ttbar) are defined. \
  Each target process needs some specifications:
    parameter | type | description
    ---|---|---
    `non_closure` | `string` | multiple non closure corrections can be specified indicated by the variable to correct on (e.g. "lep_pt")
    `DR_SR` | `string` | this correction can be specified once per `target_process`

  Each correction has following specifications:
    parameter | type | description
    ---|---|---
    `var_dependence` | `string` | variable the FF correction measurement should depend on (e.g. "pt_1" for "lep_pt")
    `var_bins` | `list` | bin edges for the variable specified in `var_dependence`
    `SRlike_cuts` | `list` | event selections for the signal-like region of the target process that should be adjusted compared to the selection used in the previous FF calculation step
    `ARlike_cuts` | `list` | event selections for the application-like region of the target process that should be adjusted compared to the selection used in the previous FF calculation step
    
</details>
