# TauFakeFactors
FakeFactor framework for estimation of misidentified tau's with ROOT.

## Setup
The environment can be set up with conda via
```bash
conda env create --file environment.yaml
```

## Event preselection
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

## Fake Factor calculation
In this step the fake factors are calculated. This should be run after the preselection step.

All information for the FF calculation step is defined in a configuration file in the `configs` folder. 
The FF calculation config has the following options:

* the expected input folder structure is FILE_PATH//preselection/ERA/CHANNEL/*.root
    * `file_path`: `string` absolute path to the folder with the preselected files
    * `era`: `string` data taking era ("2018, "2017", "2016preVFP", "2016postVFP")
    * `channel`: `string` tau pair decay channels ("et", "mt", "tt")
    * `tree`: `string` name of the tree in the preselected files (same as in preselection e.g. "ntuple")
* the output folder structure is workdir/WORKDIR_NAME/CHANNEL/*outputfiles*
    * `workdir_name`: `string` relative path where the output files will be stored
* `use_embedding`: `bool` True if embedded sample should be used, False if only MC sample should be used
* `generate_json`: `bool` True if a correctionlib json file with the FFs should be produced, False otherwise
* `target_process`: `list` of processes for which FFs should ve calculated (normally for QCD, Wjets, ttbar) \
  each target process needs specifications:
    * `split_categories`: `list` of variables 
        * the FF measurement can be split based on variables in 1D or 2D (1 or 2 variables)
        * each category/variable has a `list` of orthogonal cuts (e.g. "njets" with "==1", ">=2") 
        * implemented split variables are "njets" and "nbtag"
        * at least one inclusive category needs to be specified
    * `split_categories_binedges`: `list` of bin edge values for each variable 
      * variables should be same the as in `split_categories`
      * number of bin edges should always be N(variable cuts)+1
      * is only used if `generate_json` is True
    * `SRlike_cuts`: `list` of event selections specific for the signal-like region of the target process
    * `ARlike_cuts`: `list` of event selections specific for the application-like region of the target process
    * `SR_cuts`: `list` of event selections specific for the signal region (only needed for ttbar)
    * `AR_cuts`: `list` of event selections specific for the application region (only needed for ttbar) \
      implemented event selections are:
        * `tau_pair_sign`: `string` two options "same" or "opposite"
        * `nbtag`: `string` number of b-tagged jets (e.g. ">=1")
        * `lep_mt`: `string` threshold for the transverse mass of the lepton (e/mu) + MET pair in GeV (e.g. "<50")
        * `no_extra_lep`: `bool` True if no other leptons than the tau pair are allowed, False if other leptons should be present 
        * `had_tau_id_vs_jet`: `string` to select events above a working point (e.g. "Tight") \
        or `list` to select events between two working points (e.g. ["VVVLoose","Tight"])
    * `var_dependence`: `string` variable the FF measurement should depend on (normally pt of the hadronic tau e.g. "pt_2")
    * `var_bins`: `list` of bin edges for the variable specified in `var_dependence`
* `process_fractions`: `list` of specifications for the calculation of the process fractions
    * `processes`: `list` of sample names (processes) for which the fractions should be stored in the correctionlib json
    * `split_categories`: see `target_process` (only in 1D)
    * `AR_cuts`: see `target_process`
    * `SR_cuts`: see `target_process`, (optional) not needed for the fraction calculation
  

To run the FF calculation step execute the python script and specify the config file name:
```bash
python ff_calculation.py --config CONFIG_NAME 
```

## Fake Factor validation
