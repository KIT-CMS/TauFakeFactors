# TauFakeFactors
FakeFactor framework for the estimation of jets misidentified taus with pyROOT.

## Setup
Clone the repository via
```bash
git clone --recurse-submodules https://github.com/KIT-CMS/TauFakeFactors.git
```
The environment can be set up with conda via
```bash
conda env create --file environment.yaml
```

General definitions like paths for all steps of the fake factor measurements should be defined in a `configs/ANALYSIS/ERA/common_settings.yaml` file (file name should be always stay the same).

The expected ntuple folder structure is NTUPLE_PATH/ERA/SAMPLE_TAG/CHANNEL/*.root
parameter | type | description
  ---|---|---
  `ntuple_path` | `string` | absolute path to the folder with the n-tuples on the dcache, a remote path is expected like "root://cmsxrootd-kit.gridka.de//store/user/USER/..."
  `tree` | `string` | name of the tree in the n-tuple files ("ntuple" in CROWN)
  `era` | `string` | data taking era (e.g. "2018, "2017", "2016preVFP", "2016postVFP")
  `tau_vs_jet_wps` | `list` | list of tau ID vsJet working points to be written out in the preselection step (e.g. ["Medium", "VVVLoose"])
  `tau_vs_jet_wgt_wps` | `list` | list of tau ID vsJet working point scale factors to be written out in the preselection step (e.g. ["Medium"])

The output folder structure is OUTPUT_PATH/preselection/ERA/CHANNEL/*.root
  parameter | type | description
  ---|---|---
  `output_path` | `string` | absolute path where the files with the preselected events will be stored, a local path is expected like "/ceph/USER/..."
  `file_path` | `string` | absolute path to the folder with the preselected files (should be the same as `output_path`) to be used for the fake factor calculation
  `workdir_name` | `string` | relative path where the fake factor measurement output files will be stored; folder is produced in `workdir/`


<details>
<summary>

## Event preselection

</summary>

This framework is designed for n-tuples produced with CROWN as input. 
All information for the preselection step is defined in configuration files in the `configs/ANALYSIS/ERA/` folder using the `common_settings.yaml` file and a more specific config file. 

The preselection config has the following parameters:

* parameter | type | description
    ---|---|---
    `channel` | `string` | tau pair decay channels ("et", "mt", "tt")

* In `processes` all the processes are defined that should be preprocessed. \
  The names are also used for the output file naming after the processing. \
  Each process needs two specifications:
    subparameter | type | description
    ---|---|---
    `tau_gen_modes` | `list` | split of the events corresponding to the origin of the hadronic tau
    `samples` | `list` | list of all sample tags corresponding to the specific process
  
  The `tau_gen_modes` have following modes:
    subparameter | type | description
    ---|---|---
    `T` | `string` | genuine tau
    `J` | `string` | jet misidentified as a tau
    `L` | `string` | lepton misidentified as a tau
    `all` | `string` | if no split should be performed

* In `event_selection`, parameter for all selections that should be applied are defined. \
  This is basically a dictionary of cuts where the key is the name of a cut and the value is the cut itself as a string e.g. `had_tau_pt: "pt_2 > 30"`. The name of a cut is not really important, it is only used as an output information in the terminal. A cut can only use variables which are in the ntuples.

* In `mc_weights` all weights that should be applied for simulated samples are defined. \
  There are two types of weights.
  1. Similar to `event_selection`, a weight can directly be specified and is then applied to all samples in the same way e.g. `lep_id: "id_wgt_mu_1"`
  2. But some weights are either sample specific or need additional information. Currently implemented options are:
      subparameter | type | description
      ---|---|---
      `generator` | `string` | The normal generator weight is applied to all samples, if they aren't specified in the `"stitching"` sub-group. Stitching weights might be needed for DY+jets or W+jets, depending on which samples are used for them. 
      `lumi` | `string` | luminosity scaling, this depends on the era and uses the `era` parameter of the config to get the correct weight, so basically it's not relevant what is in the string
      `Z_pt_reweight` | `string` | reweighting of the Z boson pt, the weight in the ntuple is used and only applied to DY+jets
      `Top_pt_reweight` | `string` | reweighting of the top quark pt, the weight in the ntuple is used and only applied to ttbar

* In `emb_weights` all weights that should be applied for embedded samples are defined. \
  Like for `event_selection` a weight can directly be specified and is then applied to all samples the same way e.g. `single_trigger: "trg_wgt_single_mu24ormu27"`
* In `output_features` the to be saved/needed features for the later calculations are listed.

Scale factors for b-tagging and tau ID vs jet are applied on the fly during the FF calculation step. 

To run the preselection step, execute the python script and specify the config file (relative path possible):
```bash
python preselection.py --config-file configs/PATH/CONFIG.yaml
```
Further there are additional optional parameters: 
1. `--nthreads=SOME_INTEGER` to define the number of threads for the multiprocessing pool to run the sample processing in parallel. Default value is 8 (this should normally cover running all of the samples in parallel).
2.  `--ncores=SOME_INTEGER` to define the number of cores that should be used for each pool thread to speed up the ROOT dataframe calculation. Default value is 4.

</details>

<details>
<summary>

## Fake Factor calculation

</summary>

In this step the fake factors are calculated. This should be run after the preselection step.

All information for the FF calculation step is defined in a configuration file in the `configs/ANALYSIS/ERA/` folder using the `common_settings.yaml` and a more specific config file. \
The FF calculation config has the following parameters:

* General options for the calculation:
    parameter | type | description
    ---|---|---
    `channel` | `string` | tau pair decay channels ("et", "mt", "tt")
    `use_embedding` | `bool` | True if embedded sample should be used, False if only MC sample should be used
    `use_center_of_mass_bins` | `bool` | Changes the x-data that is entering FF and correction calculation. If set then a center of mass value is used for the x-data, calculated from events entering the corresponding bin. If not set, the bin centers are used. Default is set to True. <br> <br>This will not affect FF and correction calculation that are set to `"binwise"` (the x-data values although displayed in plots are not used)

* In `target_processes` the processes for which FFs should be calculated (normally for QCD, Wjets, ttbar) are defined. \
  Each target process needs some specifications:
    parameter | type | description
    ---|---|---
    `split_categories` | `dict` | names of variables for the fake factor measurement in different phase space regions <ul><li>the FF measurement can be split based on variables in 1D or 2D (1 or 2 variables)</li><li>each category/variable has a `list` of orthogonal cuts (e.g. "njets" with "==1", ">=2")</li><li> "njets", "nbtag", "tau_decaymode_2" or "deltaR_ditaupair" are already possible, other variables should be added during preprocessing step accordingly </li><li>at least one inclusive category needs to be specified (assuming variable is written out in preselection step)</li></ul>
    `split_categories_binedges` | `dict` | bin edge values for each `split_categories` variable <ul><li>number of bin edges should always be N(variable cuts)+1</li></ul>
    `SRlike_cuts` | `dict` | event selections for the signal-like region of the target process
    `ARlike_cuts` | `dict` | event selections for the application-like region of the target process
    `SR_cuts` | `dict` | event selections for the signal region (normally only needed for ttbar)
    `AR_cuts` | `dict` | event selections for the application region (normally only needed for ttbar)
    `var_dependence` | `string` | variable the FF measurement should depend on (normally pt of the hadronic tau e.g. `"pt_2"`)
    `var_bins` | `list` or `dict[list]` | bin edges for the variable specified in `var_dependence`. <br><br> Can either be a list representing the binning or a dictionary of lists, where keys correspond to the string representations of split categories defined in `split_categories`. <br><br> In the case of two split categories, the dictionary can be nested. If not all second split category elements share the first split category binning, the binning for the affected category must be specified separately. When using split binning, at least the first split category's bin edges must be fully defined.
    `fit_options` | `list` | a list of polynomials that should be considered for the fake factor fits can be defined with this parameter (default: `["poly_1"]`); futher it is possible to specify `"binwise"` which means that the histograms are written out directly without a fit 
    `limit_kwargs` | `dict` | this dictionary allows to define how the fitted function and its uncertainty are handled (also outside the measurement range); the default is that outside the range the fake factor functions stays constant but the up and down variations still increase/decrease; additionally negative fake factor values are not allowed and if present are set to 0
     
    Event selections can be defined the same way as in the preselection step `event_selection`. Only the tau vs jet ID cut is special because the name should always be `had_tau_id_vs_jet` (or `had_tau_id_vs_jet_*` in tt channel), this is needed to read out the working points from the cut string and apply the correct tau vs jet ID weights.

* In `process_fractions` specifications for the calculation of the process fractions are defined.
    parameter | type | description
    ---|---|---
    `processes` | `list` | sample names (from the preprocessing step) of the processes for which the fractions should be stored in the correctionlib json, the sum of fractions of the specified samples is 1.
    `split_categories` | `dict` | see `target_processes` (only in 1D)
    `AR_cuts` | `list` | see `target_processes`
    `SR_cuts` | `list` | see `target_processes`, (optional) not needed for the fraction calculation

    **Note:** When using split binning for process fraction calculations, the `var_bins` parameter can also be defined in the same manner as for `target_processes`.
  
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

All information for the FF correction calculation step is defined in a configuration file in the `configs/ANALYSIS/ERA/` folder using the `common_settings.yaml` and a more specific config file. Additional information is loaded from the used config in the previous FF calculation step (this is done automatically). \
The FF correction config has the following parameters:

* The expected input folder structure is workdir/WORKDIR_NAME/ERA/fake_factors/CHANNEL/*
    parameter | type | description
    ---|---|---
    `channel` | `string` | tau pair decay channels ("et", "mt", "tt")

* In `target_processes` the processes for which FF corrections should be calculated (normally for QCD, Wjets, ttbar) are defined.
* `split_categories` can be set for `non_closure` and `DR_SR` corrections (1D only) accordingly. `var_bins` declaration follow the specifications named in `FakeFactor calculation` section \
  Each target process needs some specifications:
    parameter | type | description
    ---|---|---
    `non_closure` | `dict` | one or two non closure corrections can be specified indicated by the variable the correction should be calculated for (e.g. `pt_1`), also more than one closure correction are allowed and are calculated while already applying the correction that were already measured
    `DR_SR` | `dict` | this correction should be specified only once per process in `target_processes`

  Each correction has following specifications:
    parameter | type | description
    ---|---|---
    `var_dependence` | `string` | variable the FF correction measurement should depend on (e.g. `"pt_1"`)
    `split_categories` | `dict` | Optional, analogous to `FakeFactor calculation` (only 1D).
    `var_bins` | `list` or `dict[list]` | Analogous to `var_bins` in `FakeFactor calculation`
    `write_corrections` | `str` | "smoothed" or "binwise"; Definition of written out correction. "smoothed" applies a gaussian density kernel.
    `SRlike_cuts` | `dict` | event selections for the signal-like region of the target process that should be replaced compared to the selection used in the previous FF calculation step
    `ARlike_cuts` | `dict` | event selections for the application-like region of the target process that should be replaced compared to the selection used in the previous FF calculation step
    `AR_SR_cuts` | `dict` | event selections for a switch from the determination region to the signal/application region, this is only relevant for `DR_SR` corrections
    `non_closure` | `dict` | this is only relevant for `DR_SR` corrections, since for this corrections additional fake factors are calculated. It's possible to calculated and apply non closure corrections to these fake factors before calculating the actual DR to SR correction.
    
To run the FF correction step, execute the python script and specify the config file (relative path possible):
```bash
python ff_corrections.py --config-file PATH/CONFIG.yaml 
```
There are two optional parameters `--skip-DRtoSR-ffs` and `--only-main-corrections`. The correction caclulation is done in 3 steps.\
The first step is to calculate additional fake factors which are needed for the final DR to SR correction. If this is already done, this step can be skipped using `--skip-DRtoSR-ffs`.\
The second step is to calculate non closure corrections for these additional DR to SR fake factors. If both steps are already done they can be skipped by using `--only-main-corrections`.\
The last step is to calculate all the specified corrections for the main fake factors.
</details>

## Hints

* check out `configs/general_definitions.py`, this file has many relevant definition for plotting (dictionaries for names) and correctionlib output information
* check `ntuple_path` and `output_path` (preselection) and `file_path` and `workdir_name` (fake factors, corrections) in the used config files to avoid wrong inputs or outputs
