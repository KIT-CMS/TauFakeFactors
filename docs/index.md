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

## General configuration
General definitions like paths for all steps of the fake factor measurements should be defined in a `configs/ANALYSIS/ERA/common_settings.yaml` file (file name should be always stay the same).

The expected ntuple folder structure is NTUPLE_PATH/ERA/SAMPLE_TAG/CHANNEL/*.root

  parameter | type | description
  ---|---|---
  `ntuple_path` | `string` | absolute path to the folder with the n-tuples on the dcache, a remote path is expected like "root://cmsdcache-kit-disk.gridka.de//store/user/USER/..."
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

## Event preselection
This framework is designed for n-tuples produced with CROWN as input. 
All information for the preselection step is defined in configuration files in the `configs/ANALYSIS/ERA/` folder using the `common_settings.yaml` file and a more specific config file. 

The preselection config has the following parameters:

  parameter | type | description
  ---|---|---
  `channel` | `string` | tau pair decay channels ("et", "mt", "tt")
  `processes` | `dict` | process parameters are explained below
  `event_selection` | `dict` | with this parameter all selections that should be applied are defined. <br>This is basically a dictionary of cuts where the key is the name of a cut and the value is the cut itself as a string e.g. `had_tau_pt: "pt_2 > 30"`. The name of a cut is not really important, it is only used as an output information in the terminal. A cut can only use variables which are in the ntuples.
  `mc_weights` | `dict` | weight parameter are defined below
  `emb_weights` | `dict` | all weights that should be applied for embedded samples are defined. <br>Like for `event_selection` a weight can directly be specified and is then applied to all samples the same way e.g. `single_trigger: "trg_wgt_single_mu24ormu27"`
  `output_features` | `list` | the to be saved/needed features for the later calculations are listed 

In `processes` all the processes are defined that should be preprocessed. <br>
The names are also used for the output file naming after the processing. <br>
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

In `mc_weights` all weights that should be applied for simulated samples are defined. <br>
There are two types of weights.

1. Similar to `event_selection`, a weight can directly be specified and is then applied to all samples in the same way e.g. `lep_id: "id_wgt_mu_1"`
2. But some weights are either sample specific or need additional information. Currently implemented options are:

    subparameter | type | description
    ---|---|---
    `generator` | `string` | The normal generator weight is applied to all samples, if they aren't specified in the `"stitching"` sub-group. Stitching weights might be needed for DY+jets or W+jets, depending on which samples are used for them. 
    `lumi` | `string` | luminosity scaling, this depends on the era and uses the `era` parameter of the config to get the correct weight, so basically it's not relevant what is in the string
    `Z_pt_reweight` | `string` | reweighting of the Z boson pt, the weight in the ntuple is used and only applied to DY+jets
    `Top_pt_reweight` | `string` | reweighting of the top quark pt, the weight in the ntuple is used and only applied to ttbar

Scale factors for b-tagging and tau ID vs jet are applied on the fly during the FF calculation step. 

To run the preselection step, execute the python script and specify the config file (relative path possible):
```bash
python preselection.py --config-file configs/PATH/CONFIG.yaml
```
Further there are additional optional parameters: 

1. `--nthreads=SOME_INTEGER` to define the number of threads for the multiprocessing pool to run the sample processing in parallel. Default value is 8 (this should normally cover running all of the samples in parallel).
2.  `--ncores=SOME_INTEGER` to define the number of cores that should be used for each pool thread to speed up the ROOT dataframe calculation. Default value is 2.

## Fake Factor calculation
In this step the fake factors are calculated. This should be run after the preselection step.

All information for the FF calculation step is defined in a configuration file in the `configs/ANALYSIS/ERA/` folder using the `common_settings.yaml` and a more specific config file. <br>
The FF calculation config has the following parameters:

General options for the calculation:

  parameter | type | description
  ---|---|---
  `channel` | `string` | tau pair decay channels ("et", "mt", "tt")
  `use_embedding` | `bool` | True if embedded sample should be used, False if only MC sample should be used
  `use_center_of_mass_bins` | `bool` | Changes the x-data that is entering FF and correction calculation. If set then a center of mass value is used for the x-data, calculated from events entering the corresponding bin. If not set, the bin centers are used. Default is set to True. <br> <br>This will not affect FF and correction calculation that are set to `"binwise"` (the x-data values although displayed in plots are not used)

In `target_processes` the processes for which FFs should be calculated (normally for QCD, Wjets, ttbar) are defined. <br>
Each target process needs some specifications:

  parameter | type | description
  ---|---|---
  `split_categories` | `dict` | names of variables for the fake factor measurement in different phase space regions <ul><li>the FF measurement can be split based on variables in 1D or 2D (1 or 2 variables)</li><li>each category/variable has a `list` of orthogonal cuts (e.g. "njets" with "==1", ">=2")</li><li> "njets", "nbtag", "tau_decaymode_2" or "deltaR_ditaupair" are already possible, other variables should be added during preprocessing step accordingly </li><li>at least one inclusive category needs to be specified (assuming variable is written out in preselection step)</li></ul> If a continous variable is used a window can be defined as `">=lower#&&#<upper"` accordingly.
  `split_categories_binedges` | `dict` | bin edge values for each `split_categories` variable. <br> The number of bin edges should always be N(variable cuts)+1
  `SRlike_cuts` | `dict` | event selections for the signal-like region of the target process
  `ARlike_cuts` | `dict` | event selections for the application-like region of the target process
  `SR_cuts` | `dict` | event selections for the signal region (normally only needed for ttbar)
  `AR_cuts` | `dict` | event selections for the application region (normally only needed for ttbar)
  `var_dependence` | `string` | variable the FF measurement should depend on (normally pt of the hadronic tau e.g. `"pt_2"`)
  `var_bins` | `list` or `dict[list]` | bin edges for the variable specified in `var_dependence`. <br> Can either be a list representing the binning or a dictionary of lists, where keys correspond to the string representations of split categories defined in `split_categories`. <br> In the case of two split categories, the dictionary can be nested. If not all second split category elements share the first split category binning, the binning for the affected category must be specified separately. When using split binning, at least the first split category's bin edges must be fully defined.
  `fit_options` | `list` | a list of polynomials that should be considered for the fake factor fits can be defined with this parameter (default: `["poly_1"]`); futher it is possible to specify `"binwise"` which means that the histograms are written out directly without a fit 
  `limit_kwargs` | `dict` | this dictionary allows to define how the fitted function and its uncertainty are handled (also outside the measurement range); the default is that outside the range the fake factor functions stays constant but the up and down variations still increase/decrease; additionally negative fake factor values are not allowed and if present are set to 0
     
Event selections can be defined the same way as in the preselection step `event_selection`. Only the tau vs jet ID cut is special because the name should always be `had_tau_id_vs_jet` (or `had_tau_id_vs_jet_*` in tt channel), this is needed to read out the working points from the cut string and apply the correct tau vs jet ID weights.

In `process_fractions` specifications for the calculation of the process fractions are defined.

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

## Fake Factor corrections
In this step the corrections for the fake factors are calculated. This should be run after the FF calculation step.

Currently two different correction types are implemented: 

1. non closure correction depending on a specific variable
2. DR to SR interpolation correction depending on a specific variable

All information for the FF correction calculation step is defined in a configuration file in the `configs/ANALYSIS/ERA/` folder using the `common_settings.yaml` and a more specific config file. Additional information is loaded from the used config in the previous FF calculation step (this is done automatically). <br>
The FF correction config has the following parameters:

The expected input folder structure is workdir/WORKDIR_NAME/ERA/fake_factors/CHANNEL/*

  parameter | type | description
  ---|---|---
  `channel` | `string` | tau pair decay channels ("et", "mt", "tt")

In `target_processes` the processes for which FF corrections should be calculated (normally for QCD, Wjets, ttbar) are defined.

`split_categories` can be set for `non_closure` and `DR_SR` corrections (1D only) accordingly. `var_bins` declaration follow the specifications named in `FakeFactor calculation` section <br>
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
  `correction_option` | `str` | Definition of how a correction should be measured. Options are `"smoothed"` (applies a gaussian density kernel), `"binwise"` or both can also be combined into e.g. `"binwise#[0,]+smoothed"`, `"binwise#[-1,-2,]+smoothed"` or `"binwise#[0,]#[-1,]+smoothed"` where the nth-bin(s) is computed using binwise method and the rest in smoothed manner. For left bins, positive bin index is used. For last bins, the count is done using negative integers, starting at -1 (like in the example). If both options are provided, then the binwise calculation is applied for the stated left and right bins.
  `bandwidth` | `float` | if `correction_option` includes `"smoothed"` this value can be set to adjust for the bandwidth used during smoothing procedure (has no effect on the result in case of `"binwise"`). If not set the default value of histogram range divided by 5 is chosen. <br> Can either be a float value representing the bandwidth for all corrections or a dictionary of float values, where keys correspond to the string representations of split categories defined in `split_categories`. 
  `SRlike_cuts` | `dict` | event selections for the signal-like region of the target process that should be replaced compared to the selection used in the previous FF calculation step
  `ARlike_cuts` | `dict` | event selections for the application-like region of the target process that should be replaced compared to the selection used in the previous FF calculation step
  `AR_SR_cuts` | `dict` | event selections for a switch from the determination region to the signal/application region, this is only relevant for `DR_SR` corrections
  `non_closure` | `dict` | this is only relevant for `DR_SR` corrections, since for this corrections additional fake factors are calculated. It's possible to calculated and apply non closure corrections to these fake factors before calculating the actual DR to SR correction.
    
To run the FF correction step, execute the python script and specify the config file (relative path possible):
```bash
python ff_corrections.py --config-file PATH/CONFIG.yaml 
```
There are two optional parameters `--skip-DRtoSR-ffs` and `--only-main-corrections`. The correction caclulation is done in 3 steps. <br>
The first step is to calculate additional fake factors which are needed for the final DR to SR correction. If this is already done, this step can be skipped using `--skip-DRtoSR-ffs`. <br>
The second step is to calculate non closure corrections for these additional DR to SR fake factors. If both steps are already done they can be skipped by using `--only-main-corrections`. <br>
The last step is to calculate all the specified corrections for the main fake factors.

## Hints
* check out `configs/general_definitions.py`, this file has many relevant definition for plotting (dictionaries for names) and correctionlib output information
* check `ntuple_path` and `output_path` (preselection) and `file_path` and `workdir_name` (fake factors, corrections) in the used config files to avoid wrong inputs or outputs

## Automated Equipopulated Binning

The framework includes a utility script, `adjust_binning.py`, for an automatic calculation and update binning definitions within the configuration in case of equipopulated binning approach of specified variable (`var_dependence`) across various categories for a more robust statistical treatment in each analysis region. The general intention is to have an equipopulated data distribution in the SRlike region given possible additional splittings (i.e. in `njets`).

The script is intended to be run berofe the fake factor and correction calculation by using

```bash
python adjust_binning.py --config config.yaml --cut-config onfig_with_region_cuts.yaml --processes QCD Wjets --cut-region ARlike
```

| Argument | Type | Description |
|---|---|---|
| `--config` | `string` | Path to the YAML configuration file that will be modified. |
| `--cuts-config` | `string` | (Optional) Path to a separate YAML file containing cut definitions. This is useful for correction YAML files where cuts are not specified directly for each correction but then are sourced from the cuts-config file, which can just be the fake factor YAML configuration. |
| `--processes` | `list[string]` | A list of processes to adjust (e.g., `QCD`, `Wjets`, `ttbar`, `process_fractions`). The default is set to include all processes. |
| `--cut-region` | `string` | The region cut to be applied when determining event counts (e.g., `SRlike`, `ARlike`). `process_fractions` always use the `AR_cut`. |
| `--dry-run` | `bool` | Preview the binning calculations and changes without modifying the config file. |

### Configuration

To enable this feature, you must add a `equipopulated_binning_options` block to the relevant process or correction section in your configuration file. The script uses this block to generate the new binning and will populate the `var_bins` key and update `split_categories` and `split_categories_binedges` if they are left blank i.e. for continuous variables such as `pt_1`.

For example, a manually binned configuration for the `QCD` process:

```yaml
# Before: Manual binning for QCD process
target_processes:
  QCD:
    split_categories:
      njets: ["==0", "==1", ">=2"]
    split_categories_binedges:
      njets: [-0.5, 0.5, 1.5, 12.5]
    var_dependence: pt_2
    var_bins: 
      "==0": [30.0, 31.6, 33.5, 35.9, 38.9, 43.6, 51.6, 150.0]
      "==1": [30.0, 32.2, 35.3, 40.3, 50.9, 150.0]
      ">=2": [30.0, 32.7, 36.5, 42.5, 57.0, 150.0]
```

Can be replaced with a configuration for automated binning:

```yaml
# After: Configuration for automated binning
target_processes:
  QCD:
    split_categories:
      njets: ["==0", "==1", ">=2"]
    split_categories_binedges:
      njets: [-0.5, 0.5, 1.5, 12.5]
    var_dependence: pt_2
    var_bins: {} # This will be filled by the script
    equipopulated_binning_options:
      variable_config:
        pt_2:
          min: 30.0
          max: 150.0
          rounding: 2
      n_bins:  # this part is necessarily only needed for continuous variables
        njets: 3  # number of bins for njets ==0, ==1, >=2
      var_dependence_n_bins: [7, 5, 5]
```

The `equipopulated_binning_options` block has the following structure:

| Parameter | Type | Description |
|---|---|---|
| `variable_config` | `dict` | Defines parameters for variables used in binning. See sub-parameters below. |
| `n_bins` | `dict` | Specifies the desired number of bins per variable in `split_categories` if those are not set there in case of continuous variables, i.e. four equidistant bins in `pt_1`: `pt_1: 4` |
| `var_dependence_n_bins` |  `int` or `list[int]` or `list[list[int]]` | Number of bins used for the variable defined in `var_dependence`, either as int (used by all categories) or a (nested) list of integers defining number of bins per (nested) category created.|

The `variable_config` for each variable contains:

| Sub-parameter | Type | Description |
|---|---|---|
| `min` | `float` | The lower bound for the variable's range. The first bin edge will be set to this value. This is applied as a cut before the calculation of equipopulated binning starts.  |
| `max` | `float` | The upper bound for the variable's range. The last bin edge will be set to this value. This is applied as a cut before the calculation of equipopulated binning starts. |
| `rounding` | `int` | The number of decimal places to which the calculated bin edges will be rounded to and written into the configuration file. |


### Advanced Splitting Strategies

The script supports various one- and two-dimensional splitting strategies, determined by the structure of `split_categories`. The order of variables defines the nesting hierarchy.

1. One-Dimensional Splitting

<ul>
<li>
<details>
<summary><strong>By a discrete variable</strong> (e.g., <code>njets</code>): Categories are explicitly listed.</summary>

Number of bins of `var_dependence` in each category are set explicitly
```yaml
split_categories:
  njets: ["==0", "==1", ">=2"]
split_categories_binedges:
  njets: [-0.5, 0.5, 1.5, 12.5]
equipopulated_binning_options:
  variable_config:
    pt_2:
      min: 30.0
      max: 150.0
      rounding: 4
  var_dependence_n_bins: [7, 5, 5]
```
Same number of bins of `var_dependence` in each category
```yaml
split_categories:
  njets: ["==0", "==1", ">=2"]
split_categories_binedges:
  njets: [-0.5, 0.5, 1.5, 12.5]
equipopulated_binning_options:
  variable_config:
    pt_2:
      min: 30.0
      max: 150.0
      rounding: 4
  var_dependence_n_bins: 7
```
</details>
</li>
<li>
<details>
<summary><strong>By a continuous variable</strong> (e.g., <code>pt_1</code>): An empty list <code>[]</code> is used as a placeholder.</summary>


The script will first create equipopulated bins for `pt_1` and then bin `var_dependence` within each of those new `pt_1` categories.

```yaml
split_categories:
  pt_1: []
split_categories_binedges:
  pt_1: []
equipopulated_binning_options:
  variable_config:
    pt_2: 
      min: 30.0
      max: 150.0
      rounding: 4
    pt_1: 
      min: 0.0
      max: 1000.0
      rounding: 2
  n_bins:
    pt_1: 4 # Creates 4 categories for pt_1, which are then used for pt_2 binning
  var_dependence_n_bins: [7, 7, 6, 5] # number of bins per pt_1 split 
```
</details>
</li>
</ul>


1. Two-Dimensional (Nested) Splitting

<ul>
<li>
<details>
<summary><strong>Two discrete variables</strong> (e.g., <code>tau_decaymode_2</code> then <code>njets</code>).</summary>


You can also merge categories using the <code>#||#</code> syntax.

```yaml
split_categories:
  tau_decaymode_2: ["==0", "==1", "==10", "==11"]
  njets: ["==0", "==1", ">=2"]
split_categories_binedges:
  tau_decaymode_2: [-0.5, 0.5, 9.5, 11.5]
  njets: [-0.5, 0.5, 1.5, 12.5]
equipopulated_binning_options:
  variable_config:
    pt_2: 
      min: 30.0
      max: 150.0
      rounding: 4
  var_dependence_n_bins: [[9, 8, 7], [7, 7, 7], 6, 6]  # here 6 == [6, 6, 6] 

```
</details>
</li>
<li>
<details>
<summary><strong>Discrete then continuous variable</strong> (e.g., <code>njets</code> then <code>pt_1</code>).</summary>


pt_1 is binned equipopulated within each njets category.

```yaml
split_categories:
  njets: ["==0", "==1", ">=2"]
  pt_1: [] # Placeholder for equipopulated split
split_categories_binedges:
  njets: [-0.5, 0.5, 1.5, 12.5]
  pt_1: [] # Placeholder for equipopulated split
equipopulated_binning_options:
  variable_config:
    pt_2: 
      min: 30.0
      max: 150.0
      rounding: 4
    pt_1: 
      min: 0.0
      max: 3000.0
      rounding: 2
  n_bins:
    pt_1: 4   # Create 4 pt_1 bins inside each njets bin
  var_dependence_n_bins: [[9, 8, 7, 7], 10, 6]  # logic analogously to Two discrete variables example
```
</details>
</li>
<li>
<details>
<summary><strong>Two continuous variables</strong> (e.g., <code>deltaR_ditaupair</code> then <code>pt_1</code>).</summary>

```yaml
split_categories:
  deltaR_ditaupair: []
  pt_1: []
split_categories_binedges:
  deltaR_ditaupair: []
  pt_1: []
equipopulated_binning_options:
  variable_config:
    pt_2: 
        min: 30.0
        max: 150.0
        rounding: 4
    deltaR_ditaupair:
      min: 0.0
      max: 5.0
      rounding: 3
    pt_1:
      min: 0.0
      max: 3000.0
      rounding: 2
  n_bins:
    deltaR_ditaupair: 2 # Create 2 bins for the first split
    pt_1: 4             # Create 4 bins for the second split
  var_dependence_n_bins: [[10, 10, 5, 5], [8, 8, 4, 4]]
```
</details>
</li>
</ul>
