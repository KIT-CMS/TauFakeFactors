# Fake Factor calculation
In this step the fake factors are calculated. This should be run after the preselection step.

### Configuration
All information for the FF calculation step is defined in a configuration file in the `configs/ANALYSIS/ERA/` folder using the `common_settings.yaml` and a more specific config file. The `common_settings.yaml` has to be named like that and is used for all steps of the fake factor estimation (`preselection`, `FF calculation`, `FF corrections`). <br>
The FF calculation config has the following parameters:

General options for the calculation:

  parameter | type | description
  ---|---|---
  `channel` | `string` | tau pair decay channels ("et", "mt", "tt")
  `use_embedding` | `bool` | `true` if embedded sample should be used, `false` if only MC sample should be used
  `use_center_of_mass_bins` | `bool` | Changes the x-data that is entering FF and correction calculation. If set then a center of mass value is used for the x-data, calculated from events entering the corresponding bin. If not set, the bin centers are used. Default is set to True. <br> <br>This will not affect FF and correction calculation that are set to `"binwise"` (the x-data values although displayed in plots are not used)
  `stat_sigma` | `float` | This parameter defines the number of standard deviations to be considered when determining the uncertainties for fit or smoothing parameters, which are then stored in the correctionlib files. Default is `1.0`.

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
  `fit_option` | `list` | a list of polynomials that should be considered for the fake factor fits can be defined with this parameter (default: `["poly_1"]`); besides polynominal fits, it is possible to use `"smoothed"` (applies a gaussian density kernel), `"binwise"` or `"skip"`. First two can also be combined into e.g. `"binwise#[0,]+smoothed"`, `"binwise#[-1,-2,]+smoothed"` or `"binwise#[0,]#[-1,]+smoothed"` where the nth-bin(s) is computed using binwise method and the rest in smoothed manner. For left bins, positive bin index is used. For last bins, the count is done using negative integers, starting at -1 (like in the example). If both options are provided, then the binwise calculation is applied for the stated left and right bins. The `"skip"` option simply returns 1.0 for all events; this is intended for subspaces where a variable implies no events (e.g. a jet variable when `njets=0`), ensuring compatibility with other fake factors in the chain.
  `bandwidth` | `float` | if `fit_option` includes `"smoothed"` this value can be set to adjust for the bandwidth used during smoothing procedure (has no effect on the result in case of `"binwise"`). If not set the default value of histogram range divided by 5 is chosen. <br> Can either be a float value representing the bandwidth for all fake factors or a dictionary of float values, where keys correspond to the string representations of split categories defined in `split_categories`. 
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

### Running calculations
To run the FF calculation step, execute the python script and specify the config file (relative path possible):
```bash
python ff_calculation.py --config-file PATH/CONFIG.yaml
```