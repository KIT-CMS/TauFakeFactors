# Fake Factor calculation
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