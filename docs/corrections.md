# Fake Factor corrections
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
  `chain_DR_SR_to_non_closure` | `bool` | Option to chain `DR_SR` correction (computed first) to the following `non_closure` corrections

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
