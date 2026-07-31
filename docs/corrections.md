# Fake Factor corrections
In this step the corrections for the fake factors are calculated. This should be run after the FF calculation step.

### Configuration
Currently two different correction types are implemented: 

1. non closure correction depending on a specific variable
2. DR to SR interpolation correction depending on a specific variable

All information for the FF correction calculation step is defined in a configuration file in the `configs/ANALYSIS/ERA/` folder using the `common_settings.yaml` and a more specific config file. The `common_settings.yaml` has to be named like that and is used for all steps of the fake factor estimation (`preselection`, `FF calculation`, `FF corrections`). Additional information is loaded from the used config in the previous FF calculation step (this is done automatically). <br>
The FF correction config has the following parameters:

The expected input folder structure is workdir/WORKDIR_NAME/ERA/fake_factors/CHANNEL/*

  parameter | type | description
  ---|---|---
  `channel` | `string` | tau pair decay channels ("et", "mt", "tt")
  `correction_tag` | `string` | this parameter can be set to save a correction calculation run to a separate folder, this allow to have multiple correction calulations for a single fake factor calculation. The default folder is `default`.

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
  `split_categories` | `dict` | Optional, analogous to `Fake Factor Calculation` (only 1D).
  `var_bins` | `list` or `dict[list]` | Analogous to `var_bins` described in `Fake Factor Calculation`
  `correction_option` | `str` | Definition of how a correction should be measured. Options are `"smoothed"` (applies a gaussian density kernel), `"binwise"`, `"skip"`. First two can also be combined into e.g. `"binwise#[0,]+smoothed"`, `"binwise#[-1,-2,]+smoothed"` or `"binwise#[0,]#[-1,]+smoothed"` where the nth-bin(s) is computed using binwise method and the rest in smoothed manner. For left bins, positive bin index is used. For last bins, the count is done using negative integers, starting at -1 (like in the example). If both options are provided, then the binwise calculation is applied for the stated left and right bins. The `"skip"` option simply returns 1.0 for all events; this is intended for subspaces where a variable implies no events (e.g. a jet variable when `njets=0`), ensuring compatibility with other corrections in the chain.
  `bandwidth` | `float` | if `correction_option` includes `"smoothed"` this value can be set to adjust for the bandwidth used during smoothing procedure (has no effect on the result in case of `"binwise"`). If not set the default value of histogram range divided by 5 is chosen. <br> Can either be a float value representing the bandwidth for all corrections or a dictionary of float values, where keys correspond to the string representations of split categories defined in `split_categories`. 
  `SRlike_cuts` | `dict` | event selections for the signal-like region of the target process that should be replaced compared to the selection used in the previous FF calculation step
  `ARlike_cuts` | `dict` | event selections for the application-like region of the target process that should be replaced compared to the selection used in the previous FF calculation step
  `AR_SR_cuts` | `dict` | event selections for a switch from the determination region to the signal/application region, this is only relevant for `DR_SR` corrections

For the DR to SR correction some special parameters can be set. These have to be set as parameters of the `DR_SR` parameter:

  parameter | type | description
  ---|---|---
  `non_closure` | `dict` | since for this corrections additional fake factors are calculated. It's possible to calculated and apply non closure corrections to these fake factors before calculating the actual DR to SR correction. This parameter can be set up as discribe for non closure corrections above.
  `use_embedding` | `bool` | by setting this parameter for the `DR_SR` correction it overwrites the same parameter that is set in `common_settings.yaml`. This can be used to switch between embedding and MC only for this one `DR_SR` correction. e.g. this was historically done for QCD because in the anti isolated region the modeling of embedding is bad.
  `use_orthogonal_fake_factors` | `bool` | if set to `true` new fake factors are calculated in an separate region and also used to calculate the `DR_SR` correction. If `false` the already calculated fake factors in the previous calculation step are used.
  `compute_orthogonal_fake_factors_using_data` | `bool` | this parameter defines if the `DR_SR` correction should be calculated on data or on MC only. e.g. this is historically used for W+jets because there is no good additional control region, therefore, the region is extended to the signal region and the `DR_SR` correction is calculated only using W+jets MC.

### Statistical check
As an option it is possible to check using a sliding-window compatibility test whether a computed correction is statistically consistent with 1.0 (i.e., a flat correction, meaning no genuine fake factor bias). The test works as follows:

1. **Error selection**: Per-bin uncertainties are taken either as the fully propagated errors (data + MC statistical, added in quadrature) or as the data-only "MC-suppressed" errors — depending on the parameter `use_suppressed_mc_errors_for_correction_selection`.
2. **Sliding-window scan**: Pulls $(y_i - 1)/\sigma_i$ are computed for each bin. The test scans all contiguous windows of 1 to $N$ bins and finds the most significant deviation (minimum p-value):

$$z_\text{window} = \frac{\left|\sum_{i} \text{pull}_i\right|}{\sqrt{N_\text{window}}}$$

3. **Look-elsewhere correction**: A Bonferroni-style correction is applied to account for the scan over multiple windows:

$$p_\text{shape} = 1 - (1 - p_\text{min})^{N_\text{bins}}$$

4. **Auto-skipping**: If `p_shape > skip_corrections_p_value` (correction is compatible with 1) **and** `skip_corrections_compatible_to_one` is `true`, the correction is replaced with a flat 1.0. Bandwidth/shape variations are also set to 1.0. Statistical variations are collapsed to a single inclusive uncertainty. MC-shift variations are either fit to a constant or also set to 1.0.

Three parameter can be set to define this check:

  parameter | type | description
  ---|---|---
  `skip_corrections_compatible_to_one` | `bool` | Master switch for automatic skipping of corrections compatible with 1.0. If `True` and a correction passes the p-value threshold (see below), it is replaced by a flat 1.0 correction and all shape/bandwidth uncertainties are also flattened. Statistical uncertainties are collapsed into a single inclusive value. Defaults to `false`.
  `skip_corrections_p_value` | `float` | P-value threshold for the compatibility-with-1 test. A correction is considered compatible with 1.0 (and thus auto-skipped) if its shape p-value is above this threshold. Higher values are more aggressive (skip more corrections). Only has an effect when `skip_corrections_compatible_to_one` is set to `true`. Defaults to `0.05`.
  `use_suppressed_mc_errors_for_correction_selection` | `bool` | Controls which per-bin uncertainties are used in the compatibility test. If `True`, only the data statistical uncertainties are used — MC statistical errors from the MC subtraction step are ignored ("suppressed"). If `False`, the fully propagated errors (data + MC statistical uncertainties added in quadrature) are used. Using the suppressed errors avoids MC-dominated samples from artificially increasing the uncertainty and causing genuine corrections to be auto-skipped. Defaults to `true`.

### Running calculations
To run the FF correction step, execute the python script and specify the config file (relative path possible):
```bash
python ff_corrections.py --config-file PATH/CONFIG.yaml 
```
There are two optional parameters `--skip-DRtoSR-ffs` and `--only-main-corrections`. The correction caclulation is done in 3 steps. <br>
The first step is to calculate additional fake factors which are needed for the final DR to SR correction. If this is already done, this step can be skipped using `--skip-DRtoSR-ffs`. <br>
The second step is to calculate non closure corrections for these additional DR to SR fake factors. If both steps are already done they can be skipped by using `--only-main-corrections`. <br>
The last step is to calculate all the specified corrections for the main fake factors.
