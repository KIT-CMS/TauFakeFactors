# Automated Equipopulated Binning

The framework includes a utility script, `adjust_binning.py`, for an automatic calculation and update binning definitions within the configuration in case of equipopulated binning approach of specified variable (`var_dependence`) across various categories for a more robust statistical treatment in each analysis region. The general intention is to have an equipopulated data distribution in the SRlike region given possible additional splittings (i.e. in `njets`).

The script is intended to be run berofe the fake factor and correction calculation by using

```bash
python adjust_binning.py --config config.yaml --cut-config config_with_region_cuts.yaml --processes QCD Wjets --cut-region ARlike
```

| Argument | Type | Description |
|---|---|---|
| `--config` | `string` | Path to the YAML configuration file that will be modified. |
| `--cuts-config` | `string` | (Optional) Path to a separate YAML file containing cut definitions. This is useful for correction YAML files where cuts are not specified directly for each correction but then are sourced from the cuts-config file, which can just be the fake factor YAML configuration. |
| `--processes` | `list[string]` | A list of processes to adjust (e.g., `QCD`, `Wjets`, `ttbar`, `process_fractions`). The default is set to include all processes. |
| `--cut-region` | `string` | The region cut to be applied when determining event counts (e.g., `SRlike`, `ARlike`). `process_fractions` always use the `AR_cut`. |
| `--dry-run` | `bool` | Preview the binning calculations and changes without modifying the config file. |

## Configuration

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


## Advanced Splitting Strategies

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


`pt_1` is binned equipopulated within each `njets` category.

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