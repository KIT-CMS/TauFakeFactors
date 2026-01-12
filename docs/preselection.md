# Event preselection
This framework is designed for n-tuples (and friend trees) produced with CROWN as input. 
All information for the preselection step is defined in configuration files in the `configs/ANALYSIS/ERA/` folder using the `common_settings.yaml` file and a more specific config file. 

The preselection config has the following parameters:

  parameter | type | description
  ---|---|---
  `channel` | `string` | tau pair decay channels ("et", "mt", "tt")
  `processes` | `dict` | process parameters are explained below
  `column_definitions` | `dict` | in this section, new columns can be defined
  based on a given `ROOT` expression. <br>The keys of the dictionary correspond
  to the name of the defined column. The values are dictionaries itself, with
  the `expression` key defining the `ROOT` expression for defining the column and
  the optional entry `exclude_processes` containing a list of processes for
      which the column should not be added. An example is given below.
  `event_selection` | `dict` | with this parameter all selections that should be applied are defined. <br>This is basically a dictionary of cuts where the key is the name of a cut and the value is the cut itself as a string e.g. `had_tau_pt: "pt_2 > 30"`. The name of a cut is not really important, it is only used as an output information in the terminal. A cut can only use variables which are in the ntuples.
  `mc_weights` | `dict` | weight parameter are defined below
  `emb_weights` | `dict` | all weights that should be applied for embedded samples are defined. <br>Like for `event_selection` a weight can directly be specified and is then applied to all samples the same way e.g. `single_trigger: "trg_wgt_single_mu24ormu27"`
  `output_features` | `list` | the to be saved/needed features for the later calculations are listed, this also includes features from friend trees if `friends` were specified in `common_settings.yaml`

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

In `column_definitions`, new columns. An example entry could look like this:

```yaml
column_definitions:
    nbtag:
        expression: n_bjets
    btag_weight:
        expression: id_wgt_bjet_pnet_shape
        exclude_processes:
        - data
    jj_deltaR:
        expression: ROOT::VecOps::DeltaR(jeta_1, jeta_2, jphi_1, jphi_2)
```

The key `expression` is required and can contain any valid `ROOT` expression.
The entry `exclude_processes` is optional. This list can contain process names
from the `processes` section of this configuration. By default, the new columns
are defined for all processes.

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
