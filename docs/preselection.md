# Event preselection
This framework is designed for n-tuples (and friend trees) produced with CROWN as input. 
All information for the preselection step should be defined in configuration files in the `configs/ANALYSIS/ERA/` folder using the `common_settings.yaml` file and a more specific config file. The `common_settings.yaml` has to be named like that and is used for all steps of the fake factor estimation (`preselection`, `FF calculation`, `FF corrections`).

The preselection config has the following parameters:

  parameter | type | description
  ---|---|---
  `channel` | `string` | tau pair decay channels ("et", "mt", "tt")
  `processes` | `dict` | process parameters are explained below
  `column_definitions` | `dict` | in this section, new columns can be defined based on a given `ROOT` expression. <br>The keys of the dictionary correspond to the name of the defined column. The values are dictionaries itself, with the `expression` key defining the `ROOT` expression for defining the column. Optional entries `processes` and `exclude_processes` allow to target specific processes, the entry `allow_redefine` can be used to enable the use of the `ROOT.RDataFrame.Redefine` function for overwriting already existing columns. For a more detailed description, see below.
  `event_selection` | `dict` | with this parameter all selections that should be applied are defined. <br>This is basically a dictionary of cuts where the key is the name of a cut and the value is the cut itself as a string e.g. `had_tau_pt: "pt_2 > 30"`. The name of a cut is not really important, it is only used as an output information in the terminal. A cut can only use variables which are in the ntuples.
  `friend_weight` | `string` | optional name of the CROWN friend column used as the event weight, default `weight` (all MC weights come from CROWN's selection friend)
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

In `column_definitions`, new columns can be added to the output `ntuples` by
using `ROOT` expression. An example entry could look like this:

```yaml
column_definitions:
    nbtag:
        expression: n_bjets
        processes:
        - ttbar
        - DY
    btag_weight:
        expression: id_wgt_bjet_pnet_shape
        exclude_processes:
        - data
        allow_redefine: True
    jj_deltaR:
        expression: ROOT::VecOps::DeltaR(jeta_1, jeta_2, jphi_1, jphi_2)
```

The key `expression` is required and can contain any valid `ROOT` expression.

The entry `exclude_processes` is optional. Column definitions are performed for
all processes except the ones given in this list. The entry `processes` is also
optional. The column definition is performed only for processes in this list.
The lists `processes` and `exclude_processes` can contain the names from the
`processes` section of this configuration. By default, the new columns are
defined for all processes. To write the new columns to the output file, you have
to explicitly add the columns to the `output_features` list. Note that you can
only set `processes` or `exclude_processes` for a column, but not both at the
same time.

If the key `allow_redefine` is set to `True`, the `ROOT.RDataFrame.Redefine`
function is used if a column with the same name has been found in the
`RDataFrame`. The values in this column are then overwritten by the expression
given for the new column.

The event weight is CROWN's friend column `weight` (or the one named in `friend_weight`); no weight is built in the preselection.

To run the preselection step, execute the python script and specify the config file (relative path possible):
```bash
python preselection.py --config-file configs/PATH/CONFIG.yaml
```
Further there are additional optional parameters: 

1. `--nthreads=SOME_INTEGER` to define the number of threads for the multiprocessing pool to run the sample processing in parallel. Default value is 8 (this should normally cover running all of the samples in parallel).
2.  `--ncores=SOME_INTEGER` to define the number of cores that should be used for each pool thread to speed up the ROOT dataframe calculation. Default value is 2.
