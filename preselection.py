"""
Script for preprocessing n-tuples for the fake factor calculation
"""

import argparse
import logging
import os
from io import StringIO
from typing import Dict, List, Tuple, Union

import ROOT
from wurlitzer import STDOUT, pipes

import helper.filters as filters
import helper.functions as func
import CustomLogging as logging_helper

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config-file",
    default=None,
    help="Path to the config file which contains information for the preselection step.",
)
parser.add_argument(
    "--nthreads",
    default=8,
    help="Number of threads to use for the multiprocessing pool in the preselection step. (default: 8)",
)
parser.add_argument(
    "--ncores",
    default=2,
    help="Number of cores to use for each pool the preselection step. (default: 2)",
)
parser.add_argument(
    "--disable-multiprocessing",
    action="store_true",
    help="Flag to disable multiprocessing for debugging purposes.",
)
parser.add_argument(
    "--log-level",
    default="INFO",
    help="Logging level to use. (default: INFO)",
)


def resolve_output_features(
    rdf: "ROOT.RDataFrame",
    output_features: List[str],
    log: logging.Logger,
) -> List[str]:
    """
    Checks the requested output features against the columns actually available on a data frame.

    Missing physics columns are a hard error, as before. Missing CROWN selection mask columns
    (presel_mask, sel_os/sel_ss, ff_*) are only warned about and dropped, so that configs which
    already request the masks can still be run on older n-tuples which do not contain them.

    Args:
        rdf: root DataFrame that is going to be snapshotted
        output_features: Requested output columns
        log: Logger to report dropped columns

    Return:
        The list of output features that are actually available
    """
    available = {str(x).replace("ntuple.", "") for x in rdf.GetColumnNames()}
    missing = [x for x in output_features if x not in available]

    droppable = [x for x in missing if func.is_selection_mask_column(x)]
    required = [x for x in missing if x not in droppable]

    if required:
        raise ValueError(f"Missing columns: {required}")
    if droppable:
        log.warning(
            f"The following CROWN selection mask columns are not available on the input "
            f"n-tuple and are dropped from the output: {droppable}"
        )

    return [x for x in output_features if x not in droppable]


@logging_helper.LogDecorator().grouped_logs(extractor=lambda args: f"preselection.{args[0]}")
def run_sample_preselection(args: Tuple[str, Dict[str, Union[Dict, List, str]], int, str, str]) -> Tuple[str, str]:
    """
    This function can be used for multiprocessing. It runs the preselection step for a specified process.

    Args:
        args: Tuple with a process name, a configuration for this process, the output directory,
            number of threads to use, a sample name and a tau gen. level mode (e.g. "T", "J", "L", "all")

    Return:
        A tuple with the tau gen. level mode and the name of the output file
    """
    process, config, output_path, ncores, sample, tau_gen_mode = args
    log = logging.getLogger(f"preselection.{process}")
    ROOT.EnableImplicitMT(ncores)

    # loading ntuple files
    ntuple_list = func.get_ntuples(config=config, process=process, sample=sample)
    if not ntuple_list:
        log.info(f"WARNING: Sample {sample} has no non-empty files. Skipping...")
        return ()
    chain = ROOT.TChain(config["tree"])

    for ntuple in ntuple_list:
        chain.Add(ntuple)

    for friend_tag in config.get("friends", []) or []:
        try:
            friend_list = func.get_friend_ntuples(
                config=config, process=process, sample=sample,
                ntuple_list=ntuple_list, tag=friend_tag,
            )
        except FileNotFoundError as e:
            log.warning(
                f"Friend tree {friend_tag!r} not found, skipping it for sample {sample!r} "
                f"(falling back to cut strings wherever it would have carried a mask): {e}"
            )
            continue
        fchain = ROOT.TChain(config["tree"])
        for friend_file in friend_list:
            fchain.Add(friend_file)
        chain.AddFriend(fchain)
        log.info(f"Attached friend tree {friend_tag!r} ({len(friend_list)} files).")

    rdf = ROOT.RDataFrame(chain)

    if func.rdf_is_empty(rdf=rdf):
        log.info(f"WARNING: Sample {sample} is empty. Skipping...")
        return ()

    # gen_match_1 does not exist on the input n-tuple for "et"/"mt" 
    if "gen_match_1" not in rdf.GetColumnNames():
        rdf = rdf.Define("gen_match_1", "-1")

    # get column definitions from config and declare definitions on the RDataFrame
    column_definitions = config.get("column_definitions", {})
    if column_definitions:
        rdf = func.define_columns(rdf, column_definitions, process)

    # apply analysis specific event filters. Configs whose input ntuples carry CROWN's
    # presel_mask branch declare a single mask
    selection_conf = config["event_selection"]
    for cut in selection_conf:
        rdf = rdf.Filter(f"({selection_conf[cut]})", f"cut on {cut}")

    if process == "embedding":
        rdf = filters.emb_tau_gen_match(rdf=rdf, channel=config["channel"])

    # CROWN's friend weight (all SFs, 1 for data/embedding); `friend_weight` picks another column, e.g. without the btag SF
    friend_weight = config.get("friend_weight", "weight")
    if friend_weight not in {str(x).replace("ntuple.", "") for x in rdf.GetColumnNames()}:
        raise ValueError(f"No CROWN friend column {friend_weight} on the input n-tuple of {sample}, attach the selection friend.")
    if friend_weight != "weight":
        rdf = rdf.Redefine("weight", friend_weight)

    if process == "embedding":
        emb_weight_conf = config["emb_weights"]
        for weight in emb_weight_conf:
            rdf = rdf.Redefine("weight", f"weight * ({emb_weight_conf[weight]})")

    # default values for some output variables which are not defined in data, embedding; will not be used in FF calculation
    if process in ["data", "embedding"]:
        if "btag_weight" not in rdf.GetColumnNames():
            rdf = rdf.Define("btag_weight", "1.")
        for idx in ["2", "1"] if config["channel"] == "tt" else ["2"]:
            for wp in config.get("tau_vs_jet_wgt_wps", []):
                if f"id_wgt_tau_vsJet_{wp}_{idx}" not in rdf.GetColumnNames():
                    rdf = rdf.Define(f"id_wgt_tau_vsJet_{wp}_{idx}", "1.")

    # for data set gen_match to -1 
    if process == "data":
        rdf = rdf.Define("gen_match_2", "-1.")

    # splitting data frame based on the tau origin (genuine, jet fake, lepton fake)
    tmp_rdf = rdf
    if tau_gen_mode != "all":
        tmp_rdf = filters.tau_origin_split(
            rdf=tmp_rdf, channel=config["channel"], tau_gen_mode=tau_gen_mode
        )

    # redirecting C++ stdout for Report() to python stdout
    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        tmp_rdf.Report().Print()
    log.debug(out.getvalue())
    log.debug("-" * 50)

    tmp_file_name = func.get_output_name(
        path=output_path, process=sample, tau_gen_mode=tau_gen_mode
    )
    # check for empty data frame -> only save/calculate if event number is not zero
    if tmp_rdf.Count().GetValue() != 0:
        log.info(f"The current data frame will be saved to {tmp_file_name}")
        features = resolve_output_features(
            rdf=tmp_rdf, output_features=output_features, log=log
        )

        tmp_rdf.Snapshot(config["tree"], tmp_file_name, features)
        log.info("-" * 50)
    else:
        log.info("No events left after filters. Data frame will not be saved.")
        log.info("-" * 50)
        return ()

    return (tau_gen_mode, tmp_file_name)


@logging_helper.LogDecorator().grouped_logs(extractor=lambda args: f"preselection.{args[0]}")
def run_preselection(args: Tuple[str, Dict[str, Union[Dict, List, str]], str, int]) -> None:
    """
    This function can be used for multiprocessing. It runs the preselection step for a specified process.

    Args:
        args: Tuple with a process name, a configuration for this process, the output directory and
            number of threads to use

    Return:
        None
    """
    process, config, output_path, ncores = args
    log = logging.getLogger(f"preselection.{process}")

    log.info(f"Processing process: {process}")
    # bookkeeping of samples files due to splitting based on the tau origin (genuine, jet fake, lepton fake)
    process_file_dict = dict()
    for tau_gen_mode in config["processes"][process]["tau_gen_modes"]:
        process_file_dict[tau_gen_mode] = list()
    process_samples = func.resolve_process_samples(config=config, process=process)
    log.info(f"Considered samples for process {process}: {process_samples}")

    # going through all contributing samples for the process
    args_list = [
        (process, config, output_path, ncores, sample, tau_gen_mode) for tau_gen_mode in config["processes"][process]["tau_gen_modes"] for sample in process_samples
    ]

    # Bounded independently from the outer (per-process) pool in __main__ below 
    _inner_cap = os.environ.get("TFF_PRESEL_INNER_WORKERS")
    results = func.optional_process_pool(
        args_list=args_list,
        function=run_sample_preselection,
        max_workers=int(_inner_cap) if _inner_cap else None,
    )
    for files in results:
        if files:
            process_file_dict[files[0]].append(files[1])

    # combining all files of a process and tau origin
    for tau_gen_mode in config["processes"][process]["tau_gen_modes"]:
        out_file_name = func.get_output_name(
            path=output_path, process=process, tau_gen_mode=tau_gen_mode
        )
        # combining sample files to a single process file, if there are any
        if len(process_file_dict[tau_gen_mode]) != 0:
            sum_rdf = ROOT.RDataFrame(config["tree"], process_file_dict[tau_gen_mode])
            log.info(
                f"The processed files for the {process} process are concatenated. The data frame will be saved to {out_file_name}"
            )
            sum_rdf.Snapshot(
                config["tree"],
                out_file_name,
                resolve_output_features(rdf=sum_rdf, output_features=output_features, log=log),
            )
            log.info("-" * 50)
        else:
            log.info(
                f"No processed files for the {process} process. An empty data frame will be saved to {out_file_name}"
            )
            # create an empty root file and save it
            f = ROOT.TFile(out_file_name, "RECREATE")
            t = ROOT.TTree(config["tree"], config["tree"])
            t.Write()
            f.Close()
            log.info("-" * 50)

        # delete not needed temporary sample files after combination
        for rf in process_file_dict[tau_gen_mode]:
            os.remove(rf)
        log.info("-" * 50)


if __name__ == "__main__":
    args = parser.parse_args()

    func.RuntimeVariables.USE_MULTIPROCESSING = not args.disable_multiprocessing

    # loading of the chosen config file
    config = func.load_config(args.config_file)

    # define output path for the preselected samples
    output_path = os.path.join(
        config["output_path"], "preselection", config["era"], config["channel"]
    )
    func.check_path(path=output_path)

    logging_helper.LOG_LEVEL = getattr(logging, args.log_level.upper(), logging.INFO)
    log = logging_helper.setup_logging(
        output_file=output_path + "/preselection.log",
        logger=logging.getLogger("preselection"),
        level=logging_helper.LOG_LEVEL,
    )

    # get needed features for fake factor calculation, from TAUER's common/variables.yaml
    # (single source of truth, see func.load_output_features) rather than a per-config copy
    output_features = func.load_output_features(config)

    # per-config additions on top of the shared variables.yaml list
    for feature in config.get("output_features", []) or []:
        if feature not in output_features:
            output_features.append(feature)

    if config["channel"] not in ["mm", "ee", "em"]:
        for wp in config["tau_vs_jet_wps"]:
            output_features.append("id_tau_vsJet_" + wp + "_2")
        for wp in config.get("tau_vs_jet_wgt_wps", []):
            output_features.append("id_wgt_tau_vsJet_" + wp + "_2")

    if config["channel"] == "tt":
        for wp in config["tau_vs_jet_wps"]:
            output_features.append("id_tau_vsJet_" + wp + "_1")
        for wp in config.get("tau_vs_jet_wgt_wps", []):
            output_features.append("id_wgt_tau_vsJet_" + wp + "_1")

    log.info(f"Used output features: {output_features}")

    # going through all wanted processes and run the preselection function with a pool of 8 workers
    args_list = [(process, config, output_path, int(args.ncores)) for process in config["processes"]]

    _outer_cap = os.environ.get("TFF_PRESEL_OUTER_WORKERS")
    func.optional_process_pool(
        args_list=args_list,
        function=run_preselection,
        max_workers=int(_outer_cap) if _outer_cap else None,
    )

    # dumping config to output directory for documentation
    with open(output_path + "/config.yaml", "w") as config_file:
        func.configured_yaml.dump(config, config_file)

    log.info("Preselection finished.")
