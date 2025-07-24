"""
Script for preprocessing n-tuples for the fake factor calculation
"""

import argparse
import json
import logging
import os
from io import StringIO
from typing import Dict, List, Tuple, Union

import ROOT
from wurlitzer import STDOUT, pipes

import helper.filters as filters
import helper.functions as func
import helper.weights as weights

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
    help="Number of cores to use for each pool the preselection step. (default: 4)",
)
parser.add_argument(
    "--disable-multiprocessing",
    action="store_true",
    help="Flag to disable multiprocessing for debugging purposes.",
)


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
    chain = ROOT.TChain(config["tree"])

    for ntuple in ntuple_list:
        chain.Add(ntuple)

    if "friends" in config:
        for friend in config["friends"]:
            friend_list = []
            for ntuple in ntuple_list:
                friend_list.append(
                    ntuple.replace("CROWNRun", "CROWNFriends/" + friend)
                )
            fchain = ROOT.TChain(config["tree"])
            for friend in friend_list:
                fchain.Add(friend)
            chain.AddFriend(fchain)

    rdf = ROOT.RDataFrame(chain)
    
    if func.rdf_is_empty(rdf=rdf):
        log.info(f"WARNING: Sample {sample} is empty. Skipping...")
        return ()

    # apply analysis specific event filters
    selection_conf = config["event_selection"]
    for cut in selection_conf:
        rdf = rdf.Filter(f"({selection_conf[cut]})", f"cut on {cut}")

    if process == "embedding":
        rdf = filters.emb_boostedtau_gen_match(rdf=rdf, channel=config["channel"])

    # calculate event weights
    rdf = rdf.Define("weight", "1.")

    mc_weight_conf = config["mc_weights"]
    if process not in ["data", "embedding"]:
        for weight in mc_weight_conf:
            if weight == "generator":
                # calculating generator weight (including cross section weight)
                if "stitching" in mc_weight_conf["generator"]:
                    if process in mc_weight_conf["generator"]["stitching"]:
                        rdf = weights.stitching_gen_weight(
                            rdf=rdf,
                            era=config["era"],
                            process=process,
                            sample_info=datasets[sample],
                        )
                    else:
                        rdf = weights.gen_weight(
                            rdf=rdf, sample_info=datasets[sample]
                        )
                else:
                    rdf = weights.gen_weight(rdf=rdf, sample_info=datasets[sample])
            elif weight == "lumi":
                rdf = weights.lumi_weight(rdf=rdf, era=config["era"])
            elif weight == "Z_pt_reweighting":
                if process in ["DYjets"]:
                    rdf = rdf.Redefine(
                        "weight", f"weight * ({mc_weight_conf[weight]})"
                    )
            elif weight == "Top_pt_reweighting":
                if process == "ttbar":
                    rdf = rdf.Redefine(
                        "weight", f"weight * ({mc_weight_conf[weight]})"
                    )
            else:
                rdf = rdf.Redefine("weight", f"weight * ({mc_weight_conf[weight]})")

    if process == "embedding":
        emb_weight_conf = config["emb_weights"]
        for weight in emb_weight_conf:
            rdf = rdf.Redefine("weight", f"weight * ({emb_weight_conf[weight]})")

    # default values for some output variables which are not defined in data, embedding; will not be used in FF calculation
    if process in ["data", "embedding"]:
        if "btag_weight" not in rdf.GetColumnNames():
            rdf = rdf.Define("btag_weight", "1.")
        for wp in config["tau_iso_wgt_wps"]:
            weightname = "id_wgt_boostedtau_iso_" + wp + "_2"
            if weightname not in rdf.GetColumnNames():
                rdf = rdf.Define(weightname, "1.")
        if config["channel"] == "tt":
            for wp in config["tau_iso_wgt_wps"]:
                weightname = "id_wgt_boostedtau_iso_" + wp + "_1"
                if weightname not in rdf.GetColumnNames():
                    rdf = rdf.Define(weightname, "1.")

    # for data set gen_match to -1
    if process == "data":
        rdf = rdf.Define("gen_match_2", "-1.")
        if config["channel"] == "tt":
            rdf = rdf.Define("gen_match_1", "-1.")

    # splitting data frame based on the tau origin (genuine, jet fake, lepton fake)
    tmp_rdf = rdf
    if tau_gen_mode != "all":
        tmp_rdf = filters.boostedtau_origin_split(
            rdf=tmp_rdf, channel=config["channel"], tau_gen_mode=tau_gen_mode
        )

    # redirecting C++ stdout for Report() to python stdout
    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        tmp_rdf.Report().Print()
    log.debug(out.getvalue())
    log.debug("-" * 50)

    # WARNING: cross check this function is something changes in the list of output features
    tmp_rdf = func.rename_boosted_variables(
        rdf=tmp_rdf, channel=config["channel"]
    )

    tmp_file_name = func.get_output_name(
        path=output_path, process=sample, tau_gen_mode=tau_gen_mode
    )
    # check for empty data frame -> only save/calculate if event number is not zero
    if tmp_rdf.Count().GetValue() != 0:
        log.info(f"The current data frame will be saved to {tmp_file_name}")
        cols = tmp_rdf.GetColumnNames()
        cols_with_friends = [str(x).replace("ntuple.", "") for x in cols]
        missing_cols = [x for x in output_features if x not in cols_with_friends]
        if len(missing_cols) != 0:
            raise ValueError(f"Missing columns: {missing_cols}")

        tmp_rdf.Snapshot(config["tree"], tmp_file_name, output_features)
        log.info("-" * 50)
    else:
        log.info("No events left after filters. Data frame will not be saved.")
        log.info("-" * 50)
        return ()

    return (tau_gen_mode, tmp_file_name)


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
    log.info(
        f"Considered samples for process {process}: {config['processes'][process]['samples']}"
    )

    # going through all contributing samples for the process
    args_list = [
        (process, config, output_path, ncores, sample, tau_gen_mode) for tau_gen_mode in config["processes"][process]["tau_gen_modes"] for sample in config["processes"][process]["samples"]
    ]

    results = func.optional_process_pool(
        args_list=args_list,
        function=run_sample_preselection,
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
            sum_rdf.Snapshot(config["tree"], out_file_name, output_features)
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

    # loading general dataset info file for xsec and event number
    with open("datasets/datasets.json", "r") as file:
        datasets = json.load(file)

    # define output path for the preselected samples
    output_path = os.path.join(
        config["output_path"], "preselection", config["era"], config["channel"]
    )
    func.check_path(path=output_path)

    func.setup_logger(
        log_file=output_path + "/preselection.log",
        log_name="preselection",
        log_level=logging.INFO,
        subcategories=config["processes"],
    )

    # get needed features for fake factor calculation
    output_features = config["output_features"]

    for wp in config["tau_iso_wps"]:
        output_features.append("id_boostedtau_iso_" + wp + "_2")
    for wp in config["tau_iso_wgt_wps"]:
        output_features.append("id_wgt_boostedtau_iso_" + wp + "_2")

    if config["channel"] == "tt":
        for wp in config["tau_iso_wps"]:
            output_features.append("id_boostedtau_iso_" + wp + "_1")
        for wp in config["tau_iso_wgt_wps"]:
            output_features.append("id_wgt_boostedtau_iso_" + wp + "_1")

    # going through all wanted processes and run the preselection function with a pool of 8 workers
    args_list = [(process, config, output_path, int(args.ncores)) for process in config["processes"]]

    func.optional_process_pool(
        args_list=args_list,
        function=run_preselection,
    )

    # dumping config to output directory for documentation
    with open(output_path + "/config.yaml", "w") as config_file:
        func.configured_yaml.dump(config, config_file)
