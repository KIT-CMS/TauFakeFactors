"""
Collection of helpful functions for other scripts
"""

import sys
import os
import glob
import array
import numpy as np
import logging
from typing import List, Dict, Union, Any

import ROOT
from XRootD import client


def check_path(path: str) -> None:
    '''
    This function checks if a given path exist. If not, this path is created.

    Args:
        path: path to check, as a string
    
    Return:
        None
    '''
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def setup_logger(log_file: str, log_name: str, subcategories: Union[List[str], None] = None) -> None:
    '''
    Setting up all relevant loggers and handlers.

    Args:
        log_file: Name of the file the logging information will be stored in
        log_name: General name of the logger
        subcategories: List of different sub logger names e.g. can be used to differentiate between processes (default: None)
    
    Return:
        None
    '''
    # create file handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(logging.INFO)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter and add it to the handlers
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    if subcategories is not None:
        for cat in subcategories:
            log = logging.getLogger(f"{log_name}.{cat}")
            log.setLevel(logging.DEBUG)
            # add the handlers to logger
            log.addHandler(ch)
            log.addHandler(fh)
    else:
        log = logging.getLogger(f"{log_name}")
        log.setLevel(logging.DEBUG)
        # add the handlers to logger
        log.addHandler(ch)
        log.addHandler(fh)


def get_ntuples(config: Dict, process: str, sample: str) -> List[str]:
    '''
    This function generates a list of paths of all ntuples for a specific sample of a process.

    Args:
        config: Dictionary with the configuration information for the preselection
        process: General name of a process e.g. "ttbar"
        sample: Exact name of the folder where the ntuples are stored e.g. "SingleMuon_Run2018A-UL2018"
    
    Return:
        List of file paths
    '''
    log = logging.getLogger(f"preselection.{process}")
    sample_path = os.path.join(
        config["ntuple_path"], config["era"], sample, config["channel"]
    )
    log.info(
        f"The following files are loaded for era: {config['era']}, channel: {config['channel']}, sample {sample}"
    )
    # now check, if the files exist
    selected_files = check_inputfiles(path=sample_path, process=process, tree=config["tree"])

    log.info("-" * 50)
    for file in selected_files: 
        log.info(file)
    log.info("-" * 50)

    return selected_files


def check_inputfiles(path: str, process: str, tree: str) -> List[str]:
    '''
    Additional function to check if the input files are empty. If yes, they are skipped.

    Args:
        path: Path where files should be checked
        process: General name of a process e.g. "ttbar"
        tree: Name of the tree in the input files

    Return:
        List of file paths with not empty files
    '''
    log = logging.getLogger(f"preselection.{process}")

    fsname = "root://cmsxrootd-kit.gridka.de/"
    xrdclient = client.FileSystem(fsname)
    status, listing = xrdclient.dirlist(path.replace(fsname, ""))

    if not status.ok:
        log.info(f"Error: {status.message}")
        sys.exit(1)

    selected_files = []
    for f in listing:
        if f.name.endswith(".root"):
            # check if file is empty
            if check_for_empty_tree(file_path=os.path.join(path, f.name), tree=tree):
                log.info(f"File {f.name} is empty. Skipping.")
                continue
            selected_files.append(os.path.join(path, f.name))

    return selected_files


def check_for_empty_tree(file_path: str, tree: str) -> bool:
    '''
    Checking if the tree in a file is empty.

    Args:
        file_path: Path of the input file
        tree: Name of the tree in the input file
    
    Return:
        Boolean which is true if the file is empty, false otherwise.
    '''
    f = ROOT.TFile.Open(file_path)
    t = f.Get(tree)

    return bool(t.GetEntries() == 0)


def rdf_is_empty(rdf: ROOT.RDataFrame) -> bool:
    '''
    Function to check if a root DataFrame is empty.

    Args:
        rdf: root DataFrame object
    
    Return:
        Boolean which is true if the root DataFrame is empty, false otherwise.
    '''
    try:
        cols = rdf.GetColumnNames()
        if len(cols) == 0:
            return True
    except:
        return True
    
    return False


def get_output_name(path: str, process: str, tau_gen_mode: str, idx: Union[int, None] = None) -> str:
    '''
    Function to generate a file output path where the file will be stored.

    Args:
        path: Path to the folder where the file should be stored at 
        process: Name of the process the file correspond to
        tau_gen_mode: Specifying the applied tau pair origin selection 
        idx: index counter, needed if a process has more than one data sample

    Return:
        String with the file path
    '''
    if tau_gen_mode == "all":
        tau_gen_mode = ""
    else:
        tau_gen_mode = "_" + tau_gen_mode

    if idx is not None:
        return os.path.join(path, f"{process}{tau_gen_mode}_{idx}.root")
    else:
        return os.path.join(path, f"{process}{tau_gen_mode}.root")


def rename_boosted_variables(rdf: Any, channel: str) -> Any:
    '''
    Function to redefine variables to the boosted tau pair information. Redefining only variables 
    which are written out for the fake factor measurement. Due to the hardcoded naming and redifinitions 
    this function needs to be adjusted if something changes in the list of output variables.

    Args:
        rdf: root DataFrame
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"

    Return:
        root DataFrame with redefined variables
    '''
    rdf = rdf.Redefine("njets", "njets_boosted")
    rdf = rdf.Redefine("nbtag", "nbtag_boosted")
    rdf = rdf.Redefine("metphi", "metphi_boosted")
    rdf = rdf.Redefine("pt_1", "boosted_pt_1")
    rdf = rdf.Redefine("q_1", "boosted_q_1")
    rdf = rdf.Redefine("pt_2", "boosted_pt_2")
    rdf = rdf.Redefine("q_2", "boosted_q_2")
    rdf = rdf.Redefine("mt_1", "boosted_mt_1")
    rdf = rdf.Redefine("deltaR_ditaupair", "boosted_deltaR_ditaupair")
    rdf = rdf.Redefine("m_vis", "boosted_m_vis")

    if "boosted_gen_match_2" in rdf.GetColumnNames():
        rdf = rdf.Redefine("gen_match_2", "boosted_gen_match_2")
    else:
        rdf = rdf.Define("boosted_gen_match_2", "-1.")
        rdf = rdf.Redefine("gen_match_2", "boosted_gen_match_2")

    if "btag_weight_boosted" in rdf.GetColumnNames():
        rdf = rdf.Redefine("btag_weight", "btag_weight_boosted")
    else:
        rdf = rdf.Define("btag_weight_boosted", "1.")
        rdf = rdf.Redefine("btag_weight", "btag_weight_boosted")

    if channel == "tt":
        if "boosted_gen_match_1" in rdf.GetColumnNames():
            rdf = rdf.Redefine("gen_match_1", "boosted_gen_match_1")
        else:
            rdf = rdf.Define("boosted_gen_match_1", "-1.")
            rdf = rdf.Redefine("gen_match_1", "boosted_gen_match_1")
    
    if channel == "et":
        rdf = rdf.Redefine("extraelec_veto", "boosted_extraelec_veto")
    if channel == "mt":
        rdf = rdf.Redefine("extramuon_veto", "boosted_extramuon_veto")

    return rdf


def get_samples(config: Dict[str, Union[str, Dict, List]]) -> List[str]:
    '''
    Function to get a list of all sample paths which will be used for the fake factor calculation. 
    This function assumes that the preselection step was already finished and takes into account 
    if embedded events or MC events should be used.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation

    Return:
        List of all paths to the relevant samples
    '''
    log = logging.getLogger("ff_calc")

    general_sample_path = os.path.join(
        config["file_path"], "preselection", config["era"], config["channel"], "*.root"
    )
    log.info(
        f"The following files are loaded for era: {config['era']}, channel: {config['channel']} from {general_sample_path}"
    )
    log.info("-" * 50)
    sample_paths = glob.glob(general_sample_path)
    tmp_list = glob.glob(general_sample_path)

    for f in tmp_list:
        sample = f.rsplit("/")[-1].rsplit(".")[0]
        if config["use_embedding"] and "_T" in sample:
            sample_paths.remove(f)
        elif not config["use_embedding"] and sample == "embedding":
            sample_paths.remove(f)

    for f in sample_paths:
        log.info(f)
    log.info("-" * 50)

    return sample_paths


def check_categories(config: Dict[str, Union[str, Dict, List]]) -> None:
    '''
    A function to check if the categories in the configuration file are defined properly. 
    Categories are defined by to parameters, the first one are the orthogonal cuts which split the data 
    into the categories and the second one is a binning corresponding to this categories. 
    The binning is needed to write a correct correctionlib file. The number of categories and bins have to match.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation

    Return:
        None
    '''
    if "target_processes" in config:
        for process in config["target_processes"]:
            categories = config["target_processes"][process]["split_categories"]
            category_edges = config["target_processes"][process]["split_categories_binedges"]
            for cat in categories:
                if len(categories[cat]) != (len(category_edges[cat]) - 1):
                    raise Exception(
                        f"Categories and binning for the categories does not match up for {process}, {cat}."
                    )
                
    if "process_fractions" in config:
        fraction_categories = config["process_fractions"]["split_categories"]
        fraction_categories_edges = config["process_fractions"]["split_categories_binedges"]
        for cat in fraction_categories:
            if len(fraction_categories[cat]) != (len(fraction_categories_edges[cat]) - 1):
                raise Exception(
                    "Categories and binning for the categories does not match up for {cat} for fractions."
                )











def modify_config(config, process, corr_config, for_AR_SR=False):
    if "SRlike_cuts" in corr_config:
        for mod in corr_config["SRlike_cuts"]:
            config["target_processes"][process]["SRlike_cuts"][mod] = corr_config[
                "SRlike_cuts"
            ][mod]
    if "ARlike_cuts" in corr_config:
        for mod in corr_config["ARlike_cuts"]:
            config["target_processes"][process]["ARlike_cuts"][mod] = corr_config[
                "ARlike_cuts"
            ][mod]
    if "AR_SR_cuts" in corr_config and for_AR_SR:
        for mod in corr_config["AR_SR_cuts"]:
            config["target_processes"][process]["ARlike_cuts"][mod] = corr_config[
                "AR_SR_cuts"
            ][mod]
            config["target_processes"][process]["SRlike_cuts"][mod] = corr_config[
                "AR_SR_cuts"
            ][mod]



def smooth_function(ff_hist, bin_edges):
    hist_bins = ff_hist.GetNbinsX()

    x = list()
    y = list()
    error_y_up = list()
    error_y_down = list()
    for nbin in range(hist_bins):
        x.append(ff_hist.GetBinCenter(nbin + 1))
        y.append(ff_hist.GetBinContent(nbin + 1))
        error_y_up.append(ff_hist.GetBinErrorUp(nbin + 1))
        error_y_down.append(ff_hist.GetBinErrorLow(nbin + 1))

    x = array.array("d", x)
    y = array.array("d", y)
    error_y_up = array.array("d", error_y_up)
    error_y_down = array.array("d", error_y_down)

    # sampling values for y based on a normal distribution with the measured statistical uncertainty
    n_samples = 20
    sampled_y = list()
    for idx in range(len(x)):
        sampled_y.append(np.random.normal(y[idx], error_y_up[idx], n_samples))
    sampled_y = np.array(sampled_y)
    sampled_y[sampled_y < 0.0] = 0.0

    # calculate widths
    fit_y_binned = list()

    n_bins = 100 * hist_bins
    for i in range(n_bins):
        fit_y_binned.append(list())

    eval_bin_edges, bin_step = np.linspace(
        bin_edges[0], bin_edges[-1], (n_bins + 1), retstep=True
    )
    bin_range = bin_edges[-1] - bin_edges[0]
    bin_half = bin_step / 2.0
    smooth_x = (eval_bin_edges + bin_half)[:-1]
    smooth_x = array.array("d", smooth_x)

    for sample in range(n_samples):
        y_arr = array.array("d", sampled_y[:, sample])
        graph = ROOT.TGraphAsymmErrors(len(x), x, y_arr, 0, 0, error_y_down, error_y_up)
        gs = ROOT.TGraphSmooth("normal")
        grout = gs.SmoothKern(graph, "normal", (bin_range / 5.0), n_bins, smooth_x)
        for i in range(n_bins):
            fit_y_binned[i].append(grout.GetPointY(i))

    smooth_y = list()
    smooth_y_up = list()
    smooth_y_down = list()

    for b in fit_y_binned:
        smooth_y.append(np.mean(b))
        smooth_y_up.append(np.std(b))
        smooth_y_down.append(np.std(b))

    smooth_y = array.array("d", smooth_y)
    smooth_y_up = array.array("d", smooth_y_up)
    smooth_y_down = array.array("d", smooth_y_down)

    corr_dict = {
        "edges": np.array(eval_bin_edges),
        "y": np.array(smooth_y),
        "up": np.array(smooth_y) + np.array(smooth_y_up),
        "down": np.array(smooth_y) - np.array(smooth_y_down),
    }

    smooth_graph = ROOT.TGraphAsymmErrors(
        len(smooth_x), smooth_x, smooth_y, 0, 0, smooth_y_down, smooth_y_up
    )
    return smooth_graph, corr_dict


def calculate_non_closure_correction(SRlike, ARlike):
    corr = SRlike["data_subtracted"].Clone()

    frac = ARlike["data_subtracted"].GetMaximum() / ARlike["data"].GetMaximum()
    predicted = ARlike["data_ff"].Clone()
    predicted.Scale(frac)

    corr.Divide(predicted)
    hist_bins = corr.GetNbinsX()
    for x in range(hist_bins):
        bincontent = corr.GetBinContent(x)
        if bincontent <= 0.0:
            corr.SetBinContent(x, 1.0)
            corr.SetBinError(x, 0.6)
    return corr, frac


def calculate_non_closure_correction_Wjets(SRlike, ARlike):
    corr = SRlike["Wjets"].Clone()
    predicted = ARlike["Wjets_ff"].Clone()
    corr.Divide(predicted)
    return corr


def calculate_non_closure_correction_ttbar(SRlike, ARlike):
    corr = SRlike["ttbar_J"].Clone()
    predicted = ARlike["ttbar_ff"].Clone()
    corr.Divide(predicted)
    return corr
