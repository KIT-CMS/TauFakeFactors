"""
Function for calculating the process fractions for the fake factors
"""

import array
import copy
import logging
from io import StringIO
from typing import Any, Dict, List, Union, Tuple

import ROOT
from wurlitzer import STDOUT, pipes

import helper.ff_functions as ff_func
import helper.plotting as plotting


def fraction_calculation(
    args: Tuple[Any, ...],
) -> Dict[str, Dict[str, Dict[str, List[float]]]]:
    """
    This function calculates the fractions of the processes for the fake factor calculation for a given category.

    Args:
        args: Tuple containing all the necessary information for the calculation of the fake factors
            split: Dictionary containing the category information
            binning: List of bin edges for the dependent variable
            config: Dictionary with all the relevant information for the fake factor calculation
            process_conf: Dictionary with all the relevant information for the fake factor calculation of the specific process
            process: Name of the process
            split_variables: List of variables that are used for the category splitting
            sample_paths: List of file paths where the samples are stored
            output_path: Path where the generated plots should be stored
            logger: Name of the logger that should be used

    Return:
        Dictionary with the category information as key and the fractions of the processes as value
    """
    (
        split,  # split: Dict[str, str],
        binning,  # binning: List[float],
        config,  # config: Dict[str, Union[str, Dict, List]],
        process_conf,  # process_conf: Dict[str, Union[str, Dict, List]],
        process,  # process: str,
        split_variables,  # split_variables: List[str],
        sample_paths,  # sample_paths: List[str],
        output_path,  # output_path: str,
        logger,  # logger: str,
        *_,  # SRlike_hists, ARlike_hists only needed for ttbar
    ) = args

    log = logging.getLogger(logger)

    AR_hists = dict()
    SR_hists = dict()

    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        log.info(f"Processing {sample} for the {', '.join([f'{var} {split[var]}' for var in split_variables])} category.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for application region
        region_conf = copy.deepcopy(process_conf["AR_cuts"])
        rdf_AR = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=split,
            region_cuts=region_conf,
        )

        log.info("Filtering events for the fraction calculation in the application region.")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_AR.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # event filter for signal region; this is not needed for the FF calculation, just for control plots
        region_conf = copy.deepcopy(process_conf["SR_cuts"])
        rdf_SR = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=split,
            region_cuts=region_conf,
        )

        log.info("Filtering events for the fraction calculation in the signal region.")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SR.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # get binning of the dependent variable
        xbinning = array.array("d", binning)
        nbinsx = len(binning) - 1

        # making the histograms
        h = rdf_AR.Histo1D(
            (process_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            process_conf["var_dependence"],
            "weight",
        )
        AR_hists[sample] = h.GetValue()

        h = rdf_SR.Histo1D(
            (process_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            process_conf["var_dependence"],
            "weight",
        )
        SR_hists[sample] = h.GetValue()

    # calculate QCD estimation; here directly estimated as difference between mc and data without SS/OS
    AR_hists["QCD"] = ff_func.QCD_SS_estimate(hists=AR_hists)
    SR_hists["QCD"] = ff_func.QCD_SS_estimate(hists=SR_hists)

    frac_hists = dict()

    for p in config[process]["processes"]:
        frac_hists[p] = ff_func.calc_fraction(
            hists=AR_hists,
            target=p,
            processes=config[process]["processes"],
        )
    frac_hists = ff_func.add_fraction_variations(
        hists=frac_hists,
        processes=config[process]["processes"],
    )

    SR_frac_hists = dict()

    for p in config[process]["processes"]:
        SR_frac_hists[p] = ff_func.calc_fraction(
            hists=SR_hists,
            target=p,
            processes=config[process]["processes"],
        )
    SR_frac_hists = ff_func.add_fraction_variations(
        hists=SR_frac_hists,
        processes=config[process]["processes"],
    )

    plotting.plot_fractions(
        variable=process_conf["var_dependence"],
        hists=frac_hists["nominal"],
        era=config["era"],
        channel=config["channel"],
        region="AR",
        fraction_name=process,
        processes=config[process]["processes"],
        category=split,
        output_path=output_path,
        logger=logger,
        save_data=True,
    )
    plotting.plot_fractions(
        variable=process_conf["var_dependence"],
        hists=SR_frac_hists["nominal"],
        era=config["era"],
        channel=config["channel"],
        region="SR",
        fraction_name=process,
        processes=config[process]["processes"],
        category=split,
        output_path=output_path,
        logger=logger,
        save_data=True,
    )

    # producing some control plots
    for _hist, _region in [
        (SR_hists, "SR"),
        (AR_hists, "AR"),
    ]:
        for yscale, save_data in zip(["linear", "log"], [True, False]):
            plotting.plot_data_mc_ratio(
                variable=process_conf["var_dependence"],
                hists=_hist,
                era=config["era"],
                channel=config["channel"],
                process=process,
                region=_region,
                data="data",
                samples=ff_func.controlplot_samples(config["use_embedding"]),
                category=split,
                output_path=output_path,
                logger=logger,
                yscale=yscale,
                save_data=save_data,
            )
    log.info("-" * 50)

    return {f"{split_variables[0]}#{split[split_variables[0]]}": frac_hists}
