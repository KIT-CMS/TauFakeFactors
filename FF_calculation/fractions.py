"""
Function for calculating the process fractions for the fake factors
"""

import array
import copy
import logging
from io import StringIO
from typing import Any, Dict, List, Union

import ROOT
from wurlitzer import STDOUT, pipes

import helper.ff_functions as func
import helper.plotting as plotting


def fraction_calculation(
    config: Dict[str, Union[str, Dict, List]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    logger: str,
    **kwargs: Dict[str, Any],
) -> Dict[str, Dict[str, Dict[str, List[float]]]]:
    """
    This function calculates fractions of processes for the application of fake factors.
    In general fake factors are calculated for more than one process and therefore an estimation
    for the contribution of each process is needed.

    Args:
        config: A dictionary with all the relevant information for the fraction calculation
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        process: This is relevant for fractions because for the tt channel two different fractions are calculated, one for each hadronic tau
        logger: Name of the logger that should be used

    Return:
        Dictionary where the categories are defined as keys and the values are fractions of a process
    """
    log = logging.getLogger(logger)

    # init histogram dict for the fraction calculation
    AR_hists = dict()
    SR_hists = dict()
    fractions = dict()

    # get config information for the fraction calculation
    process_conf = config[process]

    split_variables, split_combinations = func.get_split_combinations(
        categories=process_conf["split_categories"]
    )

    # splitting between different categories
    for split in split_combinations:
        for sample_path in sample_paths:
            # getting the name of the process from the sample path
            sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
            log.info(
                f"Processing {sample} for the {', '.join([f'{var} {split[var]}' for var in split_variables])} category."
            )
            log.info("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)

            # event filter for application region
            region_conf = copy.deepcopy(process_conf["AR_cuts"])
            rdf_AR = func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=split,
                region_cuts=region_conf,
            )

            log.info(
                "Filtering events for the fraction calculation in the application region."
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_AR.Report().Print()
            log.info(out.getvalue())
            log.info("-" * 50)

            # event filter for signal region; this is not needed for the FF calculation, just for control plots
            region_conf = copy.deepcopy(process_conf["SR_cuts"])
            rdf_SR = func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=split,
                region_cuts=region_conf,
            )

            log.info(
                "Filtering events for the fraction calculation in the signal region."
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_SR.Report().Print()
            log.info(out.getvalue())
            log.info("-" * 50)

            # get binning of the dependent variable
            xbinning = array.array("d", process_conf["var_bins"])
            nbinsx = len(process_conf["var_bins"]) - 1

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
        AR_hists["QCD"] = func.QCD_SS_estimate(hists=AR_hists)
        SR_hists["QCD"] = func.QCD_SS_estimate(hists=SR_hists)

        frac_hists = dict()

        for p in config[process]["processes"]:
            frac_hists[p] = func.calc_fraction(
                hists=AR_hists,
                target=p,
                processes=config[process]["processes"],
            )
        frac_hists = func.add_fraction_variations(
            hists=frac_hists,
            processes=config[process]["processes"],
        )

        SR_frac_hists = dict()

        for p in config[process]["processes"]:
            SR_frac_hists[p] = func.calc_fraction(
                hists=SR_hists,
                target=p,
                processes=config[process]["processes"],
            )
        SR_frac_hists = func.add_fraction_variations(
            hists=SR_frac_hists,
            processes=config[process]["processes"],
        )

        try:
            cat_string = f"{split_variables[0]}#{split[split_variables[0]]}"
        except:
            raise Exception(
                "Category splitting for fractions is only defined up to 1 dimensions."
            )

        fractions[cat_string] = frac_hists

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
        )

        data = "data"
        if config["use_embedding"]:
            samples = [
                "QCD",
                "diboson_J",
                "diboson_L",
                "Wjets",
                "ttbar_J",
                "ttbar_L",
                "DYjets_J",
                "DYjets_L",
                "ST_J",
                "ST_L",
                "embedding",
            ]
        else:
            samples = [
                "QCD",
                "diboson_J",
                "diboson_L",
                "diboson_T",
                "Wjets",
                "ttbar_J",
                "ttbar_L",
                "ttbar_T",
                "DYjets_J",
                "DYjets_L",
                "DYjets_T",
                "ST_J",
                "ST_L",
                "ST_T",
            ]
        
        for _hist, _region, _data, _samples in [
            (SR_hists, "SR", data, samples),
            (AR_hists, "AR", data, samples),
        ]:
            plotting.plot_data_mc_ratio(
                variable=process_conf["var_dependence"],
                hists=_hist,
                era=config["era"],
                channel=config["channel"],
                process=process,
                region=_region,
                data=_data,
                samples=_samples,
                category=split,
                output_path=output_path,
                logger=logger,
            )
        log.info("-" * 50)

    # transform histograms to fraction values for correctionlib
    fractions = func.get_yields_from_hists(
        hists=fractions,
        processes=config[process]["processes"],
    )

    return fractions
