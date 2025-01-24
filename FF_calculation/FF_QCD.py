"""
Function for calculating fake factors for the QCD process
"""

import array
import copy
import logging
from io import StringIO
from typing import Any, Dict, List, Union

import numpy as np
import ROOT
from wurlitzer import STDOUT, pipes

import helper.ff_functions as func
import helper.plotting as plotting
from helper.ff_evaluators import FakeFactorCorrectionEvaluator, FakeFactorEvaluator


def calculation_QCD_FFs(
    config: Dict[str, Union[str, Dict, List]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    logger: str,
    **kwargs: Dict[str, Any],
) -> Dict[str, Union[Dict[str, str], Dict[str, Dict[str, str]]]]:
    """
    This function calculates fake factors for QCD.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        process: This is relevant for QCD because for the tt channel two different QCD fake factors are calculated, one for each hadronic tau
        logger: Name of the logger that should be used

    Return:
        Dictionary where the categories are defined as keys and and the values are the fitted functions (including variations)
        e.g. corrlib_expressions[CATEGORY_1][CATEGORY_2][VARIATION] if dimension of categories is 2
    """
    log = logging.getLogger(logger)

    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()
    # init dictionary for the FF functions for correctionlib
    corrlib_expressions = dict()

    # get QCD specific config information
    process_conf = config["target_processes"][process]
    split_variables, split_combinations = func.get_split_combinations(
        categories=process_conf["split_categories"]
    )

    # splitting between different categories
    for split in split_combinations:
        for sample_path in sample_paths:
            # getting the name of the process from the sample path
            sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
            log.info(
                f"Processing {sample} for the {', '.join(['{} {}'.format(var, split[var]) for var in split_variables])} category."
            )
            log.info("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)

            # event filter for QCD signal-like region
            region_conf = copy.deepcopy(process_conf["SRlike_cuts"])
            rdf_SRlike = func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=split,
                region_cuts=region_conf,
            )

            log.info(
                f"Filtering events for the signal-like region. Target process: {process}"
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_SRlike.Report().Print()
            log.info(out.getvalue())
            log.info("-" * 50)

            # event filter for QCD application-like region
            region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
            rdf_ARlike = func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=split,
                region_cuts=region_conf,
            )

            log.info(
                f"Filtering events for the application-like region. Target process: {process}"
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_ARlike.Report().Print()
            log.info(out.getvalue())
            log.info("-" * 50)

            # get binning of the dependent variable
            xbinning = array.array("d", process_conf["var_bins"])
            nbinsx = len(process_conf["var_bins"]) - 1

            # making the histograms
            h = rdf_SRlike.Histo1D(
                (process_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
                process_conf["var_dependence"],
                "weight",
            )
            SRlike_hists[sample] = h.GetValue()

            h = rdf_ARlike.Histo1D(
                (process_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
                process_conf["var_dependence"],
                "weight",
            )
            ARlike_hists[sample] = h.GetValue()

        # calculate QCD enriched data by subtraction all the background samples
        SRlike_hists["data_subtracted"] = SRlike_hists["data"].Clone()
        ARlike_hists["data_subtracted"] = ARlike_hists["data"].Clone()
        SRlike_hists["data_subtracted_up"] = SRlike_hists["data"].Clone()
        ARlike_hists["data_subtracted_up"] = ARlike_hists["data"].Clone()
        SRlike_hists["data_subtracted_down"] = SRlike_hists["data"].Clone()
        ARlike_hists["data_subtracted_down"] = ARlike_hists["data"].Clone()

        for hist in SRlike_hists:
            if hist not in [
                "data",
                "data_subtracted",
                "data_subtracted_up",
                "data_subtracted_down",
                "QCD",
            ]:
                SRlike_hists["data_subtracted"].Add(SRlike_hists[hist], -1)
                SRlike_hists["data_subtracted_up"].Add(SRlike_hists[hist], -0.93)  # TODO: ask whats this magic numbers?
                SRlike_hists["data_subtracted_down"].Add(SRlike_hists[hist], -1.07)  # Answer: Historical reasons
        for hist in ARlike_hists:
            if hist not in [
                "data",
                "data_subtracted",
                "data_subtracted_up",
                "data_subtracted_down",
                "QCD",
            ]:
                ARlike_hists["data_subtracted"].Add(ARlike_hists[hist], -1)
                ARlike_hists["data_subtracted_up"].Add(ARlike_hists[hist], -0.93)
                ARlike_hists["data_subtracted_down"].Add(ARlike_hists[hist], -1.07)

        # Start of the FF calculation
        FF_hist, FF_hist_up, FF_hist_down = func.calculate_QCD_FF(
            SRlike=SRlike_hists, ARlike=ARlike_hists
        )
        # performing the fit and calculating the fit uncertainties
        fit_graphs, corrlib_exp = func.fit_function(
            ff_hists=[FF_hist.Clone(), FF_hist_up, FF_hist_down],
            bin_edges=process_conf["var_bins"],
            logger=logger,
            fit_option=process_conf.get("fit_option", "poly_1"),
        )

        plotting.plot_FFs(
            variable=process_conf["var_dependence"],
            ff_ratio=FF_hist,
            uncertainties=fit_graphs,
            era=config["era"],
            channel=config["channel"],
            process=process,
            category=split,
            output_path=output_path,
            logger=logger,
            draw_option=process_conf.get("fit_option", "fit")
        )

        if len(split) == 1:
            corrlib_expressions[
                f"{split_variables[0]}#{split[split_variables[0]]}"
            ] = corrlib_exp
        elif len(split) == 2:
            if (
                f"{split_variables[0]}#{split[split_variables[0]]}"
                not in corrlib_expressions
            ):
                corrlib_expressions[
                    f"{split_variables[0]}#{split[split_variables[0]]}"
                ] = dict()
            corrlib_expressions[f"{split_variables[0]}#{split[split_variables[0]]}"][
                f"{split_variables[1]}#{split[split_variables[1]]}"
            ] = corrlib_exp
        else:
            raise Exception("Category splitting is only defined up to 2 dimensions.")

        # producing some control plots
        data = "data"
        if config["use_embedding"]:
            samples = [
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

        plotting.plot_data_mc(
            variable=process_conf["var_dependence"],
            hists=SRlike_hists,
            era=config["era"],
            channel=config["channel"],
            process=process,
            region="SR_like",
            data=data,
            samples=samples,
            category=split,
            output_path=output_path,
            logger=logger,
        )
        plotting.plot_data_mc(
            variable=process_conf["var_dependence"],
            hists=ARlike_hists,
            era=config["era"],
            channel=config["channel"],
            process=process,
            region="AR_like",
            data=data,
            samples=samples,
            category=split,
            output_path=output_path,
            logger=logger,
        )
        log.info("-" * 50)

    return corrlib_expressions


def non_closure_correction(
    config: Dict[str, Union[str, Dict, List]],
    corr_config: Dict[str, Union[str, Dict]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    closure_variable: str,
    evaluator: FakeFactorEvaluator,
    corr_evaluator: FakeFactorCorrectionEvaluator,
    for_DRtoSR: bool,
    logger: str,
    **kwargs: Dict[str, Any],
) -> Dict[str, np.ndarray]:
    """
    This function calculates non closure corrections for fake factors for QCD.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        corr_config: A dictionary with all the relevant information for calculating corrections to the measured fake factors
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        process: This is relevant for QCD because for the tt channel two different QCD fake factors are calculated, one for each hadronic tau
        closure_variable: Name of the variable the correction is calculated for
        evaluator: Evaluator with QCD fake factors
        corr_evaluator: Evaluator with corrections to QCD fake factors
        for_DRtoSR: If True closure correction for the DR to SR correction fake factors will be calculated, if False for the general fake factors
        logger: Name of the logger that should be used

    Return:
        Dictionary of arrays with information about the smoothed function values to be stored with correctionlib (nominal and variations)
    """
    log = logging.getLogger(logger)

    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    # get process specific config information
    process_conf = copy.deepcopy(config["target_processes"][process])
    if for_DRtoSR:
        correction_conf = corr_config["target_processes"][process]["DR_SR"][
            "non_closure"
        ][closure_variable]
    else:
        correction_conf = corr_config["target_processes"][process]["non_closure"][
            closure_variable
        ]

    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        log.info(f"Processing {sample} for the non closure correction for {process}.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for QCD signal-like region
        region_conf = copy.deepcopy(process_conf["SRlike_cuts"])
        rdf_SRlike = func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
            region_cuts=region_conf,
        )

        log.info(
            f"Filtering events for the signal-like region. Target process: {process}"
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # event filter for QCD application-like region
        region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
            region_cuts=region_conf,
        )

        log.info(
            f"Filtering events for the application-like region. Target process: {process}"
        )

        # evaluate the measured fake factors for the specific processes
        if sample == "data":
            if config["channel"] != "tt" or process == "QCD_subleading":
                rdf_ARlike = evaluator.evaluate_subleading_lep_pt_njets(rdf=rdf_ARlike)
            else:
                rdf_ARlike = evaluator.evaluate_leading_lep_pt_njets(rdf=rdf_ARlike)
            if corr_evaluator == None:
                rdf_ARlike = rdf_ARlike.Define(
                    "weight_ff", f"weight * {process}_fake_factor"
                )
            else:
                if config["channel"] != "tt" or process == "QCD_subleading":
                    rdf_ARlike = corr_evaluator.evaluate_leading_lep_pt(rdf=rdf_ARlike)
                else:
                    rdf_ARlike = corr_evaluator.evaluate_subleading_lep_pt(
                        rdf=rdf_ARlike
                    )
                rdf_ARlike = rdf_ARlike.Define(
                    "weight_ff", f"weight * {process}_fake_factor * {process}_ff_corr"
                )

        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # get binning of the dependent variable
        xbinning = array.array("d", correction_conf["var_bins"])
        nbinsx = len(correction_conf["var_bins"]) - 1

        # making the histograms
        h = rdf_SRlike.Histo1D(
            (correction_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            correction_conf["var_dependence"],
            "weight",
        )
        SRlike_hists[sample] = h.GetValue()

        h = rdf_ARlike.Histo1D(
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        ARlike_hists[sample] = h.GetValue()

        if sample == "data":
            h = rdf_ARlike.Histo1D(
                (
                    correction_conf["var_dependence"],
                    f"{sample}_ff",
                    nbinsx,
                    xbinning,
                ),
                correction_conf["var_dependence"],
                "weight_ff",
            )
            ARlike_hists["data_ff"] = h.GetValue()

    SRlike_hists["data_subtracted"] = SRlike_hists["data"].Clone()
    ARlike_hists["data_subtracted"] = ARlike_hists["data"].Clone()

    for hist in SRlike_hists:
        if hist not in ["data", "data_subtracted", "QCD"]:
            SRlike_hists["data_subtracted"].Add(SRlike_hists[hist], -1)
    for hist in ARlike_hists:
        if hist not in ["data", "data_subtracted", "data_ff", "QCD"]:
            ARlike_hists["data_subtracted"].Add(ARlike_hists[hist], -1)

    correction_hist, process_fraction = func.calculate_non_closure_correction(
        SRlike=SRlike_hists,
        ARlike=ARlike_hists,
    )

    smoothed_graph, correction_dict = func.smooth_function(
        hist=correction_hist.Clone(), bin_edges=correction_conf["var_bins"]
    )

    if for_DRtoSR:
        add_str = "_for_DRtoSR"
    else:
        add_str = ""

    plotting.plot_correction(
        variable=correction_conf["var_dependence"],
        corr_hist=correction_hist,
        corr_graph=smoothed_graph,
        corr_name="non_closure_" + closure_variable + add_str,
        era=config["era"],
        channel=config["channel"],
        process=process,
        output_path=output_path,
        logger=logger,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SRlike_hists["data_subtracted"].Clone()
    plot_hists["data_ff"] = ARlike_hists["data_ff"].Clone()
    plot_hists["data_ff"].Scale(process_fraction)

    data = "data_subtracted"
    samples = ["data_ff"]

    plotting.plot_data_mc(
        variable=correction_conf["var_dependence"],
        hists=plot_hists,
        era=config["era"],
        channel=config["channel"],
        process=process,
        region="non_closure_" + closure_variable + add_str,
        data=data,
        samples=samples,
        category={"incl": ""},
        output_path=output_path,
        logger=logger,
    )

    data = "data"
    if config["use_embedding"]:
        samples = [
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

    plotting.plot_data_mc(
        variable=correction_conf["var_dependence"],
        hists=SRlike_hists,
        era=config["era"],
        channel=config["channel"],
        process=process,
        region="non_closure_" + closure_variable + add_str + "_SRlike_hist",
        data=data,
        samples=samples,
        category={"incl": ""},
        output_path=output_path,
        logger=logger,
    )

    return correction_dict


def DR_SR_correction(
    config: Dict[str, Union[str, Dict, List]],
    corr_config: Dict[str, Union[str, Dict]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    evaluator: FakeFactorEvaluator,
    corr_evaluator: FakeFactorCorrectionEvaluator,
    logger: str,
) -> Dict[str, np.ndarray]:
    """
    This function calculates DR to SR corrections for fake factors for QCD.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        corr_config: A dictionary with all the relevant information for calculating corrections to the measured fake factors
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        process: This is relevant for QCD because for the tt channel two different QCD fake factors are calculated, one for each hadronic tau
        evaluator: Evaluator with QCD fake factors
        corr_evaluator: Evaluator with corrections to QCD fake factors
        logger: Name of the logger that should be used

    Return:
        Dictionary of arrays with information about the smoothed function values to be stored with correctionlib (nominal and variations)
    """
    log = logging.getLogger(logger)

    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    # get process specific config information
    process_conf = copy.deepcopy(config["target_processes"][process])
    correction_conf = corr_config["target_processes"][process]["DR_SR"]

    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        log.info(f"Processing {sample} for the DR to SR correction for {process}.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for QCD signal-like region
        region_conf = copy.deepcopy(process_conf["SRlike_cuts"])
        rdf_SRlike = func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
            region_cuts=region_conf,
        )

        log.info(
            f"Filtering events for the signal-like region. Target process: {process}"
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # event filter for QCD application-like region
        region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
            region_cuts=region_conf,
        )

        log.info(
            f"Filtering events for the application-like region. Target process: {process}"
        )

        # evaluate the measured fake factors for the specific processes
        if sample == "data":
            if config["channel"] != "tt" or process == "QCD_subleading":
                rdf_ARlike = evaluator.evaluate_subleading_lep_pt_njets(rdf=rdf_ARlike)
                rdf_ARlike = corr_evaluator.evaluate_leading_lep_pt(rdf=rdf_ARlike)
            else:
                rdf_ARlike = evaluator.evaluate_leading_lep_pt_njets(rdf=rdf_ARlike)
                rdf_ARlike = corr_evaluator.evaluate_subleading_lep_pt(rdf=rdf_ARlike)
            rdf_ARlike = rdf_ARlike.Define(
                "weight_ff", f"weight * {process}_fake_factor * {process}_ff_corr"
            )

        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # get binning of the dependent variable
        xbinning = array.array("d", correction_conf["var_bins"])
        nbinsx = len(correction_conf["var_bins"]) - 1

        # making the histograms
        h = rdf_SRlike.Histo1D(
            (correction_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            correction_conf["var_dependence"],
            "weight",
        )
        SRlike_hists[sample] = h.GetValue()

        h = rdf_ARlike.Histo1D(
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        ARlike_hists[sample] = h.GetValue()

        if sample == "data":
            h = rdf_ARlike.Histo1D(
                (
                    correction_conf["var_dependence"],
                    f"{sample}_ff",
                    nbinsx,
                    xbinning,
                ),
                correction_conf["var_dependence"],
                "weight_ff",
            )
            ARlike_hists["data_ff"] = h.GetValue()

    SRlike_hists["data_subtracted"] = SRlike_hists["data"].Clone()
    ARlike_hists["data_subtracted"] = ARlike_hists["data"].Clone()

    for hist in SRlike_hists:
        if hist not in ["data", "data_subtracted", "QCD"]:
            SRlike_hists["data_subtracted"].Add(SRlike_hists[hist], -1)
    for hist in ARlike_hists:
        if hist not in ["data", "data_subtracted", "data_ff", "QCD"]:
            ARlike_hists["data_subtracted"].Add(ARlike_hists[hist], -1)

    correction_hist, process_fraction = func.calculate_non_closure_correction(
        SRlike=SRlike_hists,
        ARlike=ARlike_hists,
    )

    smoothed_graph, correction_dict = func.smooth_function(
        hist=correction_hist.Clone(),
        bin_edges=correction_conf["var_bins"],
    )

    plotting.plot_correction(
        variable=correction_conf["var_dependence"],
        corr_hist=correction_hist,
        corr_graph=smoothed_graph,
        corr_name="DR_SR",
        era=config["era"],
        channel=config["channel"],
        process=process,
        output_path=output_path,
        logger=logger,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SRlike_hists["data_subtracted"].Clone()
    plot_hists["data_ff"] = ARlike_hists["data_ff"].Clone()
    plot_hists["data_ff"].Scale(process_fraction)

    data = "data_subtracted"
    samples = ["data_ff"]

    plotting.plot_data_mc(
        variable=correction_conf["var_dependence"],
        hists=plot_hists,
        era=config["era"],
        channel=config["channel"],
        process=process,
        region="DR_SR",
        data=data,
        samples=samples,
        category={"incl": ""},
        output_path=output_path,
        logger=logger,
    )

    data = "data"
    if config["use_embedding"]:
        samples = [
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

    plotting.plot_data_mc(
        variable=correction_conf["var_dependence"],
        hists=SRlike_hists,
        era=config["era"],
        channel=config["channel"],
        process=process,
        region="DR_SR" + "_SRlike_hist",
        data=data,
        samples=samples,
        category={"incl": ""},
        output_path=output_path,
        logger=logger,
    )

    return correction_dict
