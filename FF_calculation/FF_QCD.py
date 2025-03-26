"""
Function for calculating fake factors for the QCD process
"""

import array
import copy
import logging
from io import StringIO
from typing import Any, Dict, Tuple, Union

import numpy as np
import ROOT
from wurlitzer import STDOUT, pipes

import configs.general_definitions as gd
import helper.ff_functions as ff_func
import helper.plotting as plotting
from helper.functions import RuntimeVariables


def calculation_QCD_FFs(
    args: Tuple[Any, ...],
) -> Dict[str, Union[str, Dict[str, str]]]:
    """
    This function calculates fake factors for the QCD process for a specific category (split).

    Args:
        args: Tuple containing all the necessary information for the calculation of the fake factors
            splitting: SplitQuantitiesContainer, contains the splitting information
            config: Dictionary with all the relevant information for the fake factor calculation
            process_conf: Dictionary with all the relevant information for the fake factor calculation of the specific process
            process: Name of the process
            sample_paths: List of file paths where the samples are stored
            output_path: Path where the generated plots should be stored
            logger: Name of the logger that should be used

    Return:
        Dictionary with the category information as keys and the fitted functions (including variations) as values
    """
    (
        splitting,  # splitting: Dict[str, str],
        config,  # config: Dict[str, Union[str, Dict, List]],
        process_conf,  # process_conf: Dict[str, Union[str, Dict, List]],
        process,  # process: str,
        sample_paths,  # sample_paths: List[str],
        output_path,  # output_path: str,
        logger,  # logger: str,
        *_,  # SRlike_hists, ARlike_hists only used in ttbar calculation
    ) = args

    log = logging.getLogger(logger)

    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        log.info(f"Processing {sample} for the {', '.join(['{} {}'.format(var, splitting.split[var]) for var in splitting.variables])} category.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for QCD signal-like region
        region_conf = copy.deepcopy(process_conf["SRlike_cuts"])
        rdf_SRlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for the signal-like region. Target process: {process}")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # event filter for QCD application-like region
        region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for the application-like region. Target process: {process}")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # get binning of the dependent variable
        xbinning = array.array("d", splitting.var_bins)
        nbinsx = len(splitting.var_bins) - 1

        # making the histograms
        h = RuntimeVariables.RDataFrameWrapper(rdf_SRlike).Histo1D(
            (process_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            process_conf["var_dependence"],
            "weight",
        )
        SRlike_hists[sample] = h.GetValue()

        h = RuntimeVariables.RDataFrameWrapper(rdf_ARlike).Histo1D(
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
    FF_hist, FF_hist_up, FF_hist_down = ff_func.calculate_QCD_FF(
        SRlike=SRlike_hists,
        ARlike=ARlike_hists,
    )
    # performing the fit and calculating the fit uncertainties
    nominal_draw_obj, fit_graphs, corrlib_exp, used_fit = ff_func.fit_function(
        ff_hists=[FF_hist.Clone(), FF_hist_up, FF_hist_down],
        bin_edges=splitting.var_bins,
        logger=logger,
        fit_option=splitting.fit_option,
        limit_kwargs=splitting.limit_kwargs(hist=FF_hist),
    )

    plotting.plot_FFs(
        variable=process_conf["var_dependence"],
        ff_ratio=nominal_draw_obj,
        uncertainties=fit_graphs,
        era=config["era"],
        channel=config["channel"],
        process=process,
        category=splitting.split,
        output_path=output_path,
        logger=logger,
        draw_option=used_fit,
        save_data=True,
    )

    # producing some control plots
    for _hist, _region in [
        (SRlike_hists, "SR_like"),
        (ARlike_hists, "AR_like"),
    ]:
        for yscale, save_data in zip(["linear", "log"], [True, False]):
            plotting.plot_data_mc(
                variable=process_conf["var_dependence"],
                hists=_hist,
                era=config["era"],
                channel=config["channel"],
                process=process,
                region=_region,
                data="data",
                samples=ff_func.controlplot_samples(config["use_embedding"], add_qcd=False),
                category=splitting.split,
                output_path=output_path,
                logger=logger,
                yscale=yscale,
                save_data=save_data,
            )
    log.info("-" * 50)

    return ff_func.fill_corrlib_expression(corrlib_exp, splitting.variables, splitting.split)


def non_closure_correction(
    args: Tuple[Any, ...],
) -> Dict[str, np.ndarray]:
    """
    This function calculates non closure corrections for fake factors for QCD.

    Intended to be used in a multiprocessing environment.

    Args:
        args: Tuple containing all the necessary information for the calculation of the non closure corrections
            splitting: SplitQuantitiesContainer, contains the splitting information
            config: Dictionary with all the relevant information for the fake factor calculation
            correction_conf: Dictionary with all the relevant information for the correction calculation
            process: Name of the process
            closure_variable: Name of the variable the correction is calculated for
            sample_paths: List of file paths where the samples are stored
            output_path: Path where the generated plots should be stored
            logger: Name of the logger that should be used
            evaluator: Evaluator with QCD fake factors
            corr_evaluators: List of evaluators with corrections to QCD fake factors
            for_DRtoSR: If True closure correction for the DR to SR correction fake factors will be calculated, if False for the general fake factors
    """
    (
        splitting,
        config,
        correction_conf,
        process,
        closure_variable,
        sample_paths,
        output_path,
        logger,
        evaluator,
        corr_evaluators,
        for_DRtoSR,
    ) = args

    log = logging.getLogger(logger)

    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        if splitting.split is not None:
            log.info(f"Processing {sample} for the non closure correction for {process} for {', '.join([f'{var} {splitting.split[var]}' for var in splitting.variables])}.")
        else:
            log.info(f"Processing {sample} for the non closure correction for {process}.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for QCD signal-like region
        region_conf = copy.deepcopy(config["target_processes"][process]["SRlike_cuts"])
        rdf_SRlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for the signal-like region. Target process: {process}")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # event filter for QCD application-like region
        region_conf = copy.deepcopy(config["target_processes"][process]["ARlike_cuts"])
        rdf_ARlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for the application-like region. Target process: {process}")

        # evaluate the measured fake factors for the specific processes
        if sample == "data":
            rdf_ARlike = evaluator.evaluate_fake_factor(rdf=rdf_ARlike)

            # additionally evaluate the previous corrections
            corr_str = ""
            for corr_evaluator in corr_evaluators:
                rdf_ARlike = corr_evaluator.evaluate_correction(rdf=rdf_ARlike)
                corr_str += f" * {process}_ff_corr_{corr_evaluator.variable}"

            rdf_ARlike = rdf_ARlike.Define(
                "weight_ff",
                f"weight * {process}_fake_factor{corr_str}",
            )

        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # get binning of the dependent variable
        xbinning, nbinsx = array.array("d", splitting.var_bins), len(splitting.var_bins) - 1

        # making the histograms
        h = RuntimeVariables.RDataFrameWrapper(rdf_SRlike).Histo1D(
            (correction_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            correction_conf["var_dependence"],
            "weight",
        )
        SRlike_hists[sample] = h.GetValue()

        h = RuntimeVariables.RDataFrameWrapper(rdf_ARlike).Histo1D(
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        ARlike_hists[sample] = h.GetValue()

        if sample == "data":
            h = RuntimeVariables.RDataFrameWrapper(rdf_ARlike).Histo1D(
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

    correction_hist, process_fraction = ff_func.calculate_non_closure_correction(
        SRlike=SRlike_hists,
        ARlike=ARlike_hists,
    )

    nominal_draw_obj, smoothed_graph, correction_dict = ff_func.smooth_function(
        hist=correction_hist.Clone(),
        bin_edges=splitting.var_bins,
        write_corrections=splitting.write_corrections,
        bandwidth=splitting.bandwidth,
    )

    add_str = "_for_DRtoSR" if for_DRtoSR else ""

    plotting.plot_correction(
        variable=correction_conf["var_dependence"],
        corr_hist=nominal_draw_obj,
        corr_graph=smoothed_graph,
        corr_name=f"non_closure_{closure_variable}{add_str}",
        era=config["era"],
        channel=config["channel"],
        process=process,
        output_path=output_path,
        category=splitting.split,
        logger=logger,
        save_data=True,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SRlike_hists["data_subtracted"].Clone()
    plot_hists["data_ff"] = ARlike_hists["data_ff"].Clone()
    plot_hists["data_ff"].Scale(process_fraction)

    for yscale in ["linear", "log"]:
        plotting.plot_data_mc(
            variable=correction_conf["var_dependence"],
            hists=plot_hists,
            era=config["era"],
            channel=config["channel"],
            process=process,
            region=f"non_closure_{closure_variable}{add_str}",
            data="data_subtracted",
            samples=["data_ff"],
            category=splitting.split or {"incl": ""},
            output_path=output_path,
            logger=logger,
            yscale=yscale,
        )

    # producing some control plots
    for yscale, save_data in zip(["linear", "log"], [True, False]):
        plotting.plot_data_mc(
            variable=correction_conf["var_dependence"],
            hists=SRlike_hists,
            era=config["era"],
            channel=config["channel"],
            process=process,
            region=f"non_closure_{closure_variable}{add_str}_SRlike_hist",
            data="data",
            samples=ff_func.controlplot_samples(config["use_embedding"], add_qcd=False),
            category=splitting.split or {"incl": ""},
            output_path=output_path,
            logger=logger,
            yscale=yscale,
            save_data=save_data,
        )

    if splitting.split is not None:
        return ff_func.fill_corrlib_expression(correction_dict, splitting.variables, splitting.split)
    else:
        return correction_dict


def DR_SR_correction(
    args: Tuple[Any, ...],
) -> Dict[str, np.ndarray]:
    """
    This function calculates DR to SR correction for fake factors for QCD.

    Intended to be used in a multiprocessing environment.

    Args:
        args: Tuple containing all the necessary information for the calculation of the DR to SR correction
            splitting: SplitQuantitiesContainer, contains the splitting information
            config: Dictionary with all the relevant information for the fake factor calculation
            correction_conf: Dictionary with all the relevant information for the correction calculation
            process: Name of the process
            sample_paths: List of file paths where the samples are stored
            output_path: Path where the generated plots should be stored
            logger: Name of the logger that should be used
            evaluator: Evaluator with QCD fake factors
            corr_evaluators: List of evaluators with corrections to QCD fake factors
    """
    (
        splitting,
        config,
        correction_conf,
        process,
        sample_paths,
        output_path,
        logger,
        evaluator,
        corr_evaluators,
    ) = args

    log = logging.getLogger(logger)

    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        if splitting.split is not None:
            log.info(f"Processing {sample} for the DR to SR correction for {process} for {', '.join([f'{var} {splitting.split[var]}' for var in splitting.variables])}.")
        else:
            log.info(f"Processing {sample} for the DR to SR correction for {process}.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for QCD signal-like region
        region_conf = copy.deepcopy(config["target_processes"][process]["SRlike_cuts"])
        rdf_SRlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for the signal-like region. Target process: {process}")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # event filter for QCD application-like region
        region_conf = copy.deepcopy(config["target_processes"][process]["ARlike_cuts"])
        rdf_ARlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for the application-like region. Target process: {process}")

        # evaluate the measured fake factors for the specific processes
        if sample == "data":
            rdf_ARlike = evaluator.evaluate_fake_factor(rdf=rdf_ARlike)

            # additionally evaluate the previous corrections
            corr_str = ""
            for corr_evaluator in corr_evaluators:
                rdf_ARlike = corr_evaluator.evaluate_correction(rdf=rdf_ARlike)
                corr_str += f" * {process}_ff_corr_{corr_evaluator.variable}"

            rdf_ARlike = rdf_ARlike.Define(
                "weight_ff",
                f"weight * {process}_fake_factor{corr_str}",
            )

        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # get binning of the dependent variable
        xbinning, nbinsx = array.array("d", splitting.var_bins), len(splitting.var_bins) - 1

        # making the histograms
        h = RuntimeVariables.RDataFrameWrapper(rdf_SRlike).Histo1D(
            (correction_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            correction_conf["var_dependence"],
            "weight",
        )
        SRlike_hists[sample] = h.GetValue()

        h = RuntimeVariables.RDataFrameWrapper(rdf_ARlike).Histo1D(
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        ARlike_hists[sample] = h.GetValue()

        if sample == "data":
            h = RuntimeVariables.RDataFrameWrapper(rdf_ARlike).Histo1D(
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

    correction_hist, process_fraction = ff_func.calculate_non_closure_correction(
        SRlike=SRlike_hists,
        ARlike=ARlike_hists,
    )

    nominal_draw_obj, smoothed_graph, correction_dict = ff_func.smooth_function(
        hist=correction_hist.Clone(),
        bin_edges=splitting.var_bins,
        write_corrections=splitting.write_corrections,
        bandwidth=splitting.bandwidth,
    )

    plotting.plot_correction(
        variable=correction_conf["var_dependence"],
        corr_hist=nominal_draw_obj,
        corr_graph=smoothed_graph,
        corr_name="DR_SR",
        era=config["era"],
        channel=config["channel"],
        process=process,
        output_path=output_path,
        category=splitting.split,
        logger=logger,
        save_data=True,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SRlike_hists["data_subtracted"].Clone()
    plot_hists["data_ff"] = ARlike_hists["data_ff"].Clone()
    plot_hists["data_ff"].Scale(process_fraction)

    for yscale in ["linear", "log"]:
        plotting.plot_data_mc(
            variable=correction_conf["var_dependence"],
            hists=plot_hists,
            era=config["era"],
            channel=config["channel"],
            process=process,
            region="DR_SR",
            data="data_subtracted",
            samples=["data_ff"],
            category=splitting.split or {"incl": ""},
            output_path=output_path,
            logger=logger,
            yscale=yscale,
        )

    for yscale, save_data in zip(["linear", "log"], [True, False]):
        plotting.plot_data_mc(
            variable=correction_conf["var_dependence"],
            hists=SRlike_hists,
            era=config["era"],
            channel=config["channel"],
            process=process,
            region="DR_SR" + "_SRlike_hist",
            data="data",
            samples=ff_func.controlplot_samples(config["use_embedding"], add_qcd=False),
            category=splitting.split or {"incl": ""},
            output_path=output_path,
            logger=logger,
            yscale=yscale,
            save_data=save_data,
        )

    if splitting.split is not None:
        return ff_func.fill_corrlib_expression(correction_dict, splitting.variables, splitting.split)
    else:
        return correction_dict
