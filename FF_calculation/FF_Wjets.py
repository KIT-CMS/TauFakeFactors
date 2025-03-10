"""
Function for calculating fake factors for the W-jets process
"""

import array
import copy
import logging

# from collections import defaultdict
from io import StringIO
from typing import Any, Dict, List, Union, Tuple

import numpy as np
import ROOT
from wurlitzer import STDOUT, pipes

import helper.functions as func
import helper.ff_functions as ff_func
import helper.plotting as plotting
from helper.ff_evaluators import FakeFactorCorrectionEvaluator, FakeFactorEvaluator
import configs.general_definitions as gd


def controlplot_samples(
    config: Dict[str, Union[str, Dict, List]],
) -> List[str]:
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
    ]
    if config["use_embedding"]:
        samples.append("embedding")
    else:
        samples.extend(["diboson_T", "ttbar_T", "DYjets_T", "ST_T"])

    return samples


def _split_calculation_Wjets_FFs(
    args: Tuple[Any, ...],
) -> Dict[str, Union[Dict[str, str], Dict[str, Dict[str, str]]]]:
    """
    This function calculates fake factors for the WJet process for a specific category.

    Intended to be used in a multiprocessing environment.

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
        Dictionary with the category information as keys and the fitted functions (including variations) as values
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
    ) = args

    log = logging.getLogger(logger)

    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()
    # init histogram dict for QCD SS/OS estimation
    SRlike_hists_qcd = dict()
    ARlike_hists_qcd = dict()

    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        log.info(f"Processing {sample} for the {', '.join([f'{var} {split[var]}' for var in split_variables])} category.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for Wjets signal-like region
        region_conf = copy.deepcopy(process_conf["SRlike_cuts"])
        rdf_SRlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=split,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for the signal-like region. Target process: {process}")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # QCD estimation from same sign in signal-like region
        if "tau_pair_sign" in region_conf:
            region_conf["tau_pair_sign"] = "(q_1*q_2) > 0"  # same sign
        else:
            raise ValueError(f"No tau pair sign cut defined in the {process} config. Is needed for the QCD estimation.")

        rdf_SRlike_qcd = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=split,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for QCD estimation in the signal-like region. Target process: {process}")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike_qcd.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # event filter for Wjets application-like region
        region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=split,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for the application-like region. Target process: {process}")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # QCD estimation from same sign in application-like region
        if "tau_pair_sign" in region_conf:
            region_conf["tau_pair_sign"] = "(q_1*q_2) > 0"  # same sign
        else:
            raise ValueError(
                f"No tau pair sign cut defined in the {process} config. Is needed for the QCD estimation."
            )

        rdf_ARlike_qcd = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=split,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for QCD estimation in the application-like region. Target process: {process}")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike_qcd.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # get binning of the dependent variable from the computed binning
        xbinning = array.array("d", binning)
        nbinsx = len(binning) - 1

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

        # making the histograms for QCD estimation
        h_qcd = rdf_SRlike_qcd.Histo1D(
            (process_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            process_conf["var_dependence"],
            "weight",
        )
        SRlike_hists_qcd[sample] = h_qcd.GetValue()

        h_qcd = rdf_ARlike_qcd.Histo1D(
            (process_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            process_conf["var_dependence"],
            "weight",
        )
        ARlike_hists_qcd[sample] = h_qcd.GetValue()

    # calculate QCD estimation
    SRlike_hists["QCD"] = ff_func.QCD_SS_estimate(hists=SRlike_hists_qcd)
    ARlike_hists["QCD"] = ff_func.QCD_SS_estimate(hists=ARlike_hists_qcd)

    # calculate Wjets enriched data by subtraction all there backgrould sample
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
            "Wjets",
        ]:
            SRlike_hists["data_subtracted"].Add(SRlike_hists[hist], -1)
            SRlike_hists["data_subtracted_up"].Add(SRlike_hists[hist], -0.93)
            SRlike_hists["data_subtracted_down"].Add(SRlike_hists[hist], -1.07)
    for hist in ARlike_hists:
        if hist not in [
            "data",
            "data_subtracted",
            "data_subtracted_up",
            "data_subtracted_down",
            "Wjets",
        ]:
            ARlike_hists["data_subtracted"].Add(ARlike_hists[hist], -1)
            ARlike_hists["data_subtracted_up"].Add(ARlike_hists[hist], -0.93)
            ARlike_hists["data_subtracted_down"].Add(ARlike_hists[hist], -1.07)

    # Start of the FF calculation
    FF_hist, FF_hist_up, FF_hist_down = ff_func.calculate_Wjets_FF(
        SRlike=SRlike_hists, ARlike=ARlike_hists
    )
    # performing the fit and calculating the uncertainties
    fit_graphs, corrlib_exp, used_fit = ff_func.fit_function(
        ff_hists=[FF_hist.Clone(), FF_hist_up, FF_hist_down],
        bin_edges=binning,
        logger=logger,
        fit_option=process_conf.get("fit_option", gd.default_fit_options["Wjets"]),
        limit_kwargs=process_conf.get(
            "limit_kwargs",
            gd.get_default_fit_function_limit_kwargs(binning),
        ),
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
        draw_option=used_fit,
        save_data=True,
    )

    # producing some control plots
    for _hist, _region, _data, _samples in [
        (SRlike_hists, "SR_like", "data", controlplot_samples(config)),
        (ARlike_hists, "AR_like", "data", controlplot_samples(config)),
        (SRlike_hists, "SR_like", "data_subtracted", ["Wjets"]),
        (ARlike_hists, "AR_like", "data_subtracted", ["Wjets"]),
    ]:
        for yscale, save_data in zip(["linear", "log"], [True, False]):
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
                yscale=yscale,
                save_data=save_data,
            )
    log.info("-" * 50)

    return ff_func.fill_corrlib_expression(corrlib_exp, split_variables, split)


def calculation_Wjets_FFs(
    config: Dict[str, Union[str, Dict, List]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    logger: str,
) -> Dict[str, Union[Dict[str, str], Dict[str, Dict[str, str]]]]:
    """
    This function calculates fake factors for W+jets.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        logger: Name of the logger that should be used

    Return:
        Dictionary where the categories are defined as keys and and the values are the fitted functions (including variations)
        e.g. corrlib_expressions[CATEGORY_1][CATEGORY_2][VARIATION] if dimension of categories is 2
    """

    process_conf = config["target_processes"][process]

    split_variables, split_combinations, split_binnings = ff_func.get_split_combinations(
        categories=process_conf["split_categories"],
        binning=process_conf["var_bins"],
    )

    assert len(split_variables) < 3, "Category splitting is only defined up to 2 dimensions."

    results = func.optional_process_pool(
        args_list=[
            (
                split,
                binning,
                config,
                process_conf,
                process,
                split_variables,
                sample_paths,
                output_path,
                logger,
            )
            for split, binning in zip(split_combinations, split_binnings)
        ],
        function=_split_calculation_Wjets_FFs,
    )

    return ff_func.fill_corrlib_expression(results, split_variables)


def non_closure_correction(
    config: Dict[str, Union[str, Dict, List]],
    corr_config: Dict[str, Union[str, Dict]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    closure_variable: str,
    evaluator: FakeFactorEvaluator,
    corr_evaluators: List[FakeFactorCorrectionEvaluator],
    for_DRtoSR: bool,
    logger: str,
) -> Dict[str, np.ndarray]:
    """
    This function calculates non closure corrections for fake factors for W+jets.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        corr_config: A dictionary with all the relevant information for calculating corrections to the measured fake factors
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        process: This is relevant for QCD because for the tt channel two different QCD fake factors are calculated, one for each hadronic tau
        closure_variable: Name of the variable the correction is calculated for
        evaluator: Evaluator with Wjets fake factors
        corr_evaluators: List of evaluator with corrections to Wjets fake factors
        for_DRtoSR: If True closure correction for the DR to SR correction fake factors will be calculated, if False for the general fake factors
        logger: Name of the logger that should be used

    Return:
        Dictionary of arrays with information about the smoothed function values to be stored with correctionlib (nominal and variations)
    """
    log = logging.getLogger(logger)

    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    # init histogram dict for QCD SS/OS estimation
    SRlike_hists_qcd = dict()
    ARlike_hists_qcd = dict()

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

        # event filter for Wjets signal-like region
        region_conf = copy.deepcopy(process_conf["SRlike_cuts"])
        rdf_SRlike = ff_func.apply_region_filters(
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

        # QCD estimation from same sign in signal-like region
        if "tau_pair_sign" in region_conf:
            region_conf["tau_pair_sign"] = "(q_1*q_2) > 0"  # same sign
        else:
            raise ValueError(
                f"No tau pair sign cut defined in the {process} config. Is needed for the QCD estimation."
            )

        rdf_SRlike_qcd = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
            region_cuts=region_conf,
        )

        log.info(
            f"Filtering events for QCD estimation in the signal-like region. Target process: {process}"
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike_qcd.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # event filter for Wjets application-like region
        region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = ff_func.apply_region_filters(
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
            rdf_ARlike = evaluator.evaluate_fake_factor(rdf=rdf_ARlike)

            # additionally evaluate the previous corrections
            corr_str = ""
            for corr_evaluator in corr_evaluators:
                rdf_ARlike = corr_evaluator.evaluate_correction(
                    rdf=rdf_ARlike,
                )
                corr_str += f" * {process}_ff_corr_{corr_evaluator.variable}"

            rdf_ARlike = rdf_ARlike.Define(
                "weight_ff", f"weight * {process}_fake_factor{corr_str}"
            )

        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # QCD estimation from same sign in application-like region
        if "tau_pair_sign" in region_conf:
            region_conf["tau_pair_sign"] = "(q_1*q_2) > 0"  # same sign
        else:
            raise ValueError(
                f"No tau pair sign cut defined in the {process} config. Is needed for the QCD estimation."
            )

        rdf_ARlike_qcd = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
            region_cuts=region_conf,
        )

        log.info(
            f"Filtering events for QCD estimation in the application-like region. Target process: {process}"
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike_qcd.Report().Print()
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

        # making the histograms for QCD estimation
        h_qcd = rdf_SRlike_qcd.Histo1D(
            (correction_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            correction_conf["var_dependence"],
            "weight",
        )
        SRlike_hists_qcd[sample] = h_qcd.GetValue()

        h_qcd = rdf_ARlike_qcd.Histo1D(
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        ARlike_hists_qcd[sample] = h_qcd.GetValue()

    # calculate QCD estimation
    SRlike_hists["QCD"] = ff_func.QCD_SS_estimate(hists=SRlike_hists_qcd)
    ARlike_hists["QCD"] = ff_func.QCD_SS_estimate(hists=ARlike_hists_qcd)

    SRlike_hists["data_subtracted"] = SRlike_hists["data"].Clone()
    ARlike_hists["data_subtracted"] = ARlike_hists["data"].Clone()

    for hist in SRlike_hists:
        if hist not in ["data", "data_subtracted", "Wjets"]:
            SRlike_hists["data_subtracted"].Add(SRlike_hists[hist], -1)
    for hist in ARlike_hists:
        if hist not in ["data", "data_subtracted", "data_ff", "Wjets"]:
            ARlike_hists["data_subtracted"].Add(ARlike_hists[hist], -1)

    correction_hist, process_fraction = ff_func.calculate_non_closure_correction(
        SRlike=SRlike_hists,
        ARlike=ARlike_hists,
    )

    smoothed_graph, correction_dict = ff_func.smooth_function(
        hist=correction_hist.Clone(),
        bin_edges=correction_conf["var_bins"],
        write_corrections=correction_conf["write_corrections"],
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

    for yscale in ["linear", "log"]:
        plotting.plot_data_mc(
            variable=correction_conf["var_dependence"],
            hists=plot_hists,
            era=config["era"],
            channel=config["channel"],
            process=process,
            region="non_closure_" + closure_variable + add_str,
            data="data_subtracted",
            samples=["data_ff"],
            category={"incl": ""},
            output_path=output_path,
            logger=logger,
            yscale=yscale,
        )

    for yscale in ["linear", "log"]:
        plotting.plot_data_mc(
            variable=correction_conf["var_dependence"],
            hists=SRlike_hists,
            era=config["era"],
            channel=config["channel"],
            process=process,
            region="non_closure_" + closure_variable + add_str + "_SRlike",
            data="data",
            samples=controlplot_samples(config),
            category={"incl": ""},
            output_path=output_path,
            logger=logger,
            yscale=yscale,
        )

    return correction_dict


def DR_SR_correction(
    config: Dict[str, Union[str, Dict, List]],
    corr_config: Dict[str, Union[str, Dict]],
    sample_paths: List[str],
    output_path: str,
    evaluator: FakeFactorEvaluator,
    corr_evaluator: FakeFactorCorrectionEvaluator,
    logger: str,
    **kwargs,
) -> Dict[str, np.ndarray]:
    """
    This function calculates DR to SR corrections for fake factors for W+jets.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        corr_config: A dictionary with all the relevant information for calculating corrections to the measured fake factors
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        evaluator: Evaluator with Wjets fake factors
        corr_evaluator: Evaluator with corrections to Wjets fake factors
        logger: Name of the logger that should be used

    Return:
        Dictionary of arrays with information about the smoothed function values to be stored with correctionlib (nominal and variations)
    """
    log = logging.getLogger(logger)

    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    # get process specific config information
    process_conf = copy.deepcopy(config["target_processes"]["Wjets"])
    correction_conf = corr_config["target_processes"]["Wjets"]["DR_SR"]

    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        if sample == "Wjets":
            log.info(f"Processing {sample} for the DR to SR correction for Wjets.")
            log.info("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)

            # event filter for Wjets signal-like region
            region_conf = copy.deepcopy(process_conf["SRlike_cuts"])
            rdf_SRlike = ff_func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=None,
                region_cuts=region_conf,
            )

            log.info(
                "Filtering events for the signal-like region. Target process: Wjets"
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_SRlike.Report().Print()
            log.info(out.getvalue())
            log.info("-" * 50)

            # event filter for Wjets application-like region
            region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
            rdf_ARlike = ff_func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=None,
                region_cuts=region_conf,
            )

            log.info(
                "Filtering events for the application-like region. Target process: Wjets"
            )

            rdf_ARlike = evaluator.evaluate_fake_factor(rdf=rdf_ARlike)
            rdf_ARlike = corr_evaluator.evaluate_correction(rdf=rdf_ARlike)
            rdf_ARlike = rdf_ARlike.Define(
                "weight_ff",
                f"weight * Wjets_fake_factor * Wjets_ff_corr_{corr_evaluator.variable}",
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
                (
                    correction_conf["var_dependence"],
                    f"{sample}",
                    nbinsx,
                    xbinning,
                ),
                correction_conf["var_dependence"],
                "weight",
            )
            SRlike_hists[sample] = h.GetValue()

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
            ARlike_hists["Wjets_ff"] = h.GetValue()

    correction_hist = ff_func.calculate_non_closure_correction_Wjets_fromMC(
        SRlike=SRlike_hists, ARlike=ARlike_hists
    )

    smoothed_graph, correction_dict = ff_func.smooth_function(
        hist=correction_hist.Clone(),
        bin_edges=correction_conf["var_bins"],
        write_corrections=correction_conf["write_corrections"],
    )

    plotting.plot_correction(
        variable=correction_conf["var_dependence"],
        corr_hist=correction_hist,
        corr_graph=smoothed_graph,
        corr_name="DR_SR",
        era=config["era"],
        channel=config["channel"],
        process="Wjets",
        output_path=output_path,
        logger=logger,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SRlike_hists["Wjets"].Clone()
    plot_hists["data_ff"] = ARlike_hists["Wjets_ff"].Clone()

    for yscale in ["linear", "log"]:
        plotting.plot_data_mc(
            variable=correction_conf["var_dependence"],
            hists=plot_hists,
            era=config["era"],
            channel=config["channel"],
            process="Wjets",
            region="DR_SR",
            data="data_subtracted",
            samples=["data_ff"],
            category={"incl": ""},
            output_path=output_path,
            logger=logger,
            yscale=yscale,
        )

    return correction_dict
