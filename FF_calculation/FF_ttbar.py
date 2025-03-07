"""
Function for calculating fake factors for the ttbar process
"""

import concurrent.futures
import array
import copy
import logging
from io import StringIO
from typing import Any, Dict, List, Union, Tuple

import numpy as np
import ROOT
from wurlitzer import STDOUT, pipes

import helper.ff_functions as func
import helper.plotting as plotting
from helper.ff_evaluators import FakeFactorCorrectionEvaluator, FakeFactorEvaluator
import configs.general_definitions as gd


def _split_calculation_ttbar_FFs(
    args: Tuple[Any, ...],
) -> Dict[str, Union[str, Dict[str, str]]]:
    """
    This function calculates fake factors for ttbar for a given category.

    Intded to be used in a multiprocessing environment.

    Args:
        args: Tuple of arguments that are passed to the function
            split: Dictionary containing the category information
            binning: List of bin edges for the dependent variable
            config: Dictionary with all the relevant information for the fake factor calculation
            process_conf: Dictionary with all the relevant information for the fake factor calculation of the specific process
            process: Name of the process
            split_variables: List of variables that are used for the category splitting
            sample_paths: List of file paths where the samples are stored
            output_path: Path where the generated plots should be stored
            logger: Name of the logger that should be used
            SRlike_hists: Dictionary containing histograms for the signal-like region
            ARlike_hists: Dictionary containing histograms for the application-like region

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
        SRlike_hists,  # SRlike_hists: Dict[str, ROOT.TH1D],
        ARlike_hists,  # ARlike_hists: Dict[str, ROOT.TH1D],
    ) = args

    log = logging.getLogger(logger)

    # init histogram dict for FF measurement from MC
    SR_hists = dict()
    AR_hists = dict()

    corrlib_expression = dict()

    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        # FFs for ttbar from mc -> only ttbar with true misindentified jets relevant
        if sample in ["ttbar_J"]:
            log.info(f"Processing {sample} for the {', '.join([f'{var} {split[var]}' for var in split_variables])} category.")
            log.info("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)

            # event filter for ttbar signal region
            region_conf = copy.deepcopy(process_conf["SR_cuts"])
            rdf_SR = func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=split,
                region_cuts=region_conf,
            )

            log.info(f"Filtering events for the signal region. Target process: {process}")
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_SR.Report().Print()
            log.info(out.getvalue())
            log.info("-" * 50)

            # event filter for ttbar application region
            region_conf = copy.deepcopy(process_conf["AR_cuts"])
            rdf_AR = func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=split,
                region_cuts=region_conf,
            )
            log.info(f"Filtering events for the application region. Target process: {process}")
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_AR.Report().Print()
            log.info(out.getvalue())
            log.info("-" * 50)

            # get binning of the dependent variable
            xbinning = array.array("d", binning)
            nbinsx = len(binning) - 1

            # making the histograms
            h = rdf_SR.Histo1D(
                (
                    process_conf["var_dependence"],
                    f"{sample}",
                    nbinsx,
                    xbinning,
                ),
                process_conf["var_dependence"],
                "weight",
            )
            SR_hists[sample] = h.GetValue()

            h = rdf_AR.Histo1D(
                (
                    process_conf["var_dependence"],
                    f"{sample}",
                    nbinsx,
                    xbinning,
                ),
                process_conf["var_dependence"],
                "weight",
            )
            AR_hists[sample] = h.GetValue()

    # Start of the FF calculation
    FF_hist = func.calculate_ttbar_FF(
        SR=SR_hists,
        AR=AR_hists,
        SRlike=SRlike_hists,
        ARlike=ARlike_hists,
    )
    # performing the fit and calculating the uncertainties
    fit_graphs, corrlib_exp, used_fit = func.fit_function(
        ff_hists=FF_hist.Clone(),
        bin_edges=binning,
        logger=logger,
        fit_option=process_conf.get("fit_option", gd.default_fit_options["ttbar"]),
        limit_kwargs=process_conf.get(
            "limit_kwargs",
            func.Defaults.fit_function_limit_kwargs(binning),
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
        save_data=True
    )

    # doing some control plots
    data = "data"
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

    for _hist, _region, _data, _samples in [
        (SRlike_hists, "SR_like", data, samples),
        (ARlike_hists, "AR_like", data, samples),
        (SRlike_hists, "SR_like", "data_subtracted", ["ttbar_J"]),
        (ARlike_hists, "AR_like", "data_subtracted", ["ttbar_J"]),
    ]:
        for yscale, save_data in zip(["linear", "log"], [True, False]):
            plotting.plot_data_mc_ratio(
                variable="metphi",
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

    keys = [f"{var}#{split[var]}" for var in split_variables]
    if len(keys) == 1:
        corrlib_expression[keys[0]] = corrlib_exp
    elif len(keys) == 2:
        corrlib_expression.setdefault(keys[0], {})[keys[1]] = corrlib_exp
    else:
        raise Exception("Something went wrong with the category splitting.")

    return corrlib_expression


def calculation_ttbar_FFs(
    config: Dict[str, Union[str, Dict, List]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    logger: str,
) -> Dict[str, Union[Dict[str, str], Dict[str, Dict[str, str]]]]:
    """
    This function calculates fake factors for ttbar.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        process: This is relevant for ttbar because for the tt channel two different ttbar fake factors are calculated, one for each hadronic tau
        logger: Name of the logger that should be used

    Return:
        Dictionary where the categories are defined as keys and and the values are the fitted functions (including variations)
        e.g. corrlib_expressions[CATEGORY_1][CATEGORY_2][VARIATION] if dimension of categories is 2
    """
    log = logging.getLogger(logger)

    # init histogram dict for FF data correction
    SRlike_hists = dict()
    ARlike_hists = dict()
    # init histogram dict for QCD SS/OS estimation
    SRlike_hists_qcd = dict()
    ARlike_hists_qcd = dict()

    process_conf = config["target_processes"][process]

    split_variables, split_combinations, split_binnings = func.get_split_combinations(
        categories=process_conf["split_categories"],
        binning=process_conf["var_bins"],
    )

    # calculating global histograms for the data/mc scale factor
    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        log.info(f"Processing {sample} for the ttbar global data/mc scale factor.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for ttbar signal-like region
        region_conf = copy.deepcopy(process_conf["SRlike_cuts"])
        rdf_SRlike = func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
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

        rdf_SRlike_qcd = func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for QCD estimation in the signal-like region. Target process: {process}")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike_qcd.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # event filter for ttbar application-like region
        region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
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

        rdf_ARlike_qcd = func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
            region_cuts=region_conf,
        )

        log.info(f"Filtering events for QCD estimation in the application-like region. Target process: {process}")
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike_qcd.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # make yield histograms for FF data correction
        h = rdf_SRlike.Histo1D(
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        SRlike_hists[sample] = h.GetValue()

        h = rdf_ARlike.Histo1D(
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        ARlike_hists[sample] = h.GetValue()

        # make yield histograms for QCD estimation
        h_qcd = rdf_SRlike_qcd.Histo1D(
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        SRlike_hists_qcd[sample] = h_qcd.GetValue()

        h_qcd = rdf_ARlike_qcd.Histo1D(
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        ARlike_hists_qcd[sample] = h_qcd.GetValue()

    # calculate QCD estimation
    SRlike_hists["QCD"] = func.QCD_SS_estimate(hists=SRlike_hists_qcd)
    ARlike_hists["QCD"] = func.QCD_SS_estimate(hists=ARlike_hists_qcd)

    # calculate ttbar enriched data by subtraction all there backgrould sample
    SRlike_hists["data_subtracted"] = SRlike_hists["data"].Clone()
    ARlike_hists["data_subtracted"] = ARlike_hists["data"].Clone()

    for hist in SRlike_hists:
        if hist not in ["data", "data_subtracted", "ttbar_J"]:
            SRlike_hists["data_subtracted"].Add(SRlike_hists[hist], -1)
    for hist in ARlike_hists:
        if hist not in ["data", "data_subtracted", "ttbar_J"]:
            ARlike_hists["data_subtracted"].Add(ARlike_hists[hist], -1)

    # splitting between different categories for mc-based FF calculation

    args_list = [
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
            SRlike_hists,
            ARlike_hists,
        )
        for split, binning in zip(split_combinations, split_binnings)
    ]

    if len(split_combinations) == 1:
        results = [_split_calculation_ttbar_FFs(args_list[0])]
    else:
        with concurrent.futures.ProcessPoolExecutor(max_workers=len(split_combinations)) as executor:
            results = list(executor.map(_split_calculation_ttbar_FFs, args_list))

    corrlib_expressions = dict()
    for result in results:
        if len(split_variables) == 1:
            corrlib_expressions.update(result)
        elif len(split_variables) == 2:
            key = list(result.keys())[0]
            corrlib_expressions.setdefault(key, {}).update(result[key])

    return corrlib_expressions


def non_closure_correction(
    config: Dict[str, Union[str, Dict, List]],
    corr_config: Dict[str, Union[str, Dict]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    closure_variable: str,
    evaluator: FakeFactorEvaluator,
    corr_evaluators: List[FakeFactorCorrectionEvaluator],
    logger: str,
    **kwargs: Dict[str, Any],
) -> Dict[str, np.ndarray]:
    """
    This function calculates non closure corrections for fake factors for ttbar.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        corr_config: A dictionary with all the relevant information for calculating corrections to the measured fake factors
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        process: This is relevant for ttbar because for the tt channel two different ttbar fake factors are calculated, one for each hadronic tau
        closure_variable: Name of the variable the correction is calculated for
        evaluator: Evaluator with ttbar fake factors
        corr_evaluators: List of evaluators with corrections to ttbar fake factors
        logger: Name of the logger that should be used

    Return:
        Dictionary of arrays with information about the smoothed function values to be stored with correctionlib (nominal and variations)
    """
    log = logging.getLogger(logger)

    # init histogram dict for FF measurement
    SR_hists = dict()
    AR_hists = dict()

    # get process specific config information
    process_conf = copy.deepcopy(config["target_processes"][process])
    correction_conf = corr_config["target_processes"][process]["non_closure"][
        closure_variable
    ]

    for sample_path in sample_paths:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        if sample == "ttbar_J":
            log.info(
                f"Processing {sample} for the non closure correction for {process}."
            )
            log.info("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)

            # event filter for ttbar signal region
            region_conf = copy.deepcopy(process_conf["SR_cuts"])
            rdf_SR = func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=None,
                region_cuts=region_conf,
            )

            log.info(
                f"Filtering events for the signal region. Target process: {process}"
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_SR.Report().Print()
            log.info(out.getvalue())
            log.info("-" * 50)

            # event filter for ttbar application region
            region_conf = copy.deepcopy(process_conf["AR_cuts"])
            rdf_AR = func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=None,
                region_cuts=region_conf,
            )

            log.info(
                f"Filtering events for the application region. Target process: {process}"
            )

            rdf_AR = evaluator.evaluate_fake_factor(rdf=rdf_AR)

            # additionally evaluate the previous corrections
            corr_str = ""
            for corr_evaluator in corr_evaluators:
                rdf_AR = corr_evaluator.evaluate_correction(
                    rdf=rdf_AR,
                )
                corr_str += f" * {process}_ff_corr_{corr_evaluator.variable}"

            rdf_AR = rdf_AR.Define(
                "weight_ff", f"weight * {process}_fake_factor{corr_str}"
            )

            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_AR.Report().Print()
            log.info(out.getvalue())
            log.info("-" * 50)

            # get binning of the dependent variable
            xbinning = array.array("d", correction_conf["var_bins"])
            nbinsx = len(correction_conf["var_bins"]) - 1

            # making the histograms
            h = rdf_SR.Histo1D(
                (
                    correction_conf["var_dependence"],
                    f"{sample}",
                    nbinsx,
                    xbinning,
                ),
                correction_conf["var_dependence"],
                "weight",
            )
            SR_hists[sample] = h.GetValue()

            h = rdf_AR.Histo1D(
                (
                    correction_conf["var_dependence"],
                    f"{sample}",
                    nbinsx,
                    xbinning,
                ),
                correction_conf["var_dependence"],
                "weight_ff",
            )
            AR_hists["ttbar_ff"] = h.GetValue()

    correction_hist = func.calculate_non_closure_correction_ttbar_fromMC(
        SR=SR_hists, AR=AR_hists
    )

    smoothed_graph, correction_dict = func.smooth_function(
        hist=correction_hist.Clone(),
        bin_edges=correction_conf["var_bins"],
        write_corrections=correction_conf["write_corrections"],
    )

    plotting.plot_correction(
        variable=correction_conf["var_dependence"],
        corr_hist=correction_hist,
        corr_graph=smoothed_graph,
        corr_name="non_closure_" + closure_variable,
        era=config["era"],
        channel=config["channel"],
        process=process,
        output_path=output_path,
        logger=logger,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SR_hists["ttbar_J"].Clone()
    plot_hists["data_ff"] = AR_hists["ttbar_ff"].Clone()

    data = "data_subtracted"
    samples = ["data_ff"]

    for yscale in ["linear", "log"]:
        plotting.plot_data_mc(
            variable=correction_conf["var_dependence"],
            hists=plot_hists,
            era=config["era"],
            channel=config["channel"],
            process=process,
            region="non_closure_" + closure_variable,
            data=data,
            samples=samples,
            category={"incl": ""},
            output_path=output_path,
            logger=logger,
            yscale=yscale,
        )

    return correction_dict
