"""
Function for calculating fake factors for the ttbar process
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

    # init histogram dict for FF measurement from MC
    SR_hists = dict()
    AR_hists = dict()
    # init histogram dict for FF data correction
    SRlike_hists = dict()
    ARlike_hists = dict()
    # init histogram dict for QCD SS/OS estimation
    SRlike_hists_qcd = dict()
    ARlike_hists_qcd = dict()
    # init dictionary for the FF functions for correctionlib
    corrlib_expressions = dict()

    # get QCD specific config information
    process_conf = config["target_processes"][process]

    split_variables, split_combinations = func.get_split_combinations(
        categories=process_conf["split_categories"]
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
            raise ValueError(
                "No tau pair sign cut defined in the ttbar config. Is needed for the QCD estimation."
            )

        rdf_SRlike_qcd = func.apply_region_filters(
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

        # event filter for ttbar application-like region
        region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
            region_cuts=region_conf,
        )

        log.info(
            "Filtering events for the application-like region. Target process: ttbar"
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
                "No tau pair sign cut defined in the ttbar config. Is needed for the QCD estimation."
            )

        rdf_ARlike_qcd = func.apply_region_filters(
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
    for split in split_combinations:
        for sample_path in sample_paths:
            # getting the name of the process from the sample path
            sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
            # FFs for ttbar from mc -> only ttbar with true misindentified jets relevant
            if sample in ["ttbar_J"]:
                log.info(
                    f"Processing {sample} for the {', '.join([f'{var} {split[var]}' for var in split_variables])} category."
                )
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
                    category_cuts=split,
                    region_cuts=region_conf,
                )
                log.info(
                    f"Filtering events for the application region. Target process: {process}"
                )
                # redirecting C++ stdout for Report() to python stdout
                out = StringIO()
                with pipes(stdout=out, stderr=STDOUT):
                    rdf_AR.Report().Print()
                log.info(out.getvalue())
                log.info("-" * 50)

                # get binning of the dependent variable
                xbinning = array.array("d", process_conf["var_bins"])
                nbinsx = len(process_conf["var_bins"]) - 1

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
            draw_option=used_fit,
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

        # doing some control plots
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

        plotting.plot_data_mc_ratio(
            variable="metphi",
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
        plotting.plot_data_mc_ratio(
            variable="metphi",
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

        data = "data_subtracted"
        samples = ["ttbar_J"]

        plotting.plot_data_mc_ratio(
            variable="metphi",
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
        plotting.plot_data_mc_ratio(
            variable="metphi",
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
        corr_evaluator: Evaluator with corrections to ttbar fake factors
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
            log.info(f"Processing {sample} for the non closure correction for {process}.")
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
                category_cuts=None,
                region_cuts=region_conf,
            )

            log.info(
                f"Filtering events for the application region. Target process: {process}"
            )

            # evaluate the measured fake factors for the specific processes
            if config["channel"] != "tt" or process == "ttbar_subleading":
                rdf_AR = evaluator.evaluate_subleading_lep_pt_njets(rdf=rdf_AR)
            else:
                rdf_AR = evaluator.evaluate_leading_lep_pt_njets(rdf=rdf_AR)
            if corr_evaluator == None:
                rdf_AR = rdf_AR.Define("weight_ff", f"weight * {process}_fake_factor")
            else:
                if config["channel"] != "tt" or process == "QCD_subleading":
                    rdf_AR = corr_evaluator.evaluate_leading_lep_pt(rdf=rdf_AR)
                else:
                    rdf_AR = corr_evaluator.evaluate_subleading_lep_pt(rdf=rdf_AR)
                rdf_AR = rdf_AR.Define(
                    "weight_ff", f"weight * {process}_fake_factor * {process}_ff_corr"
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
    )

    return correction_dict


def non_closure_correction_tt(
    config: Dict[str, Union[str, Dict, List]],
    corr_config: Dict[str, Union[str, Dict]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    closure_variable: str,
    evaluator: FakeFactorEvaluator,
    corr_evaluator: FakeFactorCorrectionEvaluator,
    logger: str,
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
        corr_evaluator: Evaluator with corrections to ttbar fake factors
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
        log.info(f"Processing {sample} for the non closure correction for {process}.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for ttbar signal region
        region_conf = copy.deepcopy(process_conf["SRlike_cuts"])
        rdf_SR = func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=None,
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
        region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
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

        # evaluate the measured fake factors for the specific processes
        if sample == "data":
            if config["channel"] != "tt" or process == "ttbar_subleading":
                rdf_AR = evaluator.evaluate_subleading_lep_pt_njets(rdf=rdf_AR)
            else:
                rdf_AR = evaluator.evaluate_leading_lep_pt_njets(rdf=rdf_AR)
            if corr_evaluator == None:
                rdf_AR = rdf_AR.Define("weight_ff", f"weight * {process}_fake_factor")
            else:
                if config["channel"] != "tt" or process == "QCD_subleading":
                    rdf_AR = corr_evaluator.evaluate_leading_lep_pt(rdf=rdf_AR)
                else:
                    rdf_AR = corr_evaluator.evaluate_subleading_lep_pt(rdf=rdf_AR)
                rdf_AR = rdf_AR.Define(
                    "weight_ff", f"weight * {process}_fake_factor * {process}_ff_corr"
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
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        AR_hists[sample] = h.GetValue()

        if sample == "data":
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
            AR_hists["data_ff"] = h.GetValue()
        
    SR_hists["data_subtracted"] = SR_hists["data"].Clone()
    AR_hists["data_subtracted"] = AR_hists["data"].Clone()

    for hist in SR_hists:
        if hist not in ["data", "data_subtracted", "ttbar_J"]:
            SR_hists["data_subtracted"].Add(SR_hists[hist], -1)
    for hist in AR_hists:
        if hist not in ["data", "data_subtracted", "data_ff", "ttbar_J"]:
            AR_hists["data_subtracted"].Add(AR_hists[hist], -1)

    correction_hist, process_fraction = func.calculate_non_closure_correction(
        SRlike=SR_hists, ARlike=AR_hists
    )

    smoothed_graph, correction_dict = func.smooth_function(
        hist=correction_hist.Clone(),
        bin_edges=correction_conf["var_bins"],
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

    plot_hists["data_subtracted"] = SR_hists["data_subtracted"].Clone()
    plot_hists["data_ff"] = AR_hists["data_ff"].Clone()
    plot_hists["data_ff"].Scale(process_fraction)

    data = "data_subtracted"
    samples = ["data_ff"]

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
        hists=SR_hists,
        era=config["era"],
        channel=config["channel"],
        process=process,
        region="non_closure_" + closure_variable + "_SRlike_hist",
        data=data,
        samples=samples,
        category={"incl": ""},
        output_path=output_path,
        logger=logger,
    )

    return correction_dict
