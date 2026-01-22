"""
Function for calculating fake factors for the W-jets process
"""

import array
import copy
import logging
from copy import deepcopy
from typing import Any, Dict, Tuple, Union

import numpy as np
import ROOT

import helper.ff_functions as ff_func
import helper.plotting as plotting
from helper.functions import RuntimeVariables


def calculation_Wjets_FFs(
    args: Tuple[Any, ...],
) -> Dict[str, Union[Dict[str, str], Dict[str, Dict[str, str]]]]:
    """
    This function calculates fake factors for the Wjets process for a specific category (split).

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
        lock,  # lock: str, (needed for cache snapshots)
        *_,  # SRlike_hists, ARlike_hists  only used in ttbar calculation
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
        log.info(f"Processing {sample} for the {', '.join([f'{var} {splitting.split[var]}' for var in splitting.variables])} category.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for Wjets signal-like region
        log.info(f"Filtering events for the signal-like region. Target process: {process}")
        region_conf = copy.deepcopy(process_conf["SRlike_cuts"])
        rdf_SRlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
            logger=logger,
            lock=lock,
        )

        # QCD estimation from same sign in signal-like region
        if "tau_pair_sign" in region_conf:
            region_conf["tau_pair_sign"] = "(q_1*q_2) > 0"  # same sign
        else:
            raise ValueError(f"No tau pair sign cut defined in the {process} config. Is needed for the QCD estimation.")

        log.info(f"Filtering events for QCD estimation in the signal-like region. Target process: {process}")
        rdf_SRlike_qcd = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
            logger=logger,
            lock=lock,
        )

        # event filter for Wjets application-like region
        log.info(f"Filtering events for the application-like region. Target process: {process}")
        region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
            logger=logger,
            lock=lock,
        )

        # QCD estimation from same sign in application-like region
        if "tau_pair_sign" in region_conf:
            region_conf["tau_pair_sign"] = "(q_1*q_2) > 0"  # same sign
        else:
            raise ValueError(
                f"No tau pair sign cut defined in the {process} config. Is needed for the QCD estimation."
            )

        log.info(f"Filtering events for QCD estimation in the application-like region. Target process: {process}")
        rdf_ARlike_qcd = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
            logger=logger,
            lock=lock,
        )

        # get binning of the dependent variable from the computed binning
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

        # making the histograms for QCD estimation
        h_qcd = RuntimeVariables.RDataFrameWrapper(rdf_SRlike_qcd).Histo1D(
            (process_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            process_conf["var_dependence"],
            "weight",
        )
        SRlike_hists_qcd[sample] = h_qcd.GetValue()

        h_qcd = RuntimeVariables.RDataFrameWrapper(rdf_ARlike_qcd).Histo1D(
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
            SRlike_hists["data_subtracted"].Add(SRlike_hists[hist].Clone(), -1)
            SRlike_hists["data_subtracted_up"].Add(SRlike_hists[hist].Clone().AddError(1), -1)
            SRlike_hists["data_subtracted_down"].Add(SRlike_hists[hist].Clone().AddError(-1), -1)
    for hist in ARlike_hists:
        if hist not in [
            "data",
            "data_subtracted",
            "data_subtracted_up",
            "data_subtracted_down",
            "Wjets",
        ]:
            ARlike_hists["data_subtracted"].Add(ARlike_hists[hist].Clone(), -1)
            ARlike_hists["data_subtracted_up"].Add(ARlike_hists[hist].Clone().AddError(1), -1)
            ARlike_hists["data_subtracted_down"].Add(ARlike_hists[hist].Clone().AddError(-1), -1)

    # Start of the FF calculation
    FF_hist, FF_hist_up, FF_hist_down = ff_func.calculate_Wjets_FF(
        SRlike=SRlike_hists, ARlike=ARlike_hists
    )
    # performing the fit and calculating the uncertainties
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
        category=splitting.split or {"incl": ""},
        output_path=output_path,
        logger=logger,
        draw_option=used_fit,
        save_data=True,
    )

    # producing some control plots
    for _hist, _region, _data, _samples in [
        (SRlike_hists, "SR_like", "data", ff_func.controlplot_samples(config["use_embedding"])),
        (ARlike_hists, "AR_like", "data", ff_func.controlplot_samples(config["use_embedding"])),
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
                category=splitting.split or {"incl": ""},
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
    This function calculates the non closure correction for the Wjet process for a specific category.

    Intended to be used in a multiprocessing environment.

    Args:
        args: Tuple containing all the necessary information for the calculation of the non-closure correction
            splitting: SplitQuantitiesContainer, contains the splitting information
            config: Dictionary with all the relevant information for the fake factor calculation
            correction_conf: Dictionary with all the relevant information for the non-closure correction
            process: Name of the process
            closure_variable: Name of the variable dependence of the closure correction
            sample_paths: List of file paths where the samples are stored
            output_path: Path where the generated plots should be stored
            logger: Name of the logger that should be used
            evaluator: FakeFactorEvaluator object
            corr_evaluators: List of FakeFactorCorrectionEvaluator objects
            for_DRtoSR: If True the correction is calculated for the DR to SR correction
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
        lock,  # lock: str, (needed for cache snapshots)
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
        if splitting.split is not None:
            log.info(f"Processing {sample} for the non closure correction for {process} for {', '.join([f'{var} {splitting.split[var]}' for var in splitting.variables])}.")
        else:
            log.info(f"Processing {sample} for the non closure correction for {process}.")
        log.info("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for Wjets signal-like region
        log.info(f"Filtering events for the signal-like region. Target process: {process}")
        region_conf = copy.deepcopy(config["target_processes"][process]["SRlike_cuts"])
        rdf_SRlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
            logger=logger,
            lock=lock,
        )

        # QCD estimation from same sign in signal-like region
        if "tau_pair_sign" in region_conf:
            region_conf["tau_pair_sign"] = "(q_1*q_2) > 0"  # same sign
        else:
            raise ValueError(f"No tau pair sign cut defined in the {process} config. Is needed for the QCD estimation.")

        log.info(f"Filtering events for QCD estimation in the signal-like region. Target process: {process}")
        rdf_SRlike_qcd = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
            logger=logger,
            lock=lock,
        )

        # event filter for Wjets application-like region
        log.info(f"Filtering events for the application-like region. Target process: {process}")
        region_conf = copy.deepcopy(config["target_processes"][process]["ARlike_cuts"])
        rdf_ARlike = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
            logger=logger,
            lock=lock,
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
                corr_str += f" * {corr_evaluator.corr_str}"

            rdf_ARlike = rdf_ARlike.Define(
                "weight_ff", f"weight * {process}_fake_factor{corr_str}"
            )

        # QCD estimation from same sign in application-like region
        if "tau_pair_sign" in region_conf:
            region_conf["tau_pair_sign"] = "(q_1*q_2) > 0"  # same sign
        else:
            raise ValueError(f"No tau pair sign cut defined in the {process} config. Is needed for the QCD estimation.")

        log.info(f"Filtering events for QCD estimation in the application-like region. Target process: {process}")
        rdf_ARlike_qcd = ff_func.apply_region_filters(
            rdf=rdf,
            channel=config["channel"],
            sample=sample,
            category_cuts=splitting.split,
            region_cuts=region_conf,
            logger=logger,
            lock=lock,
        )

        # get binning of the dependent variable
        xbinning = array.array("d", splitting.var_bins)
        nbinsx = len(splitting.var_bins) - 1

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

        # making the histograms for QCD estimation
        h_qcd = RuntimeVariables.RDataFrameWrapper(rdf_SRlike_qcd).Histo1D(
            (correction_conf["var_dependence"], f"{sample}", nbinsx, xbinning),
            correction_conf["var_dependence"],
            "weight",
        )
        SRlike_hists_qcd[sample] = h_qcd.GetValue()

        h_qcd = RuntimeVariables.RDataFrameWrapper(rdf_ARlike_qcd).Histo1D(
            ("#phi(#slash{E}_{T})", f"{sample}", 1, -3.5, 3.5), "metphi", "weight"
        )
        ARlike_hists_qcd[sample] = h_qcd.GetValue()

    # calculate QCD estimation
    SRlike_hists["QCD"] = ff_func.QCD_SS_estimate(hists=SRlike_hists_qcd)
    ARlike_hists["QCD"] = ff_func.QCD_SS_estimate(hists=ARlike_hists_qcd)

    SRlike_hists["data_subtracted"] = SRlike_hists["data"].Clone()
    ARlike_hists["data_subtracted"] = ARlike_hists["data"].Clone()

    _pairs = [("data_subtracted", "data"), ("data", "data")]

    SRlike_hists_sub_up = {k1: SRlike_hists[k2].Clone() for k1, k2 in _pairs}
    SRlike_hists_sub_down = deepcopy(SRlike_hists_sub_up)

    ARlike_hists_sub_up = {k1: ARlike_hists[k2].Clone() for k1, k2 in _pairs + [("data_ff", "data_ff")]}
    ARlike_hists_sub_down = deepcopy(ARlike_hists_sub_up)

    for hist in SRlike_hists:
        if hist not in ["data", "data_subtracted", "Wjets"]:
            SRlike_hists["data_subtracted"].Add(SRlike_hists[hist].Clone(), -1)
            SRlike_hists_sub_up["data_subtracted"].Add(SRlike_hists[hist].Clone().AddError(1), -1)
            SRlike_hists_sub_down["data_subtracted"].Add(SRlike_hists[hist].Clone().AddError(-1), -1)
    for hist in ARlike_hists:
        if hist not in ["data", "data_subtracted", "data_ff", "Wjets"]:
            ARlike_hists["data_subtracted"].Add(ARlike_hists[hist].Clone(), -1)
            ARlike_hists_sub_up["data_subtracted"].Add(ARlike_hists[hist].Clone().AddError(1), -1)
            ARlike_hists_sub_down["data_subtracted"].Add(ARlike_hists[hist].Clone().AddError(-1), -1)

    correction_hist, process_fraction = ff_func.calculate_non_closure_correction(
        SRlike=SRlike_hists,
        ARlike=ARlike_hists,
    )

    nominal_draw_obj, smoothed_graph, correction_dict = ff_func.smooth_function(
        hist=correction_hist.Clone(),
        bin_edges=splitting.var_bins,
        correction_option=splitting.correction_option,
        bandwidth=splitting.bandwidth,
        mc_shifted_hist={
            "MCShiftUp": ff_func.calculate_non_closure_correction(SRlike_hists_sub_up, ARlike_hists_sub_up)[0].Clone(),
            "MCShiftDown": ff_func.calculate_non_closure_correction(SRlike_hists_sub_down, ARlike_hists_sub_down)[0].Clone(),
        },
    )

    add_str = "_DRtoSR" if for_DRtoSR else ""
    plotting.plot_correction(
        variable=correction_conf["var_dependence"],
        corr_hist=nominal_draw_obj,
        corr_graph=correction_dict,
        corr_name=f"non_closure_{closure_variable}{add_str}",
        era=config["era"],
        channel=config["channel"],
        process=process,
        output_path=output_path,
        logger=logger,
        category=splitting.split or {"incl": ""},
        save_data=True,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SRlike_hists["data_subtracted"].Clone()
    plot_hists["data_ff"] = ARlike_hists["data_ff"].Clone()
    plot_hists["data_ff"].Scale(process_fraction)

    for yscale, save_data in zip(["linear", "log"], [True, False]):
        plotting.plot_data_mc_ratio(
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
            save_data=save_data,
        )

    for yscale, save_data in zip(["linear", "log"], [True, False]):
        plotting.plot_data_mc_ratio(
            variable=correction_conf["var_dependence"],
            hists=SRlike_hists,
            era=config["era"],
            channel=config["channel"],
            process=process,
            region=f"non_closure_{closure_variable}{add_str}_SRlike",
            data="data",
            samples=ff_func.controlplot_samples(config["use_embedding"]),
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
    This function calculates DR to SR correction for fake factors for Wjets.

    Intended to be used in a multiprocessing environment.

    Args:
        args: Tuple containing all the necessary information for the calculation of the DR to SR correction
            splitting: SplitQuantitiesContainer, contains the splitting information
            config: Dictionary with all the relevant information for the fake factor calculation
            correction_conf: Dictionary with all the relevant information for the correction calculation
            process: Name of the process
            split_variables: List of variables that are used for the category splitting
            sample_paths: List of file paths where the samples are stored
            output_path: Path where the generated plots should be stored
            logger: Name of the logger that should be used
            evaluator: Evaluator with Wjets fake factors
            corr_evaluators: List of evaluators with corrections to Wjets fake factors
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
        if sample == "Wjets":
            if splitting.split is not None:
                log.info(f"Processing {sample} for the DR to SR correction for {process} for {', '.join([f'{var} {splitting.split[var]}' for var in splitting.variables])}.")
            else:
                log.info(f"Processing {sample} for the DR to SR correction for {process}.")
            log.info("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)

            # event filter for Wjets signal-like region
            log.info(f"Filtering events for the signal-like region. Target process: {process}")
            region_conf = copy.deepcopy(config["target_processes"][process]["SRlike_cuts"])
            rdf_SRlike = ff_func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=splitting.split,
                region_cuts=region_conf,
                logger=logger,
                lock=lock,
            )

            # event filter for Wjets application-like region
            log.info(f"Filtering events for the application-like region. Target process: {process}")
            region_conf = copy.deepcopy(config["target_processes"][process]["ARlike_cuts"])
            rdf_ARlike = ff_func.apply_region_filters(
                rdf=rdf,
                channel=config["channel"],
                sample=sample,
                category_cuts=splitting.split,
                region_cuts=region_conf,
                logger=logger,
                lock=lock,
            )

            rdf_ARlike = evaluator.evaluate_fake_factor(rdf=rdf_ARlike)
            # additionally evaluate the previous corrections
            corr_str = ""
            for corr_evaluator in corr_evaluators:
                rdf_ARlike = corr_evaluator.evaluate_correction(
                    rdf=rdf_ARlike,
                )
                corr_str += f" * {corr_evaluator.corr_str}"

            rdf_ARlike = rdf_ARlike.Define(
                "weight_ff", f"weight * {process}_fake_factor{corr_str}"
            )

            # get binning of the dependent variable
            xbinning, nbinsx = array.array("d", splitting.var_bins), len(splitting.var_bins) - 1

            # making the histograms
            h = RuntimeVariables.RDataFrameWrapper(rdf_SRlike).Histo1D(
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
            ARlike_hists["Wjets_ff"] = h.GetValue()

    correction_hist = ff_func.calculate_non_closure_correction_Wjets_fromMC(
        SRlike=SRlike_hists, ARlike=ARlike_hists
    )

    nominal_draw_obj, smoothed_graph, correction_dict = ff_func.smooth_function(
        hist=correction_hist.Clone(),
        bin_edges=splitting.var_bins,
        correction_option=splitting.correction_option,
        bandwidth=splitting.bandwidth,
    )

    plotting.plot_correction(
        variable=correction_conf["var_dependence"],
        corr_hist=nominal_draw_obj,
        corr_graph=correction_dict,
        corr_name="DR_SR",
        era=config["era"],
        channel=config["channel"],
        process="Wjets",
        output_path=output_path,
        category=splitting.split or {"incl": ""},
        logger=logger,
        save_data=True,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SRlike_hists["Wjets"].Clone()
    plot_hists["data_ff"] = ARlike_hists["Wjets_ff"].Clone()

    for yscale, save_data in zip(["linear", "log"], [True, False]):
        plotting.plot_data_mc_ratio(
            variable=correction_conf["var_dependence"],
            hists=plot_hists,
            era=config["era"],
            channel=config["channel"],
            process="Wjets",
            region="DR_SR",
            data="data_subtracted",
            samples=["data_ff"],
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
