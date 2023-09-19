"""
Function for calculating fake factors for the ttbar process
"""

import array
import copy
import ROOT
from io import StringIO
from wurlitzer import pipes, STDOUT
import logging
from typing import Union, Dict, List

import helper.ff_functions as func
import helper.plotting as plotting


def calculation_ttbar_FFs(config: Dict[str, Union[str, Dict, List]], sample_paths: List[str], output_path: str) -> Dict[str, Union[Dict[str, str], Dict[str, Dict[str, str]]]]:
    '''
    This function calculates fake factors for ttbar.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
    
    Return:
        Dictionary where the categories are defined as keys and and the values are the fitted functions (including variations)
        e.g. corrlib_expressions[CATEGORY_1][CATEGORY_2][VARIATION] if dimension of categories is 2
    '''
    log = logging.getLogger("ff_calculation")

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
    process_conf = config["target_processes"]["ttbar"]

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
        rdf_SRlike = func.apply_region_filters(rdf=rdf, channel=config["channel"], sample=sample, category_cuts=None, region_cuts=region_conf)

        log.info(
            "Filtering events for the signal-like region. Target process: ttbar"
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # QCD estimation from same sign in signal-like region
        if "tau_pair_sign" in region_conf:
            region_conf["tau_pair_sign"] = "(q_1*q_2) > 0" # same sign
        else:
            raise ValueError("No tau pair sign cut defined in the ttbar config. Is needed for the QCD estimation.")
        
        rdf_SRlike_qcd = func.apply_region_filters(rdf=rdf, channel=config["channel"], sample=sample, category_cuts=None, region_cuts=region_conf)

        log.info(
            "Filtering events for QCD estimation in the signal-like region. Target process: ttbar"
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike_qcd.Report().Print()
        log.info(out.getvalue())
        log.info("-" * 50)

        # event filter for ttbar application-like region
        region_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = func.apply_region_filters(rdf=rdf, channel=config["channel"], sample=sample, category_cuts=None, region_cuts=region_conf)

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
            region_conf["tau_pair_sign"] = "(q_1*q_2) > 0" # same sign
        else:
            raise ValueError("No tau pair sign cut defined in the ttbar config. Is needed for the QCD estimation.")
        
        rdf_ARlike_qcd = func.apply_region_filters(rdf=rdf, channel=config["channel"], sample=sample, category_cuts=None, region_cuts=region_conf)

        log.info(
            "Filtering events for QCD estimation in the application-like region. Target process: ttbar"
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
                    f"Processing {sample} for the {', '.join([f'{var} {split[var]}'for var in split_variables])} category."
                )
                log.info("-" * 50)

                rdf = ROOT.RDataFrame(config["tree"], sample_path)

                # event filter for ttbar signal region
                region_conf = copy.deepcopy(process_conf["SR_cuts"])
                rdf_SR = func.apply_region_filters(rdf=rdf, channel=config["channel"], sample=sample, category_cuts=split, region_cuts=region_conf)
                
                log.info(
                    "Filtering events for the signal region. Target process: ttbar"
                )
                # redirecting C++ stdout for Report() to python stdout
                out = StringIO()
                with pipes(stdout=out, stderr=STDOUT):
                    rdf_SR.Report().Print()
                log.info(out.getvalue())
                log.info("-" * 50)

                # event filter for ttbar application region
                region_conf = copy.deepcopy(process_conf["AR_cuts"])
                rdf_AR = func.apply_region_filters(rdf=rdf, channel=config["channel"], sample=sample, category_cuts=split, region_cuts=region_conf)
                log.info(
                    "Filtering events for the application region. Target process: ttbar"
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
            SR=SR_hists, AR=AR_hists, SRlike=SRlike_hists, ARlike=ARlike_hists
        )
        # performing the fit and calculating the uncertainties
        fit_graphs, corrlib_exp = func.fit_function(
            ff_hists=FF_hist.Clone(), 
            bin_edges=process_conf["var_bins"],
        )

        plotting.plot_FFs(
            variable=process_conf["var_dependence"], 
            ff_ratio=FF_hist, 
            uncertainties=fit_graphs, 
            era=config["era"],
            channel=config["channel"],
            process="ttbar", 
            category=split, 
            output_path=output_path,
        )

        if len(split) == 1:
            corrlib_expressions[f"{split_variables[0]}#{split[split_variables[0]]}"] = corrlib_exp
        elif len(split) == 2:
            if f"{split_variables[0]}#{split[split_variables[0]]}" not in corrlib_expressions:
                corrlib_expressions[f"{split_variables[0]}#{split[split_variables[0]]}"] = dict()
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
            process="ttbar", 
            region="SR_like",
            data=data,
            samples=samples,
            category=split,
            output_path=output_path,

        )
        plotting.plot_data_mc_ratio(
            variable="metphi",
            hists=ARlike_hists,
            era=config["era"],
            channel=config["channel"],
            process="ttbar", 
            region="AR_like",
            data=data,
            samples=samples,
            category=split,
            output_path=output_path,
        )

        data = "data_subtracted"
        samples = ["ttbar_J"]

        plotting.plot_data_mc_ratio(
            variable="metphi",
            hists=SRlike_hists,
            era=config["era"],
            channel=config["channel"],
            process="ttbar", 
            region="SR_like",
            data=data,
            samples=samples,
            category=split,
            output_path=output_path,
        )
        plotting.plot_data_mc_ratio(
            variable="metphi",
            hists=ARlike_hists,
            era=config["era"],
            channel=config["channel"],
            process="ttbar", 
            region="AR_like",
            data=data,
            samples=samples,
            category=split,
            output_path=output_path,
        )
        log.info("-" * 50)

    return corrlib_expressions


def non_closure_correction(
    config,
    corr_config,
    closure_variable,
    sample_path_list,
    save_path,
    evaluator,
    corr_evaluator,
):
    # init histogram dict for FF measurement
    SR_hists = dict()
    AR_hists = dict()

    # get process specific config information
    process_conf = copy.deepcopy(config["target_process"]["ttbar"])
    correction_conf = corr_config["target_process"]["ttbar"]["non_closure"][
        closure_variable
    ]
    boosted = True if "boosted" in process_conf["var_dependence"] else False

    for sample_path in sample_path_list:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        if sample == "ttbar_J":
            print(
                "Processing {sample} for the non closure correction for {process}.".format(
                    sample=sample,
                    process="ttbar",
                )
            )
            print("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)

            # event filter for ttbar signal region
            region_cut_conf = copy.deepcopy(process_conf["SR_cuts"])
            rdf_SR = region_filter(rdf, config["channel"], region_cut_conf, sample)
            print(
                "Filtering events for the signal region. Target process: {}\n".format(
                    "ttbar"
                )
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_SR.Report().Print()
            print(out.getvalue())
            print("-" * 50)

            # event filter for ttbar application region
            region_cut_conf = copy.deepcopy(process_conf["AR_cuts"])
            rdf_AR = region_filter(rdf, config["channel"], region_cut_conf, sample)
            print(
                "Filtering events for the application region. Target process: {}\n".format(
                    "ttbar"
                )
            )

            # evaluate the measured fake factors for the specific processes
            if not boosted:
                rdf_AR = evaluator.evaluate_subleading_lep_pt_njets(rdf_AR)
                if corr_evaluator == None:
                    rdf_AR = rdf_AR.Define("weight_ff", "weight * ttbar_fake_factor")
                else:
                    rdf_AR = corr_evaluator.evaluate_leading_lep_pt(rdf_AR)
                    rdf_AR = rdf_AR.Define(
                        "weight_ff", "weight * ttbar_fake_factor * ttbar_ff_corr"
                    )
            elif boosted:
                rdf_AR = evaluator.evaluate_subleading_boosted_lep_pt_njets(rdf_AR)
                if corr_evaluator == None:
                    rdf_AR = rdf_AR.Define("weight_ff", "weight * ttbar_fake_factor")
                else:
                    rdf_AR = corr_evaluator.evaluate_leading_boosted_lep_pt(rdf_AR)
                    rdf_AR = rdf_AR.Define(
                        "weight_ff", "weight * ttbar_fake_factor * ttbar_ff_corr"
                    )

            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_AR.Report().Print()
            print(out.getvalue())
            print("-" * 50)

            # get binning of the dependent variable
            xbinning = array.array("d", correction_conf["var_bins"])
            nbinsx = len(correction_conf["var_bins"]) - 1

            # making the histograms
            h = rdf_SR.Histo1D(
                (
                    correction_conf["var_dependence"],
                    "{}".format(sample),
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
                    "{}".format(sample),
                    nbinsx,
                    xbinning,
                ),
                correction_conf["var_dependence"],
                "weight_ff",
            )
            AR_hists["ttbar_ff"] = h.GetValue()

    corr_hist = func.calculate_non_closure_correction_ttbar(SR_hists, AR_hists)

    smooth_graph, corr_def = func.smooth_function(
        corr_hist.Clone(), correction_conf["var_bins"]
    )

    plotting.plot_correction(
        corr_hist,
        smooth_graph,
        correction_conf["var_dependence"],
        "ttbar",
        "non_closure_" + closure_variable,
        config,
        save_path,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SR_hists["ttbar_J"].Clone()
    plot_hists["data_ff"] = AR_hists["ttbar_ff"].Clone()

    data = "data_subtracted"
    samples = ["data_ff"]
    plotting.plot_data_mc(
        plot_hists,
        config,
        correction_conf["var_dependence"],
        "ttbar",
        "non_closure_" + closure_variable,
        data,
        samples,
        {"incl": ""},
        save_path,
    )

    return corr_def
