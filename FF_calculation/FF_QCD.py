"""
Function for calculating fake factors for the QCD process
"""

import sys
import array
import copy
import ROOT
from io import StringIO
from wurlitzer import pipes, STDOUT

import helper.functions as func
import helper.plotting as plotting

from FF_calculation.FF_region_filters import region_filter


def calculation_QCD_FFs(config, sample_path_list, save_path):
    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()
    # init dictionary for the FF functions for correctionlib
    cs_expressions = dict()

    # get QCD specific config information
    process_conf = config["target_process"]["QCD"]

    split_vars, split_combinations = func.get_split_combinations(
        process_conf["split_categories"]
    )

    # splitting between different categories
    for split in split_combinations:
        for sample_path in sample_path_list:
            # getting the name of the process from the sample path
            sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
            print(
                "Processing {sample} for the {cat} category.".format(
                    sample=sample,
                    cat=", ".join(
                        ["{} {}".format(var, split[var]) for var in split_vars]
                    ),
                )
            )
            print("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)

            # event filter for QCD signal-like region
            region_cut_conf = {**split, **process_conf["SRlike_cuts"]}
            rdf_SRlike = region_filter(rdf, config["channel"], region_cut_conf, sample)
            print(
                "Filtering events for the signal-like region. Target process: {}\n".format(
                    "QCD"
                )
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_SRlike.Report().Print()
            print(out.getvalue())
            print("-" * 50)

            # event filter for QCD application-like region
            region_cut_conf = {**split, **process_conf["ARlike_cuts"]}
            rdf_ARlike = region_filter(rdf, config["channel"], region_cut_conf, sample)
            print(
                "Filtering events for the application-like region. Target process: {}\n".format(
                    "QCD"
                )
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_ARlike.Report().Print()
            print(out.getvalue())
            print("-" * 50)

            # get binning of the dependent variable
            xbinning = array.array("d", process_conf["var_bins"])
            nbinsx = len(process_conf["var_bins"]) - 1

            # making the histograms
            h = rdf_SRlike.Histo1D(
                (process_conf["var_dependence"], "{}".format(sample), nbinsx, xbinning),
                process_conf["var_dependence"],
                "weight",
            )
            SRlike_hists[sample] = h.GetValue()

            h = rdf_ARlike.Histo1D(
                (process_conf["var_dependence"], "{}".format(sample), nbinsx, xbinning),
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
                SRlike_hists["data_subtracted_up"].Add(SRlike_hists[hist], -0.93)
                SRlike_hists["data_subtracted_down"].Add(SRlike_hists[hist], -1.07)
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
            SRlike_hists, ARlike_hists
        )
        # performing the fit and calculating the uncertainties
        fit_graphs, cs_exp = func.fit_function(
            [FF_hist.Clone(), FF_hist_up, FF_hist_down], process_conf["var_bins"]
        )

        plotting.plot_FFs(FF_hist, fit_graphs, "QCD", config, split, save_path)

        if len(split) == 1:
            cs_expressions["{}#{}".format(split_vars[0], split[split_vars[0]])] = cs_exp
        elif len(split) == 2:
            if (
                "{}#{}".format(split_vars[0], split[split_vars[0]])
                not in cs_expressions
            ):
                cs_expressions[
                    "{}#{}".format(split_vars[0], split[split_vars[0]])
                ] = dict()
            cs_expressions["{}#{}".format(split_vars[0], split[split_vars[0]])][
                "{}#{}".format(split_vars[1], split[split_vars[1]])
            ] = cs_exp
        else:
            sys.exit("Category splitting is only defined up to 2 dimensions.")

        # doing some control plots
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
            ]

        plotting.plot_data_mc(
            SRlike_hists,
            config,
            process_conf["var_dependence"],
            "QCD",
            "SR_like",
            data,
            samples,
            split,
            save_path,
        )
        plotting.plot_data_mc(
            ARlike_hists,
            config,
            process_conf["var_dependence"],
            "QCD",
            "AR_like",
            data,
            samples,
            split,
            save_path,
        )
        print("-" * 50)

    return cs_expressions


def non_closure_correction(
    config,
    corr_config,
    closure_variable,
    sample_path_list,
    save_path,
    evaluator,
    for_DRtoSR=False,
):
    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    # get process specific config information
    process_conf = copy.deepcopy(config["target_process"]["QCD"])
    if for_DRtoSR:
        correction_conf = corr_config["target_process"]["QCD"]["DR_SR"]["non_closure"][
            closure_variable
        ]
    else:
        correction_conf = corr_config["target_process"]["QCD"]["non_closure"][
            closure_variable
        ]

    for sample_path in sample_path_list:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        print(
            "Processing {sample} for the non closure correction for {process}.".format(
                sample=sample,
                process="QCD",
            )
        )
        print("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for QCD signal-like region
        region_cut_conf = copy.deepcopy(process_conf["SRlike_cuts"])
        rdf_SRlike = region_filter(rdf, config["channel"], region_cut_conf, sample)
        print(
            "Filtering events for the signal-like region. Target process: {}\n".format(
                "QCD"
            )
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike.Report().Print()
        print(out.getvalue())
        print("-" * 50)

        # event filter for QCD application-like region
        region_cut_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = region_filter(rdf, config["channel"], region_cut_conf, sample)
        print(
            "Filtering events for the application-like region. Target process: {}\n".format(
                "QCD"
            )
        )

        # evaluate the measured fake factors for the specific processes
        if sample == "data":
            rdf_ARlike = evaluator.evaluate_tau_pt_njets(rdf_ARlike)
            rdf_ARlike = rdf_ARlike.Define("weight_ff", "weight * QCD_fake_factor")

        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike.Report().Print()
        print(out.getvalue())
        print("-" * 50)

        # get binning of the dependent variable
        xbinning = array.array("d", correction_conf["var_bins"])
        nbinsx = len(correction_conf["var_bins"]) - 1

        # making the histograms
        h = rdf_SRlike.Histo1D(
            (correction_conf["var_dependence"], "{}".format(sample), nbinsx, xbinning),
            correction_conf["var_dependence"],
            "weight",
        )
        SRlike_hists[sample] = h.GetValue()

        h = rdf_ARlike.Histo1D(
            ("#phi(#tau_{h})", "{}".format(sample), 1, -3.5, 3.5), "phi_2", "weight"
        )
        ARlike_hists[sample] = h.GetValue()

        if sample == "data":
            h = rdf_ARlike.Histo1D(
                (
                    correction_conf["var_dependence"],
                    "{}_ff".format(sample),
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

    corr_hist, proc_frac = func.calculate_non_closure_correction(
        SRlike_hists, ARlike_hists
    )

    smooth_graph, corr_def = func.smooth_function(
        corr_hist.Clone(), correction_conf["var_bins"]
    )

    if for_DRtoSR:
        add_str = "_for_DRtoSR"
    else:
        add_str = ""

    plotting.plot_correction(
        corr_hist,
        smooth_graph,
        correction_conf["var_dependence"],
        "QCD",
        "non_closure_" + closure_variable + add_str,
        config,
        save_path,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SRlike_hists["data_subtracted"].Clone()
    plot_hists["data_ff"] = ARlike_hists["data_ff"].Clone()
    plot_hists["data_ff"].Scale(proc_frac)

    data = "data_subtracted"
    samples = ["data_ff"]

    plotting.plot_data_mc(
        plot_hists,
        config,
        correction_conf["var_dependence"],
        "QCD",
        "non_closure_" + closure_variable + add_str,
        data,
        samples,
        {"incl": ""},
        save_path,
    )

    return corr_def


def DR_SR_correction(
    config, corr_config, sample_path_list, save_path, evaluator, corr_evaluator
):
    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    # get process specific config information
    process_conf = copy.deepcopy(config["target_process"]["QCD"])
    # change cuts from SRlike/ARlike to SR/AR
    process_conf["SRlike_cuts"]["tau_pair_sign"] = "opposite"
    process_conf["ARlike_cuts"]["tau_pair_sign"] = "opposite"
    correction_conf = corr_config["target_process"]["QCD"]["DR_SR"]

    for sample_path in sample_path_list:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        print(
            "Processing {sample} for the DR to SR correction for {process}.".format(
                sample=sample,
                process="QCD",
            )
        )
        print("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for QCD signal-like region
        region_cut_conf = copy.deepcopy(process_conf["SRlike_cuts"])
        rdf_SRlike = region_filter(rdf, config["channel"], region_cut_conf, sample)
        print(
            "Filtering events for the signal-like region. Target process: {}\n".format(
                "QCD"
            )
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike.Report().Print()
        print(out.getvalue())
        print("-" * 50)

        # event filter for QCD application-like region
        region_cut_conf = copy.deepcopy(process_conf["ARlike_cuts"])
        rdf_ARlike = region_filter(rdf, config["channel"], region_cut_conf, sample)
        print(
            "Filtering events for the application-like region. Target process: {}\n".format(
                "QCD"
            )
        )

        # evaluate the measured fake factors for the specific processes
        if sample == "data":
            # rdf_ARlike = func.eval_QCD_FF(rdf_ARlike, config, for_correction=True)
            rdf_ARlike = evaluator.evaluate_tau_pt_njets(rdf_ARlike)
            rdf_ARlike = corr_evaluator.evaluate_lep_pt(rdf_ARlike)
            rdf_ARlike = rdf_ARlike.Define(
                "weight_ff", "weight * QCD_fake_factor * QCD_ff_corr"
            )

        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike.Report().Print()
        print(out.getvalue())
        print("-" * 50)

        # get binning of the dependent variable
        xbinning = array.array("d", correction_conf["var_bins"])
        nbinsx = len(correction_conf["var_bins"]) - 1

        # making the histograms
        h = rdf_SRlike.Histo1D(
            (correction_conf["var_dependence"], "{}".format(sample), nbinsx, xbinning),
            correction_conf["var_dependence"],
            "weight",
        )
        SRlike_hists[sample] = h.GetValue()

        h = rdf_ARlike.Histo1D(
            ("#phi(#tau_{h})", "{}".format(sample), 1, -3.5, 3.5), "phi_2", "weight"
        )
        ARlike_hists[sample] = h.GetValue()

        if sample == "data":
            h = rdf_ARlike.Histo1D(
                (
                    correction_conf["var_dependence"],
                    "{}_ff".format(sample),
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

    corr_hist, proc_frac = func.calculate_non_closure_correction(
        SRlike_hists, ARlike_hists
    )

    smooth_graph, corr_def = func.smooth_function(
        corr_hist.Clone(), correction_conf["var_bins"]
    )

    plotting.plot_correction(
        corr_hist,
        smooth_graph,
        correction_conf["var_dependence"],
        "QCD",
        "DR_SR",
        config,
        save_path,
    )

    plot_hists = dict()

    plot_hists["data_subtracted"] = SRlike_hists["data_subtracted"].Clone()
    plot_hists["data_ff"] = ARlike_hists["data_ff"].Clone()
    plot_hists["data_ff"].Scale(proc_frac)

    data = "data_subtracted"
    samples = ["data_ff"]

    plotting.plot_data_mc(
        plot_hists,
        config,
        correction_conf["var_dependence"],
        "QCD",
        "DR_SR",
        data,
        samples,
        {"incl": ""},
        save_path,
    )

    return corr_def
