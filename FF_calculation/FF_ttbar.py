"""
Function for calculating fake factors for the ttbar process
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


def calculation_ttbar_FFs(config, sample_path_list, save_path):
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
    cs_expressions = dict()

    # get QCD specific config information
    process_conf = config["target_process"]["ttbar"]

    split_vars, split_combinations = func.get_split_combinations(
        process_conf["split_categories"]
    )

    # calculating global histograms for the data/mc scale factor
    for sample_path in sample_path_list:
        # getting the name of the process from the sample path
        sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
        print("Processing {} for the ttbar global data/mc scale factor.".format(sample))
        print("-" * 50)

        rdf = ROOT.RDataFrame(config["tree"], sample_path)

        # event filter for ttbar signal-like region
        region_cut_conf = {**process_conf["SRlike_cuts"]}
        rdf_SRlike = region_filter(rdf, config["channel"], region_cut_conf, sample)
        print(
            "Filtering events for the signal-like region. Target process: {}\n".format(
                "ttbar"
            )
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike.Report().Print()
        print(out.getvalue())
        print("-" * 50)

        # QCD estimation from same sign in signal-like region
        region_cut_conf["tau_pair_sign"] = "same"
        rdf_SRlike_qcd = region_filter(rdf, config["channel"], region_cut_conf, sample)
        print(
            "Filtering events for QCD estimation in the signal-like region. Target process: {}\n".format(
                "ttbar"
            )
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_SRlike_qcd.Report().Print()
        print(out.getvalue())
        print("-" * 50)

        # event filter for ttbar application-like region
        region_cut_conf = {**process_conf["ARlike_cuts"]}
        rdf_ARlike = region_filter(rdf, config["channel"], region_cut_conf, sample)
        print(
            "Filtering events for the application-like region. Target process: {}\n".format(
                "ttbar"
            )
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike.Report().Print()
        print(out.getvalue())
        print("-" * 50)

        # QCD estimation from same sign in application-like region
        region_cut_conf["tau_pair_sign"] = "same"
        rdf_ARlike_qcd = region_filter(rdf, config["channel"], region_cut_conf, sample)
        print(
            "Filtering events for QCD estimation in the application-like region. Target process: {}\n".format(
                "Wjets"
            )
        )
        # redirecting C++ stdout for Report() to python stdout
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            rdf_ARlike_qcd.Report().Print()
        print(out.getvalue())
        print("-" * 50)

        # make yield histograms for FF data correction
        h = rdf_SRlike.Histo1D(
            ("#phi(#tau_{h})", "{}".format(sample), 1, -3.5, 3.5), "phi_2", "weight"
        )
        SRlike_hists[sample] = h.GetValue()

        h = rdf_ARlike.Histo1D(
            ("#phi(#tau_{h})", "{}".format(sample), 1, -3.5, 3.5), "phi_2", "weight"
        )
        ARlike_hists[sample] = h.GetValue()

        # make yield histograms for QCD estimation
        h_qcd = rdf_SRlike_qcd.Histo1D(
            ("#phi(#tau_{h})", "{}".format(sample), 1, -3.5, 3.5), "phi_2", "weight"
        )
        SRlike_hists_qcd[sample] = h_qcd.GetValue()

        h_qcd = rdf_ARlike_qcd.Histo1D(
            ("#phi(#tau_{h})", "{}".format(sample), 1, -3.5, 3.5), "phi_2", "weight"
        )
        ARlike_hists_qcd[sample] = h_qcd.GetValue()

    # calculate QCD estimation
    SRlike_hists["QCD"] = func.QCD_SS_estimate(SRlike_hists_qcd)
    ARlike_hists["QCD"] = func.QCD_SS_estimate(ARlike_hists_qcd)

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
        for sample_path in sample_path_list:
            # getting the name of the process from the sample path
            sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
            # FFs for ttbar from mc -> only ttbar with true misindentified jets relevant
            if sample in ["ttbar_J"]:
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

                # event filter for ttbar signal region
                region_cut_conf = {**split, **process_conf["SR_cuts"]}
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
                region_cut_conf = {**split, **process_conf["AR_cuts"]}
                rdf_AR = region_filter(rdf, config["channel"], region_cut_conf, sample)
                print(
                    "Filtering events for the application region. Target process: {}\n".format(
                        "ttbar"
                    )
                )
                # redirecting C++ stdout for Report() to python stdout
                out = StringIO()
                with pipes(stdout=out, stderr=STDOUT):
                    rdf_AR.Report().Print()
                print(out.getvalue())
                print("-" * 50)

                # get binning of the dependent variable
                xbinning = array.array("d", process_conf["var_bins"])
                nbinsx = len(process_conf["var_bins"]) - 1

                # making the histograms
                h = rdf_SR.Histo1D(
                    (
                        process_conf["var_dependence"],
                        "{}".format(sample),
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
                        "{}".format(sample),
                        nbinsx,
                        xbinning,
                    ),
                    process_conf["var_dependence"],
                    "weight",
                )
                AR_hists[sample] = h.GetValue()

        # Start of the FF calculation
        FF_hist = func.calculate_ttbar_FF(
            SR_hists, AR_hists, SRlike_hists, ARlike_hists
        )
        # performing the fit and calculating the uncertainties
        fit_graphs, cs_exp = func.fit_function(
            FF_hist.Clone(), process_conf["var_bins"]
        )

        plotting.plot_FFs(FF_hist, fit_graphs, "ttbar", config, split, save_path)

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
            SRlike_hists,
            config,
            "phi_2",
            "ttbar",
            "SR_like",
            data,
            samples,
            split,
            save_path,
        )
        plotting.plot_data_mc_ratio(
            ARlike_hists,
            config,
            "phi_2",
            "ttbar",
            "AR_like",
            data,
            samples,
            split,
            save_path,
        )

        data = "data_subtracted"
        samples = ["ttbar_J"]

        plotting.plot_data_mc_ratio(
            SRlike_hists,
            config,
            "phi_2",
            "ttbar",
            "SR_like",
            data,
            samples,
            split,
            save_path,
        )
        plotting.plot_data_mc_ratio(
            ARlike_hists,
            config,
            "phi_2",
            "ttbar",
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
            rdf_AR = evaluator.evaluate_tau_pt_njets(rdf_AR)
            if corr_evaluator == None:
                rdf_AR = rdf_AR.Define("weight_ff", "weight * ttbar_fake_factor")
            else:
                rdf_AR = corr_evaluator.evaluate_lep_pt(rdf_AR)
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
