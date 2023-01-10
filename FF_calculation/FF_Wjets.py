"""
Function for calculating fake factors for the W-jets process
"""

import sys
import array
import ROOT
from io import StringIO
from wurlitzer import pipes, STDOUT

import helper.functions as func
import helper.plotting as plotting

from FF_calculation.FF_region_filters import region_filter


def calculation_Wjets_FFs(config, sample_path_list):
    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()
    # init histogram dict for QCD SS/OS estimation
    SRlike_hists_qcd = dict()
    ARlike_hists_qcd = dict()
    # init dictionary for the FF functions for correctionlib
    cs_expressions = dict()

    # get QCD specific config information
    process_conf = config["target_process"]["Wjets"]

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

            # event filter for Wjets signal-like region
            region_cut_conf = {**split, **process_conf["SRlike_cuts"]}
            rdf_SRlike = region_filter(rdf, config["channel"], region_cut_conf, sample)
            print(
                "Filtering events for the signal-like region. Target process: {}\n".format(
                    "Wjets"
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
            rdf_SRlike_qcd = region_filter(
                rdf, config["channel"], region_cut_conf, sample
            )
            print(
                "Filtering events for QCD estimation in the signal-like region. Target process: {}\n".format(
                    "Wjets"
                )
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_SRlike_qcd.Report().Print()
            print(out.getvalue())
            print("-" * 50)

            # event filter for Wjets application-like region
            region_cut_conf = {**split, **process_conf["ARlike_cuts"]}
            rdf_ARlike = region_filter(rdf, config["channel"], region_cut_conf, sample)
            print(
                "Filtering events for the application-like region. Target process: {}\n".format(
                    "Wjets"
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
            rdf_ARlike_qcd = region_filter(
                rdf, config["channel"], region_cut_conf, sample
            )
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

            # making the histograms for QCD estimation
            h_qcd = rdf_SRlike_qcd.Histo1D(
                (process_conf["var_dependence"], "{}".format(sample), nbinsx, xbinning),
                process_conf["var_dependence"],
                "weight",
            )
            SRlike_hists_qcd[sample] = h_qcd.GetValue()

            h_qcd = rdf_ARlike_qcd.Histo1D(
                (process_conf["var_dependence"], "{}".format(sample), nbinsx, xbinning),
                process_conf["var_dependence"],
                "weight",
            )
            ARlike_hists_qcd[sample] = h_qcd.GetValue()

        # calculate QCD estimation
        SRlike_hists["QCD"] = func.QCD_SS_estimate(SRlike_hists_qcd)
        ARlike_hists["QCD"] = func.QCD_SS_estimate(ARlike_hists_qcd)

        # calculate Wjets enriched data by subtraction all there backgrould sample
        SRlike_hists["data_subtracted"] = SRlike_hists["data"].Clone()
        ARlike_hists["data_subtracted"] = ARlike_hists["data"].Clone()
        SRlike_hists["data_subtracted_up"] = SRlike_hists["data"].Clone()
        ARlike_hists["data_subtracted_up"] = ARlike_hists["data"].Clone()
        SRlike_hists["data_subtracted_down"] = SRlike_hists["data"].Clone()
        ARlike_hists["data_subtracted_down"] = ARlike_hists["data"].Clone()

        for hist in SRlike_hists:
            if hist not in ["data", "data_subtracted", "data_subtracted_up", "data_subtracted_down", "Wjets"]:
                SRlike_hists["data_subtracted"].Add(SRlike_hists[hist], -1)
                SRlike_hists["data_subtracted_up"].Add(SRlike_hists[hist], -0.93)
                SRlike_hists["data_subtracted_down"].Add(SRlike_hists[hist], -1.07)
        for hist in ARlike_hists:
            if hist not in ["data", "data_subtracted", "data_subtracted_up", "data_subtracted_down", "Wjets"]:
                ARlike_hists["data_subtracted"].Add(ARlike_hists[hist], -1)
                ARlike_hists["data_subtracted_up"].Add(ARlike_hists[hist], -0.93)
                ARlike_hists["data_subtracted_down"].Add(ARlike_hists[hist], -1.07)

        # Start of the FF calculation
        FF_hist, FF_hist_up, FF_hist_down = func.calculate_Wjets_FF(SRlike_hists, ARlike_hists)
        # the fit is performed during the plotting
        cs_exp = plotting.plot_FFs(
            [FF_hist, FF_hist_up, FF_hist_down], "Wjets", config, split
        )  
        
        if len(split) == 1:
            cs_expressions["{}#{}".format(split_vars[0],split[split_vars[0]])] = cs_exp
        elif len(split) == 2:
            if "{}#{}".format(split_vars[0],split[split_vars[0]]) not in cs_expressions:
                cs_expressions["{}#{}".format(split_vars[0],split[split_vars[0]])] = dict()
            cs_expressions["{}#{}".format(split_vars[0],split[split_vars[0]])]["{}#{}".format(split_vars[1],split[split_vars[1]])] = cs_exp
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
            ]

        plotting.plot_data_mc_ratio(
            SRlike_hists,
            config,
            process_conf["var_dependence"],
            "Wjets",
            "SR_like",
            data,
            samples,
            split,
        )
        plotting.plot_data_mc_ratio(
            ARlike_hists,
            config,
            process_conf["var_dependence"],
            "Wjets",
            "AR_like",
            data,
            samples,
            split,
        )

        data = "data_subtracted"
        samples = ["Wjets"]

        plotting.plot_data_mc_ratio(
            SRlike_hists,
            config,
            process_conf["var_dependence"],
            "Wjets",
            "SR_like",
            data,
            samples,
            split,
        )
        plotting.plot_data_mc_ratio(
            ARlike_hists,
            config,
            process_conf["var_dependence"],
            "Wjets",
            "AR_like",
            data,
            samples,
            split,
        )
        print("-" * 50)

    return cs_expressions
