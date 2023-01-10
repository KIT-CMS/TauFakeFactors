"""
Function for calculating fake factors for the QCD process
"""

import sys
import array
import ROOT
from io import StringIO
from wurlitzer import pipes, STDOUT

import helper.functions as func
import helper.plotting as plotting

from FF_calculation.FF_region_filters import region_filter


def calculation_QCD_FFs(config, sample_path_list):
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
            if hist not in ["data", "data_subtracted", "data_subtracted_up", "data_subtracted_down", "QCD"]:
                SRlike_hists["data_subtracted"].Add(SRlike_hists[hist], -1)
                SRlike_hists["data_subtracted_up"].Add(SRlike_hists[hist], -0.93)
                SRlike_hists["data_subtracted_down"].Add(SRlike_hists[hist], -1.07)
        for hist in ARlike_hists:
            if hist not in ["data", "data_subtracted", "data_subtracted_up", "data_subtracted_down", "QCD"]:
                ARlike_hists["data_subtracted"].Add(ARlike_hists[hist], -1)
                ARlike_hists["data_subtracted_up"].Add(ARlike_hists[hist], -0.93)
                ARlike_hists["data_subtracted_down"].Add(ARlike_hists[hist], -1.07)

        # Start of the FF calculation
        FF_hist, FF_hist_up, FF_hist_down = func.calculate_QCD_FF(SRlike_hists, ARlike_hists)
        # the fit is performed during the plotting
        cs_exp = plotting.plot_FFs(
            [FF_hist, FF_hist_up, FF_hist_down], "QCD", config, split
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
        )
        print("-" * 50)

    return cs_expressions
