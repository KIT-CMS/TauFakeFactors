"""
Function for calculating the process fractions for the fake factors
"""

import array
import ROOT
from io import StringIO
from wurlitzer import pipes, STDOUT

import helper.functions as func
import helper.plotting as plotting

from FF_calculation.FF_region_filters import region_filter


def fraction_calculation(config, sample_path_list, save_path):
    # init histogram dict for the fraction calculation
    AR_hists = dict()
    SR_hists = dict()
    fractions = dict()

    # get config information for the fraction calculation
    process_conf = config["process_fractions"]

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

            # event filter for application region
            region_cut_conf = {**split, **process_conf["AR_cuts"]}
            rdf_AR = region_filter(rdf, config["channel"], region_cut_conf, sample)
            print(
                "Filtering events for the fraction calculation in the application region.\n"
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_AR.Report().Print()
            print(out.getvalue())
            print("-" * 50)

            # event filter for signal region; this is not needed for the FF calculation, just for control plots
            region_cut_conf = {**split, **process_conf["SR_cuts"]}
            rdf_SR = region_filter(rdf, config["channel"], region_cut_conf, sample)
            print(
                "Filtering events for the fraction calculation in the signal region.\n"
            )
            # redirecting C++ stdout for Report() to python stdout
            out = StringIO()
            with pipes(stdout=out, stderr=STDOUT):
                rdf_SR.Report().Print()
            print(out.getvalue())
            print("-" * 50)

            # get binning of the dependent variable
            xbinning = array.array("d", process_conf["var_bins"])
            nbinsx = len(process_conf["var_bins"]) - 1

            # making the histograms
            h = rdf_AR.Histo1D(
                (process_conf["var_dependence"], "{}".format(sample), nbinsx, xbinning),
                process_conf["var_dependence"],
                "weight",
            )
            AR_hists[sample] = h.GetValue()

            h = rdf_SR.Histo1D(
                (process_conf["var_dependence"], "{}".format(sample), nbinsx, xbinning),
                process_conf["var_dependence"],
                "weight",
            )
            SR_hists[sample] = h.GetValue()

        # calculate QCD estimation; here directly estimated as difference between mc and data without SS/OS
        AR_hists["QCD"] = func.QCD_SS_estimate(AR_hists)
        SR_hists["QCD"] = func.QCD_SS_estimate(SR_hists)

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

        frac_hists = dict()

        for p in config["process_fractions"]["processes"]:
            frac_hists[p] = func.calc_fraction(
                AR_hists, p, config["process_fractions"]["processes"]
            )
        frac_hists = func.add_fraction_variations(
            frac_hists, config["process_fractions"]["processes"]
        )

        cat = "#".join(["{}#{}".format(split[var], var) for var in split_vars])
        fractions[cat] = frac_hists

        plotting.fraction_plot(
            frac_hists["nominal"],
            config,
            process_conf["var_dependence"],
            "AR",
            config["process_fractions"]["processes"],
            split,
            save_path,
        )
        plotting.plot_data_mc_ratio(
            AR_hists,
            config,
            process_conf["var_dependence"],
            "fraction",
            "AR",
            data,
            samples,
            split,
            save_path,
        )
        plotting.plot_data_mc_ratio(
            SR_hists,
            config,
            process_conf["var_dependence"],
            "fraction",
            "SR",
            data,
            samples,
            split,
            save_path,
        )
        print("-" * 50)

    return fractions
