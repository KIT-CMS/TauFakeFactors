"""
Function for calculating fake factors for the QCD process
"""

import array
import ROOT

from helper import region_filters as filters
import helper.functions as func
import helper.plotting as plotting


def calculation_QCD_FFs(config, sample_path_list):
    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    process_conf = config["target_process"]["QCD"]

    # differentiating between N jets categories
    for njets in process_conf["njet_split_categories"]:
        for sample_path in sample_path_list:
            # getting the name of the process from the sample path
            sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
            print("Processing {} for the {} jets category.".format(sample, njets))
            print("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)

            # event filter for QCD signal-like region
            region_cut_conf = process_conf["SRlike_cuts"]
            rdf_SRlike = filters.same_opposite_sign(
                rdf, config["channel"], region_cut_conf["tau_pair_sign"]
            )
            rdf_SRlike = filters.lepton_mT(
                rdf_SRlike, config["channel"], region_cut_conf
            )
            rdf_SRlike = filters.no_extra_leptons(
                rdf_SRlike, config["channel"], region_cut_conf
            )
            rdf_SRlike = filters.tau_id_vs_jets_WP(
                rdf_SRlike, config["channel"], region_cut_conf
            )
            rdf_SRlike = filters.jet_number(rdf_SRlike, config["channel"], njets)

            print(
                "Filtering events for the signal-like region. Target process: {}\n".format(
                    "QCD"
                )
            )
            rdf_SRlike.Report().Print()
            print("-" * 50)

            # event filter for QCD application-like region
            region_cut_conf = process_conf["ARlike_cuts"]
            rdf_ARlike = filters.same_opposite_sign(
                rdf, config["channel"], region_cut_conf["tau_pair_sign"]
            )
            rdf_ARlike = filters.lepton_mT(
                rdf_ARlike, config["channel"], region_cut_conf
            )
            rdf_ARlike = filters.no_extra_leptons(
                rdf_ARlike, config["channel"], region_cut_conf
            )
            rdf_ARlike = filters.tau_id_vs_jets_between_WPs(
                rdf_ARlike, config["channel"], region_cut_conf
            )
            rdf_ARlike = filters.jet_number(rdf_ARlike, config["channel"], njets)

            print(
                "Filtering events for the application-like region. Target process: {}\n".format(
                    "QCD"
                )
            )
            rdf_ARlike.Report().Print()
            print("-" * 50)

            # get binning for tau pT
            xbinning = array.array("d", process_conf["tau_pt_bins"])
            nbinsx = len(process_conf["tau_pt_bins"]) - 1
            # make tau pT histograms
            h = rdf_SRlike.Histo1D(
                ("p_{T}(#tau_{h})", "{}".format(sample), nbinsx, xbinning),
                "pt_2",
                "weight",
            )
            SRlike_hists[sample] = h.GetValue()
            h = rdf_ARlike.Histo1D(
                ("p_{T}(#tau_{h})", "{}".format(sample), nbinsx, xbinning),
                "pt_2",
                "weight",
            )
            ARlike_hists[sample] = h.GetValue()

        # calculate Wjets enriched data by subtraction all there backgrould sample
        SRlike_hists["data_subtracted"] = SRlike_hists["data"].Clone()
        ARlike_hists["data_subtracted"] = ARlike_hists["data"].Clone()

        for hist in SRlike_hists:
            if hist not in ["data", "data_subtracted", "QCD"] and "_T" not in hist:
                SRlike_hists["data_subtracted"].Add(SRlike_hists[hist], -1)
        for hist in ARlike_hists:
            if hist not in ["data", "data_subtracted", "QCD"] and "_T" not in hist:
                ARlike_hists["data_subtracted"].Add(ARlike_hists[hist], -1)

        # Start of the calculation
        FF_hist = func.calculate_QCD_FF(SRlike_hists, ARlike_hists)
        plotting.plot_FFs(FF_hist, config, "QCD", njets)

        sig = "data"
        bkg = [
            "diboson_J",
            "diboson_L",
            "Wjets",
            "ttbar_J",
            "ttbar_L",
            "DYjets_J",
            "DYjets_L",
            "embedding",
        ]
        plotting.plot_data_mc(
            SRlike_hists, config, "SR_like", "tau_pt", "QCD", njets, sig, bkg
        )
        plotting.plot_data_mc(
            ARlike_hists, config, "AR_like", "tau_pt", "QCD", njets, sig, bkg
        )
