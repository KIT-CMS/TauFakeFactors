"""
Function for calculating fake factors for the W-jets process
"""

import array
import ROOT

from helper import region_filters as filters
import helper.functions as func
import helper.plotting as plotting

def calculation_Wjets_FFs(config, sample_path_list):
    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()
    # init histogram dict for QCD SS/OS estimation
    SRlike_hists_qcd = dict()
    ARlike_hists_qcd = dict()

    process_conf = config["target_process"]["Wjets"]

    # differentiating between N_jets categories
    for njets in process_conf["njet_split_categories"]:
        for sample_path in sample_path_list:
            # getting the name of the process from the sample path
            sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
            print("Processing {} for the {} jets category.".format(sample, njets))
            print("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)
            
            # default event filter for Wjets region
            rdf_SRlike = filters.signal_like_region_Wjets(rdf, config["channel"], njets)
            rdf_ARlike = filters.application_like_region_Wjets(rdf, config["channel"], njets)
            # split into same/opposite sign regions for QCD estimation
            rdf_SRlike_qcd = filters.same_sign(rdf_SRlike, config["channel"])
            rdf_ARlike_qcd = filters.same_sign(rdf_ARlike, config["channel"])
            rdf_SRlike = filters.opposite_sign(rdf_SRlike, config["channel"])
            rdf_ARlike = filters.opposite_sign(rdf_ARlike, config["channel"])

            print("Filtering events for the signal-like region. Target process: {}\n".format("Wjets"))
            rdf_SRlike.Report().Print()
            print("-" * 50)
            print("Filtering events for the application-like region. Target process: {}\n".format("Wjets"))
            rdf_ARlike.Report().Print()
            print("-" * 50)
            
            # get binning for tau pT
            xbinning = array.array("d", process_conf["tau_pt_bins"])
            nbinsx = len(process_conf["tau_pt_bins"]) - 1
            # make tau pT histograms
            h = rdf_SRlike.Histo1D(("p_{T}(#tau_{h})", "{}".format(sample), nbinsx, xbinning), "pt_2", "weight")
            SRlike_hists[sample] = h.GetValue()
            h = rdf_ARlike.Histo1D(("p_{T}(#tau_{h})", "{}".format(sample), nbinsx, xbinning), "pt_2", "weight")
            ARlike_hists[sample] = h.GetValue()
            # make tau pT histograms for QCD estimation
            h_qcd = rdf_SRlike_qcd.Histo1D(("p_{T}(#tau_{h})", "{}".format(sample), nbinsx, xbinning), "pt_2", "weight")
            SRlike_hists_qcd[sample] = h_qcd.GetValue()
            h_qcd = rdf_ARlike_qcd.Histo1D(("p_{T}(#tau_{h})", "{}".format(sample), nbinsx, xbinning), "pt_2", "weight")
            ARlike_hists_qcd[sample] = h_qcd.GetValue()

        # calculate QCD estimation
        SRlike_hists["QCD"] = func.QCD_SS_estimate(SRlike_hists_qcd)
        ARlike_hists["QCD"] = func.QCD_SS_estimate(ARlike_hists_qcd)

        # calculate Wjets enriched data by subtraction all there backgrould sample
        SRlike_hists["data_substracted"] = SRlike_hists["data"].Clone()
        ARlike_hists["data_substracted"] = ARlike_hists["data"].Clone()
        
        for hist in SRlike_hists:
            if hist not in ["data", "data_substracted", "Wjets"] and "_T" not in hist:
                SRlike_hists["data_substracted"].Add(SRlike_hists[hist], -1)
        for hist in ARlike_hists:
            if hist not in ["data", "data_substracted", "Wjets"] and "_T" not in hist:
                ARlike_hists["data_substracted"].Add(ARlike_hists[hist], -1)
        
        # Start of the calculation
        plotting.plot_FFs(SRlike_hists["data_substracted"], ARlike_hists["data_substracted"], config, "Wjets", njets)
        
        plotting.plot_histogram(SRlike_hists, config, "SR_like", "Wjets", njets)
        plotting.plot_histogram(ARlike_hists, config, "AR_like", "Wjets", njets)
        plotting.plot_full_histogram(SRlike_hists, config, "SR_like", "Wjets", njets)
        plotting.plot_full_histogram(ARlike_hists, config, "AR_like", "Wjets", njets)
        plotting.plot_ratio_histogram(SRlike_hists, config, "SR_like", "Wjets", njets)
        plotting.plot_ratio_histogram(ARlike_hists, config, "AR_like", "Wjets", njets)