"""
Function for calculating fake factors for the QCD process
"""

import array
import ROOT

from helper import region_filters as filters
import helper.plotting as plotting

def calculation_QCD_FFs(config, sample_path_list):
    # init histogram dict for FF measurement
    SRlike_hists = dict()
    ARlike_hists = dict()

    process_conf = config["target_process"]["QCD"]

    # differentiating between N_jets categories
    for njets in process_conf["njet_split_categories"]:
        for sample_path in sample_path_list:
            # getting the name of the process from the sample path
            sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
            print("Processing {} for the {} jets category.".format(sample, njets))
            print("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)

            # event filter for QCD region
            rdf_SRlike = filters.signal_like_region_QCD(rdf, config["channel"], njets)
            rdf_ARlike = filters.application_like_region_QCD(rdf, config["channel"], njets)

            print("Filtering events for the signal-like region. Target process: {}\n".format("QCD"))
            rdf_SRlike.Report().Print()
            print("-" * 50)
            print("Filtering events for the application-like region. Target process: {}\n".format("QCD"))
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
    
        # calculate Wjets enriched data by subtraction all there backgrould sample
        SRlike_hists["data_substracted"] = SRlike_hists["data"].Clone()
        ARlike_hists["data_substracted"] = ARlike_hists["data"].Clone()
        
        for hist in SRlike_hists:
            if hist not in ["data", "data_substracted", "QCD"] and "_T" not in hist:
                SRlike_hists["data_substracted"].Add(SRlike_hists[hist], -1)
        for hist in ARlike_hists:
            if hist not in ["data", "data_substracted", "QCD"] and "_T" not in hist:
                ARlike_hists["data_substracted"].Add(ARlike_hists[hist], -1)
        
        # Start of the calculation
        plotting.plot_FFs(SRlike_hists["data_substracted"], ARlike_hists["data_substracted"], config, "QCD", njets)