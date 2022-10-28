import ROOT
import argparse
import yaml
import glob
import sys
import array

import helper.region_filters as filters
import helper.plotting as plotting
from helper.qcd_estimation import QCD_SS_estimate

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config",
    default="hist_config",
    help="Name of a config file in configs/ which contains information for the preselection step.",
)


if __name__ == "__main__":
    args = parser.parse_args()

    # loading of the chosen config file
    with open("configs/" + args.config + ".yaml", "r") as file:
        config = yaml.load(file, yaml.FullLoader)

    # getting all the input files
    sample_path_list = glob.glob(config["file_path"] + "/preselection/" + config["era"] + "/" + config["channel"] + "/*.root")
    print(
        "The following files are loaded for era: {}, channel: {}".format(
            config["era"], config["channel"]
        )
    )
    print("-" * 50)
    for f in sample_path_list:
        print(f)
    print("-" * 50)
    
    for process in config["target_process"]:
        SRlike_hists = dict()
        ARlike_hists = dict()
        if process == "Wjets":
            SRlike_hists_qcd = dict()
            ARlike_hists_qcd = dict()
        process_conf = config["target_process"][process]
        for njets in process_conf["jet_split_categories"]:
            for sample_path in sample_path_list:
                # getting the name of the process from the sample path
                sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
                print("Processing {} for the {} jets category.".format(sample, njets))
                print("-" * 50)

                rdf = ROOT.RDataFrame(config["tree"], sample_path)
                if process == "QCD":
                    rdf_SRlike = filters.signal_like_region_QCD(rdf, config["channel"], njets)
                    rdf_ARlike = filters.application_like_region_QCD(rdf, config["channel"], njets)
                elif process == "Wjets":
                    rdf_SRlike = filters.signal_like_region_Wjets(rdf, config["channel"], njets)
                    rdf_ARlike = filters.application_like_region_Wjets(rdf, config["channel"], njets)
                    rdf_SRlike_qcd = filters.same_sign(rdf_SRlike, config["channel"])
                    rdf_ARlike_qcd = filters.same_sign(rdf_ARlike, config["channel"])
                    rdf_SRlike = filters.opposite_sign(rdf_SRlike, config["channel"])
                    rdf_ARlike = filters.opposite_sign(rdf_ARlike, config["channel"])
                else:
                    sys.exit("Target region: Such a region is not defined: {}".format(process))

                print("Filtering events for the signal-like region. Target process: {}\n".format(process))
                rdf_SRlike.Report().Print()
                print("-" * 50)
                print("Filtering events for the application-like region. Target process: {}\n".format(process))
                rdf_ARlike.Report().Print()
                print("-" * 50)

                xbinning = array.array("d", process_conf["tau_pt_bins"])
                nbinsx = len(process_conf["tau_pt_bins"]) - 1
                h = rdf_SRlike.Histo1D(("p_{T}(#tau_{h})", "e#tau_{h} channel;p_{T}(#tau_{h}) (GeV);FF_{QCD}", nbinsx, xbinning), "pt_2", "weight")
                SRlike_hists[sample] = h.GetValue()
                h = rdf_ARlike.Histo1D(("p_{T}(#tau_{h})", "e#tau_{h} channel;p_{T}(#tau_{h}) (GeV);FF_{QCD}", nbinsx, xbinning), "pt_2", "weight")
                ARlike_hists[sample] = h.GetValue()

                if process == "Wjets":
                    h_qcd = rdf_SRlike_qcd.Histo1D(("p_{T}(#tau_{h})", "e#tau_{h} channel", nbinsx, xbinning), "pt_2", "weight")
                    SRlike_hists_qcd[sample] = h_qcd.GetValue()
                    h_qcd = rdf_ARlike_qcd.Histo1D(("p_{T}(#tau_{h})", "e#tau_{h} channel;", nbinsx, xbinning), "pt_2", "weight")
                    ARlike_hists_qcd[sample] = h_qcd.GetValue()
        
            SRlike_hists["data_substracted"] = SRlike_hists["data"].Clone()
            ARlike_hists["data_substracted"] = ARlike_hists["data"].Clone()

            if process == "Wjets":
                SRlike_hists["QCD"] = QCD_SS_estimate(SRlike_hists_qcd)
                ARlike_hists["QCD"] = QCD_SS_estimate(ARlike_hists_qcd)
            
            # calculating data minus all mc backgrounds except the target process
            for hist in SRlike_hists:
                if hist not in ["data", "data_substracted", process] and "_T" not in hist:
                    SRlike_hists["data_substracted"].Add(SRlike_hists[hist], -1)
            for hist in ARlike_hists:
                if hist not in ["data", "data_substracted", process] and "_T" not in hist:
                    ARlike_hists["data_substracted"].Add(ARlike_hists[hist], -1)
            
            plotting.plot_FFs(SRlike_hists["data_substracted"], ARlike_hists["data_substracted"], config, process, njets)
            
            if process == "Wjets":
                plotting.plot_histogram(SRlike_hists, config, "SR_like", process, njets)
                plotting.plot_histogram(ARlike_hists, config, "AR_like", process, njets)
                plotting.plot_full_histogram(SRlike_hists, config, "SR_like", process, njets)
                plotting.plot_full_histogram(ARlike_hists, config, "AR_like", process, njets)
                plotting.plot_ratio_histogram(SRlike_hists, config, "SR_like", process, njets)
                plotting.plot_ratio_histogram(ARlike_hists, config, "AR_like", process, njets)