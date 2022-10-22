import ROOT
import argparse
import yaml
import glob
import sys
import array

import helper.region_filters as filters
import helper.plotting as plotting

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

    SRlike_hists = dict()
    ARlike_hists = dict()
    
    for njets in config["jet_split_categories"]:
        for sample_path in sample_path_list:
            # getting the name of the process from the sample path
            sample = sample_path.rsplit("/")[-1].rsplit(".")[0]
            print("Processing {} for the {} jets category.".format(sample, njets))
            print("-" * 50)

            rdf = ROOT.RDataFrame(config["tree"], sample_path)
            if config["target_process"] == "QCD":
                rdf_SRlike = filters.signal_like_region_QCD(rdf, config["channel"], njets)
                rdf_ARlike = filters.application_like_region_QCD(rdf, config["channel"], njets)
            else:
                sys.exit("Target region: Such a region is not defined: {}".format(config["target_process"]))

            print("Filtering events for the signal-like region. Target process: {}\n".format(config["target_process"]))
            rdf_SRlike.Report().Print()
            print("-" * 50)
            print("Filtering events for the application-like region. Target process: {}\n".format(config["target_process"]))
            rdf_ARlike.Report().Print()
            print("-" * 50)

            xbinning = array.array("d", [30,35,40,50,500])
            nbinsx = 4
            h = rdf_SRlike.Histo1D(("p_{T}(#tau_{h})", "e#tau_{h} channel;p_{T}(#tau_{h}) (GeV);FF_{QCD}", nbinsx, xbinning), "pt_2", "weight")
            SRlike_hists[sample] = h.GetValue()
            h = rdf_ARlike.Histo1D(("p_{T}(#tau_{h})", "e#tau_{h} channel;p_{T}(#tau_{h}) (GeV);FF_{QCD}", nbinsx, xbinning), "pt_2", "weight")
            ARlike_hists[sample] = h.GetValue()
    
        SRlike_hists["data_substracted"] = SRlike_hists["data"].Clone()
        ARlike_hists["data_substracted"] = ARlike_hists["data"].Clone()
        
        # calculating data minus all mc backgrounds except the target process
        for hist in SRlike_hists:
            if hist not in ["data", "data_substracted", config["target_process"]] and "_T" not in hist:
                SRlike_hists["data_substracted"].Add(SRlike_hists[hist], -1)
        for hist in ARlike_hists:
            if hist not in ["data", "data_substracted", config["target_process"]] and "_T" not in hist:
                ARlike_hists["data_substracted"].Add(ARlike_hists[hist], -1)
        
        plotting.plot_FFs(SRlike_hists["data_substracted"], ARlike_hists["data_substracted"], config["workdir_name"], config["target_process"], njets)
        plotting.plot_histogram(SRlike_hists, config["workdir_name"], config["target_process"], "SR_like", njets)
        plotting.plot_histogram(ARlike_hists, config["workdir_name"], config["target_process"], "AR_like", njets)