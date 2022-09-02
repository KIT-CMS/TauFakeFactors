import ROOT
import argparse
import yaml
import glob
import os
import sys

# import logging
# import logging.handlers
# logger = logging.getLogger()

import helper.filters as filters
import helper.weights as weights

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config",
    default="config",
    help="Name of a config file in configs/ which contains information for the preselection step.",
)


# def setup_logging(output_path):
#     logger.setLevel("INFO")
#     terminal_handler = logging.StreamHandler()
#     terminal_handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
#     logger.addHandler(terminal_handler)
#     handler = logging.handlers.WatchedFileHandler(output_path + "/preselection.log", "w")
#     handler.setFormatter(logging.Formatter(logging.BASIC_FORMAT))
#     logger.addHandler(handler)

def get_samples(config, sample):
    sample_paths = (
        config["ntuple_path"]
        + "/"
        + config["era"]
        + "/"
        + sample
        + "/"
        + config["channel"]
        + "/*.root"
    )
    print(
        "The following files are loaded for era: {}, channel: {}".format(
            config["era"], config["channel"]
        )
    )
    print("-" * 50)
    sample_path_list = glob.glob(sample_paths)
    for f in sample_path_list:
        print(f)
    print("-" * 50)
    if sample_path_list == []:
        sys.exit("Input files: No files found for {}".format(sample))

    return sample_path_list

def check_output_path(output_path):
    if not os.path.exists(output_path):
        print(
            "Output directory does not exist! Making directory {}".format(output_path)
        )
        print("-" * 50)
        os.makedirs(output_path, exist_ok=True)

def get_output_name(output_path, process, tau_gen_mode, idx=None):    
    if tau_gen_mode == "all":
        tau_gen_mode = ""
    else:
        tau_gen_mode = "_" + tau_gen_mode

    if idx is not None:
        return output_path + "/" + process + tau_gen_mode + "_" + str(idx) + ".root"
    else:
        return output_path + "/" + process + tau_gen_mode + ".root"

class Logger(object):
    def __init__(self, filename="Default.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

output_feature = [
                "weight_sf",
                "sf",
                "weight",
                "passesDLVeto",
                "passes3LVeto",
                "bpt_1",
                "bpt_2",
                "njets",
                "nbtag",
                "mjj",
                "met",
                "lep_pt",
                "lep_eta",
                "lep_phi",
                "lep_iso",
                "lep_q",
                "lep_dR",
                "lep_vvvloose",
                "lep_vvloose",
                "lep_vloose",
                "lep_loose",
                "lep_medium",
                "mt_leplep",
                "m_leplep",
                "otherLep_pt",
                "otherLep_eta",
                "otherLep_iso",
                "otherLep_q",
                "alltau_vvvlooseDNN",
                "alltau_vvlooseDNN",
                "alltau_vlooseDNN",
                "alltau_looseDNN",
                "alltau_mediumDNN",
                "alltau_tightDNN",
                "alltau_vtightDNN",
                "alltau_vvtightDNN",
                "alltau_pt",
                "alltau_eta",
                "alltau_phi",
                "alltau_q",
                "alltau_decay",
                "alltau_gen_match",
                "alltau_lepVeto",
                "alltau_mvis",
                "alltau_mt",
                "alltau_mt2",
                "alltau_svfit",
                "alltau_beta",
                "alltau_Zpt",
                "alltau_dRToLep",
                "alltau_dRToOtherLep",
                "alltau_dRToB",
                "tau_iso_ind",
                "n_iso_otherLep",
                "n_iso_lep",
                "jdeta",
                "njetingap20",
            ]


if __name__ == "__main__":
    args = parser.parse_args()

    with open("configs/" + args.config + ".yaml", "r") as file:
        config = yaml.load(file, yaml.FullLoader)
    with open("configs/datasets.yaml", "r") as file:
        datasets = yaml.load(file, yaml.FullLoader)

    output_path = config["output_path"] + "/preselection/" + config["channel"]
    check_output_path(output_path)
    # setup_logging(output_path)
    fil = open(output_path + "preselection.log", "a")
    sys.stdout = fil # Logger(output_path + "preselection.log")

    # going through all wanted processes
    for process in config["processes"]:
        process_file_list = []
        # going through all sample splitting modes, the modes are based on the tau origin (genuine, jet fake, lepton fake)
        for tau_gen_mode in config["processes"][process]["tau_gen_modes"]:
            # going through all ntuple files (samples) for the process
            for idx, sample in enumerate(config["processes"][process]["samples"]):
                rdf = ROOT.RDataFrame(config["tree"], get_samples(config, sample))

                # apply wanted event filters
                rdf = filters.had_tau_pt_cut(rdf, config["channel"])
                rdf = filters.had_tau_eta_cut(rdf, config["channel"])
                rdf = filters.had_tau_decay_mode_cut(rdf, config["channel"])
                rdf = filters.had_tau_vsLep_id_cut(rdf, config["channel"])
                rdf = filters.trigger_cut(rdf, config["channel"])

                if process in ["DYjets", "ttbar", "diboson"]:
                    rdf = filters.tau_gen_match_split(rdf, config["channel"], tau_gen_mode)
                if process == "embedding":
                    rdf = filters.emb_tau_gen_match(rdf, config["channel"])

                rdf.Report().Print()
                # verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)
                #logger.info()
                print("-" * 50)

                # check for empty data frame -> only save/calculate if event number is not zero
                if rdf.Count().GetValue() != 0:
                    # calculate event weights for plotting
                    rdf = rdf.Define("weight", "1.")
                    if process not in ["data", "embedding"]:
                        rdf = weights.iso_weight(rdf, config["channel"])
                        rdf = weights.id_weight(rdf, config["channel"])
                        rdf = weights.tau_id_vsJet_weight(rdf, config["channel"])
                        rdf = weights.tau_id_vsMu_weight(rdf, config["channel"])
                        rdf = weights.tau_id_vsEle_weight(rdf, config["channel"])
                        rdf = weights.lumi_weight(rdf, config["era"])
                        rdf = weights.pileup_weight(rdf)
                        rdf = weights.trigger_weight(rdf, config["channel"], process)
                        rdf = weights.gen_weight(rdf, datasets[sample])
                        rdf = weights.Z_pt_reweight(rdf, process)
                        rdf = weights.Top_pt_reweight(rdf, process)
                    if process == "embedding":
                        rdf = weights.emb_gen_weight(rdf)
                        rdf = weights.iso_weight(rdf, config["channel"])
                        rdf = weights.id_weight(rdf, config["channel"])
                        rdf = weights.trigger_weight(rdf, config["channel"], process)
                    rdf = rdf.Define("weight_sf", "weight")
                    rdf = rdf.Define("sf", "0.")

                    # calculate needed variables
                    if config["channel"] == "tt":
                        rdf = rdf.Define("passesDLVeto", "(id_tau_vsMu_VLoose_1 > 0.5) && (id_tau_vsEle_VVLoose_1 > 0.5) && (id_tau_vsMu_VLoose_2 > 0.5) && (id_tau_vsEle_VVLoose_2 > 0.5)")
                        rdf = rdf.Define("alltau_lepVeto", "(id_tau_vsMu_VLoose_1 > 0.5) && (id_tau_vsEle_VVLoose_1 > 0.5) && (id_tau_vsMu_VLoose_2 > 0.5) && (id_tau_vsEle_VVLoose_2 > 0.5)")
                    elif config["channel"] == "mt":
                        rdf = rdf.Define("passesDLVeto", "(dimuon_veto < 0.5)")
                        rdf = rdf.Define("alltau_lepVeto", "(id_tau_vsMu_Tight_2 > 0.5) && (id_tau_vsEle_VVLoose_2 > 0.5)")
                    elif config["channel"] == "et":
                        rdf = rdf.Define("passesDLVeto", "(dimuon_veto < 0.5)")
                        rdf = rdf.Define("alltau_lepVeto", "(id_tau_vsMu_VLoose_2 > 0.5) && (id_tau_vsEle_Tight_2 > 0.5)")
                    else:
                        print("Dilepton veto: Such a channel is not defined: {}".format(config["channel"]))

                    rdf = rdf.Define("passes3LVeto", "(extramuon_veto < 0.5) && (extraelec_veto < 0.5) && (dimuon_veto < 0.5)")

                    rdf = rdf.Define("lep_pt", "pt_1")
                    rdf = rdf.Define("lep_eta", "eta_1")
                    rdf = rdf.Define("lep_phi", "phi_1")
                    rdf = rdf.Define("lep_iso", "iso_1")
                    rdf = rdf.Define("lep_q", "q_1")

                    rdf = rdf.Define("lep_dR", "-999.")
                    rdf = rdf.Define("mt_leplep", "-999.")
                    rdf = rdf.Define("m_leplep", "-999.")

                    rdf = rdf.Define("lep_vvvloose", "0.")
                    rdf = rdf.Define("lep_vvloose", "0.")
                    rdf = rdf.Define("lep_vloose", "0.")
                    rdf = rdf.Define("lep_loose", "0.")
                    rdf = rdf.Define("lep_medium", "0.")

                    rdf = rdf.Define("otherLep_pt", "-999.")
                    rdf = rdf.Define("otherLep_eta", "-999.")
                    rdf = rdf.Define("otherLep_iso", "-999.")
                    rdf = rdf.Define("otherLep_q", "-999")

                    rdf = rdf.Define("tau_iso_ind", "0")
                    rdf = rdf.Define("n_iso_otherLep", "0")
                    rdf = rdf.Define("n_iso_lep", "iso_1 < 0.15")
                    rdf = rdf.Define("jdeta", "abs(beta_1 - beta_2)")
                    rdf = rdf.Define("njetingap20", "njets*(jpt_1<30)")

                    if config["channel"] == "tt":
                        rdf = rdf.Define("alltau_vvvlooseDNN", "std::vector<Int_t> alltau_vvvlooseDNN{id_tau_vsJet_VVVLoose_1,id_tau_vsJet_VVVLoose_2}; return alltau_vvvlooseDNN;")
                        rdf = rdf.Define("alltau_vvlooseDNN", "std::vector<Int_t> alltau_vvlooseDNN{id_tau_vsJet_VVLoose_1,id_tau_vsJet_VVLoose_2}; return alltau_vvlooseDNN;")
                        rdf = rdf.Define("alltau_vlooseDNN", "std::vector<Int_t> alltau_vlooseDNN{id_tau_vsJet_VLoose_1,id_tau_vsJet_VLoose_2}; return alltau_vlooseDNN;")
                        rdf = rdf.Define("alltau_looseDNN", "std::vector<Int_t> alltau_looseDNN{id_tau_vsJet_Loose_1,id_tau_vsJet_Loose_2}; return alltau_looseDNN;")
                        rdf = rdf.Define("alltau_mediumDNN", "std::vector<Int_t> alltau_mediumDNN{id_tau_vsJet_Medium_1,id_tau_vsJet_Medium_2}; return alltau_mediumDNN;")
                        rdf = rdf.Define("alltau_tightDNN", "std::vector<Int_t> alltau_tightDNN{id_tau_vsJet_Tight_1,id_tau_vsJet_Tight_2}; return alltau_tightDNN;")
                        rdf = rdf.Define("alltau_vtightDNN", "std::vector<Int_t> alltau_vtightDNN{id_tau_vsJet_VTight_1,id_tau_vsJet_VTight_2}; return alltau_vtightDNN;")
                        rdf = rdf.Define("alltau_vvtightDNN", "std::vector<Int_t> alltau_vvtightDNN{id_tau_vsJet_VVTight_1,id_tau_vsJet_VVTight_2}; return alltau_vvtightDNN;")

                        rdf = rdf.Define("alltau_pt", "std::vector<Double_t> alltau_pt{pt_1,pt_2}; return alltau_pt;")
                        rdf = rdf.Define("alltau_eta", "std::vector<Double_t> alltau_eta{eta_1,eta_2}; return alltau_eta;")
                        rdf = rdf.Define("alltau_phi", "std::vector<Double_t> alltau_phi{phi_1,phi_2}; return alltau_phi;")
                        rdf = rdf.Define("alltau_q", "std::vector<Int_t> alltau_q{q_1,q_2}; return alltau_q;")
                        rdf = rdf.Define("alltau_decay", "std::vector<Int_t> alltau_decay{decaymode_1,decaymode_2}; return alltau_decay;")
                        if process != "data":
                            rdf = rdf.Define("alltau_gen_match", "std::vector<Int_t> alltau_gen_match{gen_match_1,gen_match_2}; return alltau_gen_match;")
                        else:
                            rdf = rdf.Define("alltau_gen_match", "-999")
                        rdf = rdf.Define("alltau_mvis", "std::vector<Double_t> alltau_mvis{m_vis,m_vis}; return alltau_mvis;")
                        rdf = rdf.Define("alltau_mt", "std::vector<Double_t> alltau_mt{mt_1,mt_2}; return alltau_mt;")
                        rdf = rdf.Define("alltau_mt2", "std::vector<Double_t> alltau_mt2{mt_2,mt_1}; return alltau_mt2;")
                        rdf = rdf.Define("alltau_svfit", "std::vector<Double_t> alltau_svfit{0.,0.}; return alltau_svfit;")
                        rdf = rdf.Define("alltau_beta", "std::vector<Double_t> alltau_beta{-999.,-999.}; return alltau_beta;")

                        rdf = rdf.Define("v4_tau_1", "ROOT::Math::PtEtaPhiMVector v4_tau_1(pt_1,eta_1,phi_1,mass_1); return v4_tau_1;")
                        rdf = rdf.Define("v4_tau_2", "ROOT::Math::PtEtaPhiMVector v4_tau_2(pt_2,eta_2,phi_2,mass_2); return v4_tau_2;")
                        rdf = rdf.Define("v4_met", "ROOT::Math::PtEtaPhiMVector v4_met(met,0.,metphi,0.); return v4_met;")
                        rdf = rdf.Define("alltau_Zpt", "std::vector<Double_t> alltau_Zpt{(v4_tau_1 + v4_tau_2 + v4_met).Pt(),(v4_tau_1 + v4_tau_2 + v4_met).Pt()}; return alltau_Zpt;")

                        rdf = rdf.Define("alltau_dRToLep", "std::vector<Double_t> alltau_dRToLep{ROOT::VecOps::DeltaR(eta_1, eta_2, phi_1, phi_2),ROOT::VecOps::DeltaR(eta_1, eta_2, phi_1, phi_2)}; return alltau_dRToLep;")
                        rdf = rdf.Define("alltau_dRToOtherLep", "std::vector<Double_t> alltau_dRToOtherLep{1e6,1e6}; return alltau_dRToOtherLep;")

                        rdf = rdf.Define("dR_TauToB_1", "ROOT::VecOps::DeltaR(beta_1, eta_2, bphi_1, phi_2)")
                        rdf = rdf.Define("dR_TauToB_2", "ROOT::VecOps::DeltaR(beta_2, eta_2, bphi_2, phi_2)")
                        rdf = rdf.Define("alltau_dRToB", "std::vector<Double_t> alltau_dRToB{min(dR_TauToB_1,dR_TauToB_2),min(dR_TauToB_1,dR_TauToB_2)}; return alltau_dRToB;")

                    elif config["channel"] == "mt" or config["channel"] == "et":
                        rdf = rdf.Define("alltau_vvvlooseDNN", "std::vector<Int_t> alltau_vvvlooseDNN{id_tau_vsJet_VVVLoose_2}; return alltau_vvvlooseDNN;")
                        rdf = rdf.Define("alltau_vvlooseDNN", "std::vector<Int_t> alltau_vvlooseDNN{id_tau_vsJet_VVLoose_2}; return alltau_vvlooseDNN;")
                        rdf = rdf.Define("alltau_vlooseDNN", "std::vector<Int_t> alltau_vlooseDNN{id_tau_vsJet_VLoose_2}; return alltau_vlooseDNN;")
                        rdf = rdf.Define("alltau_looseDNN", "std::vector<Int_t> alltau_looseDNN{id_tau_vsJet_Loose_2}; return alltau_looseDNN;")
                        rdf = rdf.Define("alltau_mediumDNN", "std::vector<Int_t> alltau_mediumDNN{id_tau_vsJet_Medium_2}; return alltau_mediumDNN;")
                        rdf = rdf.Define("alltau_tightDNN", "std::vector<Int_t> alltau_tightDNN{id_tau_vsJet_Tight_2}; return alltau_tightDNN;")
                        rdf = rdf.Define("alltau_vtightDNN", "std::vector<Int_t> alltau_vtightDNN{id_tau_vsJet_VTight_2}; return alltau_vtightDNN;")
                        rdf = rdf.Define("alltau_vvtightDNN", "std::vector<Int_t> alltau_vvtightDNN{id_tau_vsJet_VVTight_2}; return alltau_vvtightDNN;")

                        rdf = rdf.Define("alltau_pt", "std::vector<Double_t> alltau_pt{pt_2}; return alltau_pt;")
                        rdf = rdf.Define("alltau_eta", "std::vector<Double_t> alltau_eta{eta_2}; return alltau_eta;")
                        rdf = rdf.Define("alltau_phi", "std::vector<Double_t> alltau_phi{phi_2}; return alltau_phi;")
                        rdf = rdf.Define("alltau_q", "std::vector<Int_t> alltau_q{q_2}; return alltau_q;")
                        rdf = rdf.Define("alltau_decay", "std::vector<Int_t> alltau_decay{decaymode_2}; return alltau_decay;")
                        if process != "data":
                            rdf = rdf.Define("alltau_gen_match", "std::vector<Int_t> alltau_gen_match{gen_match_2}; return alltau_gen_match;")
                        else:
                            rdf = rdf.Define("alltau_gen_match", "-999")
                        rdf = rdf.Define("alltau_mvis", "std::vector<Double_t> alltau_mvis{m_vis}; return alltau_mvis;")
                        rdf = rdf.Define("alltau_mt", "std::vector<Double_t> alltau_mt{mt_1}; return alltau_mt;")
                        rdf = rdf.Define("alltau_mt2", "std::vector<Double_t> alltau_mt2{mt_2}; return alltau_mt2;")
                        rdf = rdf.Define("alltau_svfit", "std::vector<Double_t> alltau_svfit{0.}; return alltau_svfit;")
                        rdf = rdf.Define("alltau_beta", "std::vector<Double_t> alltau_beta{-999.}; return alltau_beta;")

                        rdf = rdf.Define("v4_tau_1", "ROOT::Math::PtEtaPhiMVector v4_tau_1(pt_1,eta_1,phi_1,mass_1); return v4_tau_1;")
                        rdf = rdf.Define("v4_tau_2", "ROOT::Math::PtEtaPhiMVector v4_tau_2(pt_2,eta_2,phi_2,mass_2); return v4_tau_2;")
                        rdf = rdf.Define("v4_met", "ROOT::Math::PtEtaPhiMVector v4_met(met,0.,metphi,0.); return v4_met;")
                        rdf = rdf.Define("alltau_Zpt", "std::vector<Double_t> alltau_Zpt{(v4_tau_1 + v4_tau_2 + v4_met).Pt()}; return alltau_Zpt;")

                        rdf = rdf.Define("alltau_dRToLep", "std::vector<Double_t> alltau_dRToLep{ROOT::VecOps::DeltaR(eta_1, eta_2, phi_1, phi_2)}; return alltau_dRToLep;")
                        rdf = rdf.Define("alltau_dRToOtherLep", "std::vector<Double_t> alltau_dRToOtherLep{1e6}; return alltau_dRToOtherLep;")

                        rdf = rdf.Define("dR_TauToB_1", "ROOT::VecOps::DeltaR(beta_1, eta_2, bphi_1, phi_2)")
                        rdf = rdf.Define("dR_TauToB_2", "ROOT::VecOps::DeltaR(beta_2, eta_2, bphi_2, phi_2)")
                        rdf = rdf.Define("alltau_dRToB", "std::vector<Double_t> alltau_dRToB{min(dR_TauToB_1,dR_TauToB_2)}; return alltau_dRToB;")

                    else:
                        print("Tau collection: Such a channel is not defined: {}".format(config["channel"]))

                
                    print("The current data frame will be saved to {}".format(get_output_name(output_path, process, tau_gen_mode, idx)))
                    rdf.Snapshot(config["tree"], get_output_name(output_path, process, tau_gen_mode, idx), output_feature)
                    print("-" * 50)
                    process_file_list.append(get_output_name(output_path, process, tau_gen_mode, idx))
                
                else:
                    print("No events left after filters. Data frame will not be saved.")
                    print("-" * 50)

            # combining sample files to a single process file, if there are any
            if len(process_file_list) != 0:
                sum_rdf = ROOT.RDataFrame(config["tree"], process_file_list)
                print("The processed files for the {} process are concatenated. The data frame will be saved to {}".format(process, get_output_name(output_path, process, tau_gen_mode)))
                sum_rdf.Snapshot(config["tree"], get_output_name(output_path, process, tau_gen_mode), output_feature)
                print("-" * 50)
            else:
                print("No processed files for the {} process. An empty data frame will be saved to {}".format(process, get_output_name(output_path, process, tau_gen_mode)))
                sum_rdf.Snapshot(config["tree"], get_output_name(output_path, process, tau_gen_mode), [])
                print("-" * 50)

            # delete unneeded sample files after combination
            for rf in process_file_list:
                os.remove(rf)
            process_file_list = []
            print("-" * 50)

            with open(output_path + "/config.yaml", "w") as config_file:
                yaml.dump(config, config_file, default_flow_style=False)
