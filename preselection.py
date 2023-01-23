"""
Script for preprocessing n-tuples for the fake factor calculation
"""

import os
import sys
import argparse
from io import StringIO
from wurlitzer import pipes, STDOUT
import yaml
import ROOT

from helper.Logger import Logger
import helper.filters as filters
import helper.weights as weights
import helper.functions as func


parser = argparse.ArgumentParser()

parser.add_argument(
    "--config",
    default="config",
    help="Name of a config file in configs/ which contains information for the preselection step.",
)

output_feature = [
    "weight",
    "btag_weight",
    "njets",
    "nbtag",
    "id_tau_vsJet_Tight_2",
    "id_tau_vsJet_Medium_2",
    "id_tau_vsJet_VVVLoose_2",
    "id_wgt_tau_vsJet_Tight_2",
    "id_wgt_tau_vsJet_Medium_2",
    "id_wgt_tau_vsJet_VVVLoose_2",
    "q_1",
    "pt_2",
    "phi_2",
    "q_2",
    "gen_match_2",
    "m_vis",
    "mt_1",
    "no_extra_lep",
    "deltaR_ditaupair",
    # following variables are not directly needed for the FF calculation
    "pt_1",
    "eta_1",
    "phi_1",
    "iso_1",
    "eta_2",
    "iso_2",
    "bpt_1",
    "beta_1",
    "bphi_1",
    "bpt_2",
    "mjj",
    "met",
    "metphi",
]


if __name__ == "__main__":
    args = parser.parse_args()

    # loading of the chosen config file
    with open("configs/" + args.config + ".yaml", "r") as file:
        config = yaml.load(file, yaml.FullLoader)
    # loading general dataset info file for xsec and event number
    with open("configs/datasets.yaml", "r") as file:
        datasets = yaml.load(file, yaml.FullLoader)

    # define output path for the preselected samples
    output_path = (
        config["output_path"]
        + "/preselection/"
        + config["era"]
        + "/"
        + config["channel"]
    )
    func.check_output_path(output_path)
    
    # start output logging
    sys.stdout = Logger(output_path + "/preselection.log")

    # these variables are no defined in et, mt ntuples
    if config["channel"] == "tt":
        output_feature.append("id_tau_vsJet_Tight_1")
        output_feature.append("id_tau_vsJet_Medium_1")
        output_feature.append("id_tau_vsJet_VVVLoose_1")
        output_feature.append("id_wgt_tau_vsJet_Tight_1")
        output_feature.append("id_wgt_tau_vsJet_Medium_1")
        output_feature.append("id_wgt_tau_vsJet_VVVLoose_1")

    # going through all wanted processes
    for process in config["processes"]:

        # bookkeeping of samples files due to splitting based on the tau origin (genuine, jet fake, lepton fake)
        process_file_dict = dict()
        for tau_gen_mode in config["processes"][process]["tau_gen_modes"]:
            process_file_dict[tau_gen_mode] = list()

        # going through all contributing samples for the process
        for idx, sample in enumerate(config["processes"][process]["samples"]):
            rdf = ROOT.RDataFrame(config["tree"], func.get_ntuples(config, sample))

            # apply wanted event filters
            selection_conf = config["event_selection"]
            if "had_tau_pt" in selection_conf:
                rdf = filters.had_tau_pt_cut(rdf, config["channel"], selection_conf)
            if "had_tau_eta" in selection_conf:
                rdf = filters.had_tau_eta_cut(rdf, config["channel"], selection_conf)
            if "had_tau_decay_mode" in selection_conf:
                rdf = filters.had_tau_decay_mode_cut(
                    rdf, config["channel"], selection_conf
                )
            if "had_tau_id_vs_ele" in selection_conf:
                rdf = filters.had_tau_id_vsEle_cut(
                    rdf, config["channel"], selection_conf
                )
            if "had_tau_id_vs_mu" in selection_conf:
                rdf = filters.had_tau_id_vsMu_cut(
                    rdf, config["channel"], selection_conf
                )
            if "lep_iso" in selection_conf:
                rdf = filters.lep_iso_cut(rdf, config["channel"], selection_conf)
            if "trigger" in selection_conf:
                rdf = filters.trigger_cut(rdf, config["channel"])

            if process == "embedding":
                rdf = filters.emb_tau_gen_match(rdf, config["channel"])

            # calculate event weights for plotting
            mc_weight_conf = config["mc_weights"]
            rdf = rdf.Define("weight", "1.")
            if process not in ["data", "embedding"]:
                if "generator" in mc_weight_conf:
                    rdf = weights.gen_weight(rdf, datasets[sample])
                if "lumi" in mc_weight_conf:
                    rdf = weights.lumi_weight(rdf, config["era"])
                if "pileup" in mc_weight_conf:
                    rdf = weights.pileup_weight(rdf)
                if "lep_iso" in mc_weight_conf:
                    rdf = weights.lep_iso_weight(rdf, config["channel"])
                if "lep_id" in mc_weight_conf:
                    rdf = weights.lep_id_weight(rdf, config["channel"])
                if "had_tau_id_vs_mu" in mc_weight_conf:
                    rdf = weights.tau_id_vsMu_weight(
                        rdf, config["channel"], selection_conf
                    )
                if "had_tau_id_vs_ele" in mc_weight_conf:
                    rdf = weights.tau_id_vsEle_weight(
                        rdf, config["channel"], selection_conf
                    )
                if "trigger" in mc_weight_conf:
                    rdf = weights.trigger_weight(rdf, config["channel"], process)
                if "Z_pt_reweight" in mc_weight_conf:
                    rdf = weights.Z_pt_reweight(rdf, process)
                if "Top_pt_reweight" in mc_weight_conf:
                    rdf = weights.Top_pt_reweight(rdf, process)

            emb_weight_conf = config["emb_weights"]
            if process == "embedding":
                if "generator" in emb_weight_conf:
                    rdf = weights.emb_gen_weight(rdf)
                if "lep_iso" in emb_weight_conf:
                    rdf = weights.lep_iso_weight(rdf, config["channel"])
                if "lep_id" in emb_weight_conf:
                    rdf = weights.lep_id_weight(rdf, config["channel"])
                if "trigger" in emb_weight_conf:
                    rdf = weights.trigger_weight(rdf, config["channel"], process)

            # default values for some output variables which are not defined in data, embedding; will not be used in FF calculation
            if process == "data":
                rdf = rdf.Define("btag_weight", "1.")
                rdf = rdf.Define("id_wgt_tau_vsJet_Tight_2", "1.")
                rdf = rdf.Define("id_wgt_tau_vsJet_Medium_2", "1.")
                rdf = rdf.Define("id_wgt_tau_vsJet_VVVLoose_2", "1.")
                rdf = rdf.Define("gen_match_2", "-1.")
                if config["channel"] == "tt":
                    rdf = rdf.Define("id_wgt_tau_vsJet_Tight_1", "1.")
                    rdf = rdf.Define("id_wgt_tau_vsJet_Medium_1", "1.")
                    rdf = rdf.Define("id_wgt_tau_vsJet_VVVLoose_1", "1.")
            if process == "embedding":
                rdf = rdf.Define("btag_weight", "1.")
                # rdf = rdf.Define("id_wgt_tau_vsJet_Medium_2", "1.")
                # rdf = rdf.Define("id_wgt_tau_vsJet_VVVLoose_2", "1.")

            # calculate additional variables
            if config["channel"] == "et":
                rdf = rdf.Define(
                    "no_extra_lep",
                    "(extramuon_veto < 0.5) && (extraelec_veto < 0.5) && (dimuon_veto < 0.5)",
                )
            elif config["channel"] == "mt":
                rdf = rdf.Define(
                    "no_extra_lep",
                    "(extramuon_veto < 0.5) && (extraelec_veto < 0.5) && (dimuon_veto < 0.5)",
                )
            elif config["channel"] == "tt":
                rdf = rdf.Define(
                    "no_extra_lep",
                    "(extramuon_veto < 0.5) && (extraelec_veto < 0.5) && (dimuon_veto < 0.5)",
                )
            else:
                print(
                    "Extra lepton veto: Such a channel is not defined: {}".format(
                        config["channel"]
                    )
                )

            # splitting data frame based on the tau origin (genuine, jet fake, lepton fake)
            for tau_gen_mode in config["processes"][process]["tau_gen_modes"]:
                tmp_rdf = rdf
                if tau_gen_mode != "all":
                    tmp_rdf = filters.tau_origin_split(
                        tmp_rdf, config["channel"], tau_gen_mode
                    )

                # redirecting C++ stdout for Report() to python stdout
                out = StringIO()
                with pipes(stdout=out, stderr=STDOUT):
                    tmp_rdf.Report().Print()
                print(out.getvalue())
                print("-" * 50)

                tmp_file_name = func.get_output_name(
                    output_path, process, tau_gen_mode, idx
                )
                # check for empty data frame -> only save/calculate if event number is not zero
                if tmp_rdf.Count().GetValue() != 0:
                    print(
                        "The current data frame will be saved to {}".format(
                            tmp_file_name
                        )
                    )
                    tmp_rdf.Snapshot(config["tree"], tmp_file_name, output_feature)
                    print("-" * 50)
                    process_file_dict[tau_gen_mode].append(tmp_file_name)
                else:
                    print("No events left after filters. Data frame will not be saved.")
                    print("-" * 50)

        # combining all files of a process and tau origin
        for tau_gen_mode in config["processes"][process]["tau_gen_modes"]:
            out_file_name = func.get_output_name(output_path, process, tau_gen_mode)
            # combining sample files to a single process file, if there are any
            if len(process_file_dict[tau_gen_mode]) != 0:
                sum_rdf = ROOT.RDataFrame(
                    config["tree"], process_file_dict[tau_gen_mode]
                )
                print(
                    "The processed files for the {} process are concatenated. The data frame will be saved to {}".format(
                        process, out_file_name
                    )
                )
                sum_rdf.Snapshot(config["tree"], out_file_name, output_feature)
                print("-" * 50)
            else:
                print(
                    "No processed files for the {} process. An empty data frame will be saved to {}".format(
                        process, out_file_name
                    )
                )
                sum_rdf.Snapshot(config["tree"], out_file_name, [])
                print("-" * 50)

            # delete unneeded sample files after combination
            for rf in process_file_dict[tau_gen_mode]:
                os.remove(rf)
            process_file_list = []
            print("-" * 50)

    # dumping config to output directory for documentation
    with open(output_path + "/config.yaml", "w") as config_file:
        yaml.dump(config, config_file, default_flow_style=False)
