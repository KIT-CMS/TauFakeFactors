"""
Main script initializing the fake factor calculation
"""

import argparse
import yaml
import sys
import os
import ROOT

import FF_calculation.FF_QCD as FF_QCD
import FF_calculation.FF_Wjets as FF_Wjets
import FF_calculation.FF_ttbar as FF_ttbar
from FF_calculation.fractions import fraction_calculation
from helper.Logger import Logger
import helper.correctionlib_json as cs
import helper.functions as func

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config",
    default="hist_config",
    help="Name of a config file in configs/ which contains information for the fake factor corrections step.",
)


if __name__ == "__main__":
    args = parser.parse_args()

    # loading of the chosen config file
    with open("configs/" + args.config + ".yaml", "r") as file:
        corr_config = yaml.load(file, yaml.FullLoader)
    with open("workdir/" + corr_config["workdir_name"] + "/" + corr_config["era"] + "/fake_factors/" + corr_config["channel"] + "/config.yaml", "r") as file:
        config = yaml.load(file, yaml.FullLoader)

    save_path_ffs = "workdir/{}/{}/corrections".format(config["workdir_name"], config["era"])
    func.check_output_path(os.getcwd() + "/" + save_path_ffs)
    save_path_plots = "workdir/{}/{}/corrections/{}".format(config["workdir_name"], config["era"], config["channel"])
    func.check_output_path(os.getcwd() + "/" + save_path_plots)

    # start output logging
    sys.stdout = Logger(save_path_plots + "/ff_corrections.log")

    # getting all the input files
    sample_path_list = func.get_samples(config)

    # initializing the fake factor calculation
    fake_factors = dict()
    processes = list()

    if "target_process" in corr_config:
        for process in corr_config["target_process"]:
            if process == "QCD" and "SR_DR" in corr_config["target_process"]["QCD"]:
                print("Calculating fake factors for the SR-DR correction for the QCD process.")
                print("-" * 50)
                temp_conf = config.copy()
                temp_conf["target_process"]["QCD"]["SRlike_cuts"]["lep_iso"] = ">=0.15"
                temp_conf["target_process"]["QCD"]["ARlike_cuts"]["lep_iso"] = ">=0.15"
                cs_expressions = FF_QCD.calculation_QCD_FFs(
                    temp_conf, sample_path_list, save_path_plots
                )
                fake_factors[process] = cs_expressions
                processes.append(process)

            elif process == "Wjets" and "SR_DR" in corr_config["target_process"]["Wjets"]:
                print("Calculating fake factors for the SR-DR correction for the Wjets process.")
                print("-" * 50)
                temp_conf = config.copy()
                temp_conf["target_process"]["Wjets"]["SRlike_cuts"]["lep_mt"] = "<=70"
                temp_conf["target_process"]["Wjets"]["ARlike_cuts"]["lep_mt"] = "<=70"
                cs_expressions = FF_Wjets.calculation_Wjets_FFs(
                    temp_conf, sample_path_list, save_path_plots
                )
                fake_factors[process] = cs_expressions
                processes.append(process)

    if corr_config["generate_json"]:
        cs.generate_ff_cs_json(config, save_path_plots, fake_factors, fractions=None, processes=processes, for_corrections=True)


    # This activates implicit multi-threading
    ROOT.EnableImplicitMT()

    corrections = dict()

    if "target_process" in corr_config:
        for process in corr_config["target_process"]:
            if process == "QCD" and "non_closure" in corr_config["target_process"]["QCD"]:
                print("Calculating corrections for the QCD process.")
                print("-" * 50)
                corrections[process] = dict()
                corr = FF_QCD.non_closure_correction(config, corr_config, sample_path_list, save_path_plots)
                corrections[process]["non_closure"] = corr
            elif process == "Wjets" and "non_closure" in corr_config["target_process"]["Wjets"]:
                print("Calculating corrections for the Wjets process.")
                print("-" * 50)
                corrections[process] = dict()
                corr = FF_Wjets.non_closure_correction(config, corr_config, sample_path_list, save_path_plots)
                corrections[process]["non_closure"] = corr
            elif process == "ttbar" and "non_closure" in corr_config["target_process"]["ttbar"]:
                print("Calculating corrections for the ttbar process.")
                print("-" * 50)
                corrections[process] = dict()
                corr = FF_ttbar.non_closure_correction(config, corr_config, sample_path_list, save_path_plots)
                corrections[process]["non_closure"] = corr
            else:
                sys.exit(
                    "Target process: Correction: Such a process is not defined: {}".format(
                        process
                    )
                )

    cs.generate_corr_cs_json(corr_config, corrections, save_path_ffs)