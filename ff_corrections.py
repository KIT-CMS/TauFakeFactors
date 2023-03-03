"""
Main script initializing the fake factor calculation
"""

import argparse
import yaml
import sys
import os
import copy
import ROOT

import FF_calculation.FF_QCD as FF_QCD
import FF_calculation.FF_Wjets as FF_Wjets
import FF_calculation.FF_ttbar as FF_ttbar
from helper.Logger import Logger
import helper.correctionlib_json as cs
import helper.functions as func

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config",
    default="hist_config",
    help="Name of a config file in configs/ which contains information for the fake factor corrections step.",
)
parser.add_argument(
    "--only-corrections",
    action="store_true",
    help="Using this argument means to skip the calculation of the fake factors for additional regions and directly calculate the corrections. This is useful if you already calculated the needed fake factors.",
)


if __name__ == "__main__":
    args = parser.parse_args()

    # loading of the chosen config file
    with open("configs/" + args.config + ".yaml", "r") as file:
        corr_config = yaml.load(file, yaml.FullLoader)
    with open(
        "workdir/"
        + corr_config["workdir_name"]
        + "/"
        + corr_config["era"]
        + "/fake_factors/"
        + corr_config["channel"]
        + "/config.yaml",
        "r",
    ) as file:
        config = yaml.load(file, yaml.FullLoader)

    save_path_ffs = "workdir/{}/{}".format(
        config["workdir_name"], config["era"]
    )
    func.check_output_path(os.getcwd() + "/" + save_path_ffs)
    save_path_plots = "workdir/{}/{}/corrections/{}".format(
        config["workdir_name"], config["era"], config["channel"]
    )
    func.check_output_path(os.getcwd() + "/" + save_path_plots)

    # start output logging
    sys.stdout = Logger(save_path_plots + "/ff_corrections.log")

    # getting all the input files
    sample_path_list = func.get_samples(config)

    # initializing the fake factor calculation
    fake_factors = dict()
    processes = list()

    if not args.only_corrections:
        if "target_process" in corr_config:
            for process in corr_config["target_process"]:
                if process == "QCD" and "DR_SR" in corr_config["target_process"]["QCD"]:
                    print(
                        "Calculating fake factors for the SR-DR correction for the QCD process."
                    )
                    print("-" * 50)
                    temp_conf = copy.deepcopy(config)
                    func.modify_config(
                        temp_conf,
                        process,
                        corr_config["target_process"]["QCD"]["DR_SR"],
                    )
                    cs_expressions = FF_QCD.calculation_QCD_FFs(
                        temp_conf, sample_path_list, save_path_plots
                    )
                    fake_factors[process] = cs_expressions
                    processes.append(process)

                elif (
                    process == "Wjets"
                    and "DR_SR" in corr_config["target_process"]["Wjets"]
                ):
                    print(
                        "Calculating fake factors for the SR-DR correction for the Wjets process."
                    )
                    print("-" * 50)
                    temp_conf = copy.deepcopy(config)
                    func.modify_config(
                        temp_conf,
                        process,
                        corr_config["target_process"]["Wjets"]["DR_SR"],
                    )
                    cs_expressions = FF_Wjets.calculation_Wjets_FFs(
                        temp_conf, sample_path_list, save_path_plots
                    )
                    fake_factors[process] = cs_expressions
                    processes.append(process)

        if corr_config["generate_json"]:
            cs.generate_ff_cs_json(
                config,
                save_path_plots,
                fake_factors,
                fractions=None,
                processes=processes,
                for_corrections=True,
            )

    # This activates implicit multi-threading
    ROOT.EnableImplicitMT()

    corrections = {"QCD": dict(), "Wjets": dict(), "ttbar": dict()}

    if "target_process" in corr_config:
        for process in corr_config["target_process"]:
            if process == "QCD":
                if "non_closure" in corr_config["target_process"]["QCD"]:
                    for closure_corr in corr_config["target_process"]["QCD"][
                        "non_closure"
                    ]:
                        print(
                            "Calculating closure correction for the QCD process dependent on {}.".format(
                                closure_corr
                            )
                        )
                        print("-" * 50)
                        temp_conf = copy.deepcopy(config)
                        func.modify_config(
                            temp_conf,
                            process,
                            corr_config["target_process"]["QCD"]["non_closure"][
                                closure_corr
                            ],
                        )

                        corr = FF_QCD.non_closure_correction(
                            temp_conf,
                            corr_config,
                            closure_corr,
                            sample_path_list,
                            save_path_plots,
                        )
                        corrections[process]["non_closure_" + closure_corr] = corr

                if "DR_SR" in corr_config["target_process"]["QCD"]:
                    print("Calculating DR to SR correction for the QCD process.")
                    print("-" * 50)
                    temp_conf = copy.deepcopy(config)
                    func.modify_config(
                        temp_conf,
                        process,
                        corr_config["target_process"]["QCD"]["DR_SR"],
                    )

                    corr = FF_QCD.DR_SR_correction(
                        temp_conf, corr_config, sample_path_list, save_path_plots
                    )
                    corrections[process]["DR_SR"] = corr

            elif process == "Wjets":
                if "non_closure" in corr_config["target_process"]["Wjets"]:
                    for closure_corr in corr_config["target_process"]["Wjets"][
                        "non_closure"
                    ]:
                        print(
                            "Calculating closure correction for the Wjets process dependent on {}.".format(
                                closure_corr
                            )
                        )
                        print("-" * 50)
                        temp_conf = copy.deepcopy(config)
                        func.modify_config(
                            temp_conf,
                            process,
                            corr_config["target_process"]["Wjets"]["non_closure"][
                                closure_corr
                            ],
                        )

                        corr = FF_Wjets.non_closure_correction(
                            temp_conf,
                            corr_config,
                            closure_corr,
                            sample_path_list,
                            save_path_plots,
                        )
                        corrections[process]["non_closure_" + closure_corr] = corr

                if "DR_SR" in corr_config["target_process"]["Wjets"]:
                    print("Calculating DR to SR correction for the Wjets process.")
                    print("-" * 50)
                    temp_conf = copy.deepcopy(config)
                    func.modify_config(
                        temp_conf,
                        process,
                        corr_config["target_process"]["Wjets"]["DR_SR"],
                    )

                    corr = FF_Wjets.DR_SR_correction(
                        temp_conf, corr_config, sample_path_list, save_path_plots
                    )
                    corrections[process]["DR_SR"] = corr

            elif (
                process == "ttbar"
                and "non_closure" in corr_config["target_process"]["ttbar"]
            ):
                for closure_corr in corr_config["target_process"]["ttbar"][
                    "non_closure"
                ]:
                    print(
                        "Calculating closure correction for the ttbar process dependent on {}.".format(
                            closure_corr
                        )
                    )
                    print("-" * 50)
                    temp_conf = copy.deepcopy(config)
                    func.modify_config(
                        temp_conf,
                        process,
                        corr_config["target_process"]["ttbar"]["non_closure"][
                            closure_corr
                        ],
                    )

                    corr = FF_ttbar.non_closure_correction(
                        temp_conf,
                        corr_config,
                        closure_corr,
                        sample_path_list,
                        save_path_plots,
                    )
                    corrections[process]["non_closure_" + closure_corr] = corr

            else:
                sys.exit(
                    "Target process: Correction: Such a process is not defined: {}".format(
                        process
                    )
                )

    cs.generate_corr_cs_json(corr_config, corrections, save_path_ffs)
