"""
Main script initializing the fake factor calculation
"""

import argparse
import yaml
import sys
import os

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
    help="Name of a config file in configs/ which contains information for the fake factor calculation step.",
)
parser.add_argument(
    "--config-file",
    default=None,
    help="path to the config file",
)


if __name__ == "__main__":
    args = parser.parse_args()
    # loading of the chosen config file
    if args.config_file is not None:
        with open(args.config_file, "r") as file:
            config = yaml.load(file, yaml.FullLoader)
    else:
        # loading of the chosen config file
        with open("configs/" + args.config + ".yaml", "r") as file:
            config = yaml.load(file, yaml.FullLoader)

    save_path_ffs = "workdir/{}/{}".format(
        config["workdir_name"], config["era"]
    )
    func.check_output_path(os.getcwd() + "/" + save_path_ffs)
    save_path_plots = "workdir/{}/{}/fake_factors/{}".format(
        config["workdir_name"], config["era"], config["channel"]
    )
    func.check_output_path(os.getcwd() + "/" + save_path_plots)

    # start output logging
    sys.stdout = Logger(save_path_plots + "/ff_calculation.log")

    # getting all the input files
    sample_path_list = func.get_samples(config)
    if len(sample_path_list) == 0:
        raise Exception("No input files found!") 

    # check binning of defined categories in the config
    func.check_categories(config)

    # initializing the fake factor calculation
    fake_factors = dict()
    processes = list()

    if "target_process" in config:
        for process in config["target_process"]:
            if process == "QCD":
                print("Calculating fake factors for the QCD process.")
                print("-" * 50)
                cs_expressions = FF_QCD.calculation_QCD_FFs(
                    config, sample_path_list, save_path_plots
                )
                processes.append(process)
            elif process == "Wjets":
                print("Calculating fake factors for the Wjets process.")
                print("-" * 50)
                cs_expressions = FF_Wjets.calculation_Wjets_FFs(
                    config, sample_path_list, save_path_plots
                )
                processes.append(process)
            elif process == "ttbar":
                print("Calculating fake factors for the ttbar process.")
                print("-" * 50)
                cs_expressions = FF_ttbar.calculation_ttbar_FFs(
                    config, sample_path_list, save_path_plots
                )
                processes.append(process)
            else:
                sys.exit(
                    "Target process: Such a process is not defined: {}".format(process)
                )
            fake_factors[process] = cs_expressions

    if "process_fractions" in config:
        print("Calculating the process fractions for the FF application.")
        print("-" * 50)
        fractions = fraction_calculation(config, sample_path_list, save_path_plots)
        fractions = func.get_yields_from_hists(
            fractions, config["process_fractions"]["processes"]
        )

    if config["generate_json"]:
        cs.generate_ff_cs_json(
            config,
            save_path_ffs,
            fake_factors,
            fractions=fractions,
            processes=processes,
            for_corrections=False,
        )

    # dumping config to output directory for documentation
    with open(save_path_plots + "/config.yaml", "w") as config_file:
        yaml.dump(config, config_file, default_flow_style=False)
