"""
Main script initializing the fake factor calculation
"""

import argparse
import yaml
import sys

from FF_calculation.FF_QCD import calculation_QCD_FFs
from FF_calculation.FF_Wjets import calculation_Wjets_FFs
from FF_calculation.FF_ttbar import calculation_ttbar_FFs
from FF_calculation.fractions import fraction_calculation
from helper.correctionlib_json import generate_cs_json
import helper.functions as func

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config",
    default="hist_config",
    help="Name of a config file in configs/ which contains information for the fake factor calculation step.",
)


if __name__ == "__main__":
    args = parser.parse_args()

    # loading of the chosen config file
    with open("configs/" + args.config + ".yaml", "r") as file:
        config = yaml.load(file, yaml.FullLoader)

    # getting all the input files
    sample_path_list = func.get_samples(config)
    
    # initializing the fake factor calculation
    fake_factors = dict()
    if "target_process" in config:
        for process in config["target_process"]:
            if process == "QCD":
                print("Calculating fake factors for the QCD process.")
                print("-" * 50)
                cs_expressions = calculation_QCD_FFs(config, sample_path_list)
            elif process == "Wjets":
                print("Calculating fake factors for the Wjets process.")
                print("-" * 50)
                cs_expressions = calculation_Wjets_FFs(config, sample_path_list)
            elif process == "ttbar":
                print("Calculating fake factors for the ttbar process.")
                print("-" * 50)
                cs_expressions = calculation_ttbar_FFs(config, sample_path_list)
            else:
                sys.exit(
                    "Target process: Such a process is not defined: {}".format(process)
                )
            fake_factors[process] = cs_expressions

    if "process_fractions" in config:
        print("Calculating the process fractions for the FF application.")
        print("-" * 50)
        fractions = fraction_calculation(config, sample_path_list)
        fractions = func.get_yields_from_hists(fractions, config["process_fractions"]["processes"])
    
    if config["generate_json"]:
        generate_cs_json(config, fake_factors, fractions)