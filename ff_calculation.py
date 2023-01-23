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
    help="Name of a config file in configs/ which contains information for the fake factor calculation step.",
)
parser.add_argument(
    "--only-correction",
    action="store_true",
    help="Using this argument means to skip the calculation of the fake factors and directly calculate the corrections.",
)


if __name__ == "__main__":
    args = parser.parse_args()

    # loading of the chosen config file
    with open("configs/" + args.config + ".yaml", "r") as file:
        config = yaml.load(file, yaml.FullLoader)

    save_path = "workdir/{}/{}/".format(config["workdir_name"], config["channel"])
    func.check_output_path(os.getcwd() + "/" + save_path)

    # start output logging
    sys.stdout = Logger(save_path + "ff_calculation.log")

    # getting all the input files
    sample_path_list = func.get_samples(config)

    if not args.only_correction:
        # check binning of defined categories in the config
        func.check_categories(config)

        # initializing the fake factor calculation
        fake_factors = dict()
        if "target_process" in config:
            for process in config["target_process"]:
                if process == "QCD":
                    print("Calculating fake factors for the QCD process.")
                    print("-" * 50)
                    cs_expressions = FF_QCD.calculation_QCD_FFs(
                        config, sample_path_list
                    )
                elif process == "Wjets":
                    print("Calculating fake factors for the Wjets process.")
                    print("-" * 50)
                    cs_expressions = FF_Wjets.calculation_Wjets_FFs(
                        config, sample_path_list
                    )
                elif process == "ttbar":
                    print("Calculating fake factors for the ttbar process.")
                    print("-" * 50)
                    cs_expressions = FF_ttbar.calculation_ttbar_FFs(
                        config, sample_path_list
                    )
                else:
                    sys.exit(
                        "Target process: Such a process is not defined: {}".format(
                            process
                        )
                    )
                fake_factors[process] = cs_expressions

        if "process_fractions" in config:
            print("Calculating the process fractions for the FF application.")
            print("-" * 50)
            fractions = fraction_calculation(config, sample_path_list)
            fractions = func.get_yields_from_hists(
                fractions, config["process_fractions"]["processes"]
            )

        if config["generate_json"]:
            cs.generate_ff_cs_json(config, fake_factors, fractions, save_path)

    if config["do_corrections"]:
        # This activates implicit multi-threading
        ROOT.EnableImplicitMT()

        corrections = dict()

        if "target_process" in config:
            for process in config["target_process"]:
                if process == "QCD":
                    print("Calculating corrections for the QCD process.")
                    print("-" * 50)
                    corrections[process] = dict()
                    corr = FF_QCD.non_closure_correction(config, sample_path_list)
                    corrections[process]["non_closure"] = corr
                elif process == "Wjets":
                    print("Calculating corrections for the Wjets process.")
                    print("-" * 50)
                    corrections[process] = dict()
                    corr = FF_Wjets.non_closure_correction(config, sample_path_list)
                    corrections[process]["non_closure"] = corr
                elif process == "ttbar":
                    print("Calculating corrections for the ttbar process.")
                    print("-" * 50)
                    corrections[process] = dict()
                    corr = FF_ttbar.non_closure_correction(config, sample_path_list)
                    corrections[process]["non_closure"] = corr
                else:
                    sys.exit(
                        "Target process: Such a process is not defined: {}".format(
                            process
                        )
                    )

        cs.generate_corr_cs_json(config, corrections, save_path)
