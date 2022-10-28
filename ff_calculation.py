import ROOT
import argparse
import yaml
import glob
import sys

from FF_calculation.FF_QCD import calculation_QCD_FFs
from FF_calculation.FF_Wjets import calculation_Wjets_FFs
from FF_calculation.FF_ttbar import calculation_ttbar_FFs

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
        if process == "QCD":
            calculation_QCD_FFs(config, sample_path_list)
        elif process == "Wjets":
            calculation_Wjets_FFs(config, sample_path_list)
        elif process == "ttbar":
            calculation_ttbar_FFs(config, sample_path_list)
        else:
            sys.exit("Target process: Such a process is not defined: {}".format(process))