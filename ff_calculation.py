"""
Main script initializing the fake factor calculation
"""

import argparse
import yaml
import os
import logging

import FF_calculation.FF_QCD as FF_QCD
import FF_calculation.FF_Wjets as FF_Wjets
import FF_calculation.FF_ttbar as FF_ttbar
from FF_calculation.fractions import fraction_calculation
import helper.correctionlib_json as corrlib
import helper.functions as func

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config-file",
    default=None,
    help="Path to the config file which contains information for the fake factor calculation step.",
)


if __name__ == "__main__":
    args = parser.parse_args()

    # loading of the chosen config file
    with open(args.config_file, "r") as file:
        config = yaml.load(file, yaml.FullLoader)

    save_path_ffs = os.path.join("workdir", config["workdir_name"], config["era"]) 
    func.check_path(path=os.path.join(os.getcwd(), save_path_ffs))
    save_path_plots = os.path.join("workdir", config["workdir_name"], config["era"], "fake_factors", config["channel"])  
    func.check_path(path=os.path.join(os.getcwd(), save_path_plots))

    # start output logging
    func.setup_logger(log_file=save_path_plots+"/ff_calculation.log", log_name="ff_calculation")
    log = logging.getLogger("ff_calculation")

    # getting all the input files
    sample_paths = func.get_samples(config=config)
    if len(sample_paths) == 0:
        raise Exception("No input files found!")

    # check binning of defined categories in the config
    func.check_categories(config=config)

    # initializing the fake factor calculation
    fake_factors = dict()

    if "target_processes" in config:
        for process in config["target_processes"]:            
            if process in ["QCD", "QCD_subleading"]:
                log.info(f"Calculating fake factors for the {process} process.")
                log.info("-" * 50)
                corrlib_expressions = FF_QCD.calculation_QCD_FFs(
                    config=config, 
                    sample_paths=sample_paths, 
                    output_path=save_path_plots, 
                    process=process,
                )
            elif process == "Wjets":
                log.info("Calculating fake factors for the Wjets process.")
                log.info("-" * 50)
                corrlib_expressions = FF_Wjets.calculation_Wjets_FFs(
                    config=config, 
                    sample_paths=sample_paths, 
                    output_path=save_path_plots, 
                )
            elif process == "ttbar":
                log.info("Calculating fake factors for the ttbar process.")
                log.info("-" * 50)
                corrlib_expressions = FF_ttbar.calculation_ttbar_FFs(
                    config=config, 
                    sample_paths=sample_paths, 
                    output_path=save_path_plots, 
                )
            else:
                raise Exception(
                    f"Target process: Such a process is not known: {process}"
                )
            fake_factors[process] = corrlib_expressions

    else:
        raise Exception("No target processes are defined!")

    if "process_fractions" in config:
        log.info("Calculating the process fractions for the FF application.")
        log.info("-" * 50)
        fractions = fraction_calculation(
            config=config, 
            sample_paths=sample_paths, 
            output_path=save_path_plots,
        )
    else:
        raise Exception("The process fraction calculation is needed but not defined.")

    if config["generate_json"]:
        corrlib.generate_ff_corrlib_json(
            config=config, 
            ff_functions=fake_factors,
            fractions=fractions,
            output_path=save_path_ffs,
            for_corrections=False,
        )

    # dumping config to output directory for documentation
    with open(save_path_plots + "/config.yaml", "w") as config_file:
        yaml.dump(config, config_file, default_flow_style=False)
