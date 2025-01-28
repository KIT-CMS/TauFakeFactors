"""
Main script initializing the fake factor calculation
"""

import argparse
import yaml
import os
import logging
import concurrent.futures
from typing import Dict, Union, List, Tuple

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


def run_ff_calculation(
    args: Tuple[str, Dict[str, Union[Dict, List, str]], List[str], str]
) -> Dict:
    """
    This function can be used for multiprocessing. It runs the fake factor calculation step for a specified process.

    Args:
        args: Tuple with a process name, a configuration for this process, a list of all sample paths and a path for the output

    Return:
        Depending on the "process" either a dictionary with fake factor function expressions or a dictionary with process fraction values
    """
    process, config, sample_paths, output_path = args
    log = logging.getLogger(f"ff_calculation.{process}")

    if process in ["QCD", "QCD_subleading"]:
        log.info(f"Calculating fake factors for the {process} process.")
        log.info("-" * 50)
        result = FF_QCD.calculation_QCD_FFs(
            config=config,
            sample_paths=sample_paths,
            output_path=output_path,
            process=process,
            logger=f"ff_calculation.{process}",
        )
    elif process == "Wjets":
        log.info("Calculating fake factors for the Wjets process.")
        log.info("-" * 50)
        result = FF_Wjets.calculation_Wjets_FFs(
            config=config,
            sample_paths=sample_paths,
            output_path=output_path,
            logger=f"ff_calculation.{process}",
        )
    elif process in ["ttbar", "ttbar_subleading"]:
        log.info(f"Calculating fake factors for the {process} process.")
        log.info("-" * 50)
        result = FF_ttbar.calculation_ttbar_FFs(
            config=config,
            sample_paths=sample_paths,
            output_path=output_path,
            process=process,
            logger=f"ff_calculation.{process}",
        )
    elif process in ["process_fractions", "process_fractions_subleading"]:
        log.info(f"Calculating the process {process} for the FF application.")
        log.info("-" * 50)
        result = fraction_calculation(
            config=config,
            sample_paths=sample_paths,
            output_path=output_path,
            process=process,
            logger=f"ff_calculation.{process}",
        )
    else:
        raise Exception(f"Target process: Such a process is not known: {process}")

    return result


if __name__ == "__main__":
    args = parser.parse_args()

    # loading of the chosen config file
    config = func.load_config(args.config_file)


    save_path_ffs = os.path.join("workdir", config["workdir_name"], config["era"])
    func.check_path(path=os.path.join(os.getcwd(), save_path_ffs))
    save_path_plots = os.path.join(
        "workdir",
        config["workdir_name"],
        config["era"],
        "fake_factors",
        config["channel"],
    )
    func.check_path(path=os.path.join(os.getcwd(), save_path_plots))

    # start output logging
    subcategories = list(config["target_processes"].keys())
    if "process_fractions" in config:
        subcategories = subcategories + ["process_fractions"]
    if "process_fractions_subleading" in config:
        subcategories = subcategories + ["process_fractions_subleading"]

    func.setup_logger(
        log_file=save_path_plots + "/ff_calculation.log",
        log_name="ff_calculation",
        subcategories=subcategories,
    )

    # getting all the input files
    sample_paths = func.get_samples(config=config)
    if len(sample_paths) == 0:
        raise Exception("No input files found!")

    # check binning of defined categories in the config
    func.check_categories(config=config)

    # initializing the fake factor calculation
    fake_factors = dict()
    fractions = None
    fractions_subleading = None

    if "target_processes" in config:
        args_list = [
            (process, config, sample_paths, save_path_plots)
            for process in config["target_processes"]
        ]
        if "process_fractions" in config:
            args_list.append(("process_fractions", config, sample_paths, save_path_plots))
        if "process_fractions_subleading" in config:
            args_list.append(("process_fractions_subleading", config, sample_paths, save_path_plots))

        with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
            for args, result in zip(
                args_list, executor.map(run_ff_calculation, args_list)
            ):
                if args[0] in config["target_processes"]:
                    fake_factors[args[0]] = result
                elif args[0] == "process_fractions":
                    fractions = result
                elif args[0] == "process_fractions_subleading":
                    fractions_subleading = result
    else:
        raise Exception("No target processes are defined!")

    corrlib.generate_ff_corrlib_json(
        config=config,
        ff_functions=fake_factors,
        fractions=fractions,
        fractions_subleading=fractions_subleading,
        output_path=save_path_ffs,
        for_corrections=False,
    )

    # dumping config to output directory for documentation
    with open(save_path_plots + "/config.yaml", "w") as config_file:
        yaml.dump(config, config_file, default_flow_style=False)
