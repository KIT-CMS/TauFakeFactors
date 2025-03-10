"""
Main script initializing the fake factor calculation
"""

import argparse
import yaml
import os
import logging
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

parser.add_argument(
    "--disable-multiprocessing",
    action="store_true",
    help="Flag to disable multiprocessing for debugging purposes.",
)

FF_CALCULATION_FUNCTIONS = {
    "QCD": FF_QCD.calculation_QCD_FFs,
    "QCD_subleading": FF_QCD.calculation_QCD_FFs,
    "Wjets": FF_Wjets.calculation_Wjets_FFs,
    "ttbar": FF_ttbar.calculation_ttbar_FFs,
    "ttbar_subleading": FF_ttbar.calculation_ttbar_FFs,
    "process_fractions": fraction_calculation,
    "process_fractions_subleading": fraction_calculation,
}


def run_ff_calculation(
    args: Tuple[str, Dict[str, Union[Dict, List, str]], List[str], str]
) -> Tuple[Tuple, Dict]:
    """
    This function can be used for multiprocessing. It runs the fake factor calculation step for a specified process.

    Args:
        args: Tuple with a process name, a configuration for this process, a list of all sample paths and a path for the output

    Return:
        Depending on the "process" either a dictionary with fake factor function expressions or a dictionary with process fraction values
    """
    process, config, sample_paths, output_path = args
    log = logging.getLogger(f"ff_calculation.{process}")

    try:
        log.info(f"Calculating fake factors for the {process} process.")
        log.info("-" * 50)
        return (
            args,
            FF_CALCULATION_FUNCTIONS[process](
                config=config,
                sample_paths=sample_paths,
                output_path=output_path,
                process=process,
                logger=f"ff_calculation.{process}",
            ),
        )
    except KeyError:
        raise Exception(f"Target process: Such a process is not known: {process}")


if __name__ == "__main__":
    args = parser.parse_args()

    func.RuntimeVariables.USE_MULTIPROCESSING = not args.disable_multiprocessing

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

    # getting all the ntuple input files
    sample_paths = func.get_samples(config=config)
    if len(sample_paths) == 0:
        raise Exception("No input files found!")

    # check binning of defined categories in the config
    func.check_categories(config=config)

    # initializing the fake factor calculation
    args_list = []
    if "target_processes" in config:
        args_list.extend(
            [
                (process, config, sample_paths, save_path_plots)
                for process in config["target_processes"]
            ]
        )
    if "process_fractions" in config:
        args_list.append(
            ("process_fractions", config, sample_paths, save_path_plots)
        )
    if "process_fractions_subleading" in config:
        args_list.append(
            ("process_fractions_subleading", config, sample_paths, save_path_plots)
        )

    results = func.optional_process_pool(
        args_list=args_list,
        function=run_ff_calculation,
    )

    fake_factors, fractions, fractions_subleading = {}, None, None
    for args, result in results:
        if "target_processes" in config and args[0] in config["target_processes"]:
            fake_factors[args[0]] = result
        elif args[0] == "process_fractions":
            fractions = result
        elif args[0] == "process_fractions_subleading":
            fractions_subleading = result

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

    with open(os.path.join(save_path_plots, "done"), "w") as done_file:
        done_file.write("")
