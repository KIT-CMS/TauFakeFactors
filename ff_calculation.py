"""
Main script initializing the fake factor calculation
"""

import argparse
import logging
import os
from typing import Dict, List, Tuple, Union

import FF_calculation.FF_QCD as FF_QCD
import FF_calculation.FF_ttbar as FF_ttbar
import FF_calculation.FF_Wjets as FF_Wjets
import helper.correctionlib_json as corrlib
import helper.ff_functions as ff_func
import helper.functions as func
import CustomLogging as logging_helper
from FF_calculation.fractions import fraction_calculation
from helper.hooks_and_patches import Histo1DPatchedRDataFrame, PassThroughWrapper

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

parser.add_argument(
    "--ignore-cached-intermediary-steps",
    action="store_true",
    help="""
        Flag to use intermediary filtered ROOT RDataFrames even if cached versions are available.
    """,
)
parser.add_argument(
    "--log-level",
    default="INFO",
    help="Logging level to use. (default: INFO)",
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

FF_DATA_SCALING_FACTOR_CALCULATION_FUNCTIONS = {
    **{k: lambda *args, **kwargs: (None, None) for k in {
        "QCD",
        "QCD_subleading",
        "Wjets",
        "process_fractions",
        "process_fractions_subleading",
    }
    },  # only necessary for ttbar and ttbar_subleading
    "ttbar": FF_ttbar.calculation_FF_data_scaling_factor,
    "ttbar_subleading": FF_ttbar.calculation_FF_data_scaling_factor,
}


@logging_helper.LogDecorator().grouped_logs()
def FF_calculation(
    config: Dict[str, Union[str, Dict, List]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    logger: str,
) -> Dict[str, Union[Dict[str, str], Dict[str, Dict[str, str]]]]:
    """
    This function calculates fake factors for a given process or fractions of the processes.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        process: Process for which the fake factors are calculated
        logger: Name of the logger that should be used

    Return:
        Dictionary where the categories are defined as keys and and the values are the fitted functions (including variations)
        e.g. corrlib_expressions[CATEGORY_1][CATEGORY_2][VARIATION] if dimension of categories is 2
    """

    is_fraction = "fraction" in process
    split_limit = 1 if is_fraction else 2

    process_conf = config[process] if is_fraction else config["target_processes"][process]

    split_collections = ff_func.SplitQuantities(process_conf)

    assert len(split_collections.split_variables) <= split_limit, f"Category splitting of {process} is only defined up to {split_limit} dimensions."

    try:
        SRlike_hists, ARlike_hists = FF_DATA_SCALING_FACTOR_CALCULATION_FUNCTIONS[process](
            config=config,
            process_conf=process_conf,
            sample_paths=sample_paths,
            process=process,
            logger=logger,
        )
        results = func.optional_process_pool(
            args_list=[
                (
                    split_collection,
                    config,
                    process_conf,
                    process,
                    sample_paths,
                    output_path,
                    logger,
                    SRlike_hists,
                    ARlike_hists,
                )
                for split_collection in split_collections
            ],
            function=FF_CALCULATION_FUNCTIONS[process],
        )
    except KeyError:
        raise Exception(f"Target process: Such a process is not known: {process}")

    if is_fraction:
        fractions = dict()
        for result in results:
            fractions.update(result)

        return ff_func.get_yields_from_hists(
            hists=fractions,
            processes=config[process]["processes"],
        )
    else:
        return ff_func.fill_corrlib_expression(results, split_collections.split_variables)


@logging_helper.LogDecorator().grouped_logs(extractor=lambda args: f"ff_calculation.{args[0]}")
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

    log.info(f"Calculating fake factors for the {process} process.")
    log.info("-" * 50)
    return (
        args,
        FF_calculation(
            config=config,
            sample_paths=sample_paths,
            output_path=output_path,
            process=process,
            logger=f"ff_calculation.{process}",
        ),
    )


if __name__ == "__main__":
    args = parser.parse_args()

    func.RuntimeVariables.USE_MULTIPROCESSING = not args.disable_multiprocessing
    func.RuntimeVariables.USE_CACHED_INTERMEDIATE_STEPS = not args.ignore_cached_intermediary_steps

    # loading of the chosen config file
    config = func.load_config(args.config_file)

    func.RuntimeVariables.INPUT_FILE_PATH = os.path.join(config["output_path"], config["era"], config["channel"])

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

    logging_helper.LOG_LEVEL = getattr(logging, args.log_level.upper(), logging.INFO)
    log = logging_helper.setup_logging(
        output_file=save_path_plots + "/ff_calculation.log",
        logger=logging.getLogger("ff_calculation"),
        level=logging_helper.LOG_LEVEL,
    )

    # getting all the ntuple input files
    sample_paths = func.get_samples(config=config)
    if len(sample_paths) == 0:
        raise Exception("No input files found!")

    # check binning of defined categories in the config
    func.check_categories(config=config)

    func.RuntimeVariables.RDataFrameWrapper = PassThroughWrapper
    if config.get("use_center_of_mass_bins", True):
        func.RuntimeVariables.RDataFrameWrapper = Histo1DPatchedRDataFrame

    # initializing the fake factor calculation
    if "target_processes" in config:
        args_list = [
            (process, config, sample_paths, save_path_plots)
            for process in config["target_processes"]
        ]
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
        func.configured_yaml.dump(config, config_file)

    with open(os.path.join(save_path_plots, "done"), "w") as done_file:
        done_file.write("")

    log.info("Fake factor calculation finished successfully.")
