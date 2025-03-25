"""
Main script initializing the fake factor correction calculation
"""

import argparse
import copy
import logging
import os
from copy import deepcopy
from typing import Any, Dict, List, Tuple, Union

import numpy as np
import yaml

import FF_calculation.FF_QCD as FF_QCD
import FF_calculation.FF_ttbar as FF_ttbar
import FF_calculation.FF_Wjets as FF_Wjets
import helper.correctionlib_json as corrlib
import helper.ff_functions as ff_func
import helper.functions as func
from ff_calculation import FF_calculation
from helper.ff_evaluators import FakeFactorCorrectionEvaluator, FakeFactorEvaluator
from helper.hooks_and_patches import Histo1DPatchedRDataFrame, PassThroughWrapper

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config-file",
    default=None,
    help="Path to the config file which contains information for the fake factor corrections step.",
)
parser.add_argument(
    "--skip-DRtoSR-ffs",
    action="store_true",
    help="Using this argument means to skip the calculation of the fake factors for the DR to SR correction and directly calculate the corrections. This is useful if you already calculated the needed additional fake factors.",
)
parser.add_argument(
    "--only-main-corrections",
    action="store_true",
    help="Using this argument means to skip the calculation of the fake factors and corrections for the DR to SR correction and directly calculate the main corrections. This is useful if you already calculated the needed additional fake factors.",
)

parser.add_argument(
    "--disable-multiprocessing",
    action="store_true",
    help="Flag to disable multiprocessing for debugging purposes.",
)

NON_CLOSURE_CORRECTION_FUNCTIONS = {
    "QCD": FF_QCD.non_closure_correction,
    "QCD_subleading": FF_QCD.non_closure_correction,
    "Wjets": FF_Wjets.non_closure_correction,
    "ttbar": FF_ttbar.non_closure_correction,
    "ttbar_subleading": FF_ttbar.non_closure_correction,
}

DR_SR_CORRECTION_FUNCTIONS = {
    "QCD": FF_QCD.DR_SR_correction,
    "QCD_subleading": FF_QCD.DR_SR_correction,
    "Wjets": FF_Wjets.DR_SR_correction,
}


def non_closure_correction(
    config: Dict[str, Union[str, Dict, List]],
    corr_config: Dict[str, Union[str, Dict]],
    sample_paths: List[str],
    output_path: str,
    process: str,
    closure_variable: str,
    evaluator: FakeFactorEvaluator,
    corr_evaluators: List[FakeFactorCorrectionEvaluator],
    for_DRtoSR: bool,
    logger: str,
) -> Dict[str, np.ndarray]:
    """
    This function calculates non closure corrections for a given process and splitting
    the calculation into different categories if necessary.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        corr_config: A dictionary with all the relevant information for calculating corrections to the measured fake factors
        sample_paths: List of file paths where the samples are stored
        output_path: Path where the generated plots should be stored
        process: This is relevant for the tt channel where two different fake factors are calculated, one for each hadronic tau
        closure_variable: Name of the variable the correction is calculated for
        evaluator: Evaluator with fake factors
        corr_evaluators: List of evaluators with corrections to fake factors
        for_DRtoSR: If True closure correction for the DR to SR correction fake factors will be calculated, if False for the general fake factors
        logger: Name of the logger that should be used

    Return:
        Dictionary of arrays with information about the smoothed function values to be stored with correctionlib (nominal and variations)
    """

    # get process specific config information
    if for_DRtoSR:
        correction_conf = corr_config["target_processes"][process]["DR_SR"]["non_closure"][closure_variable]
    else:
        correction_conf = corr_config["target_processes"][process]["non_closure"][closure_variable]

    if "split_categories" in correction_conf:
        split_variables, split_combinations, split_binnings, split_bandwidths = ff_func.get_split_combinations(
            categories=correction_conf["split_categories"],
            binning=correction_conf["var_bins"],
            bandwidth=correction_conf.get("bandwidth", None),
        )
    else:
        split_variables = []
        split_combinations = [None]
        split_binnings = [correction_conf["var_bins"]]
        split_bandwidths = [correction_conf.get("bandwidth", None)]

    results = func.optional_process_pool(
        args_list=[
            (
                split,
                binning,
                bandwidth,
                config,
                correction_conf,
                process,
                closure_variable,
                split_variables,
                sample_paths,
                output_path,
                logger,
                evaluator,
                corr_evaluators,
                for_DRtoSR,
            )
            for split, binning, bandwidth in zip(
                split_combinations,
                split_binnings,
                split_bandwidths,
            )
        ],
        function=NON_CLOSURE_CORRECTION_FUNCTIONS[process],
    )

    if len(split_combinations) == 1 and split_combinations[0] is None:
        return results[0]
    else:
        return ff_func.fill_corrlib_expression(results, split_variables)


def run_non_closure_correction(
    config: Dict[str, Union[str, Dict, List]],
    corr_config: Dict[str, Union[str, Dict]],
    process_config: Dict[str, Union[str, Dict]],
    process: str,
    evaluator: FakeFactorEvaluator,
    sample_paths: List[str],
    output_path: str,
    for_DRtoSR: bool,
    logger: str,
):
    """
    Function to calculate the non closure corrections for a given process.

    Args:
        config: Dictionary with information for fake factor calculation
        corr_config: A dictionary with information for calculating corrections
        process_config: A dictionary with information for calculating corrections for the specific process
        process: Name of the process
        evaluator: Fake factor evaluator
        sample_paths: List of file paths where samples are stored
        output_path: Path where generated plots are stored
        for_DRtoSR: If True closure correction for the DR to SR correction will be calculated
        logger: Name of the logger that should be used

    Return:
        Dictionary with the process name as key and a dictionary with the corrections
    """

    log = logging.getLogger(logger)
    corrections = {process: {}}
    all_non_closure_corr_vars, correction_set = [], None
    for idx, (closure_corr, _corr_config) in enumerate(
        process_config["non_closure"].items(),
    ):
        temp_conf = copy.deepcopy(config)

        if for_DRtoSR:
            log.info(f"Calculating closure correction for the DR to SR correction of the {process} process dependent on {closure_corr}.")
            log.info("-" * 50)
            func.modify_config(
                config=temp_conf,
                corr_config=process_config,
                process=process,
                to_AR_SR=False,
            )

        if "split_categories" in _corr_config:
            split_variables = list(_corr_config["split_categories"].keys())
            assert len(split_variables) == 1, "Only one split variable is supported"
            all_non_closure_corr_vars.append((closure_corr, split_variables[0]))
        else:
            all_non_closure_corr_vars.append(closure_corr)

        func.modify_config(
            config=temp_conf,
            corr_config=process_config["non_closure"][closure_corr],
            process=process,
            to_AR_SR=False,
        )

        corr_evaluators = []

        for n in range(idx):
            assert correction_set is not None, "Correction set must be calculated first! - This should not have happened!"
            corr_evaluators.append(
                FakeFactorCorrectionEvaluator.loading_from_CorrectionSet(
                    correction=correction_set,
                    process=process,
                    corr_variable=all_non_closure_corr_vars[n],
                    for_DRtoSR=for_DRtoSR,
                    logger=f"ff_corrections.{process}",
                )
            )

        result = non_closure_correction(
            config=temp_conf,
            corr_config=corr_config,
            sample_paths=sample_paths,
            output_path=output_path,
            process=process,
            closure_variable=closure_corr,
            evaluator=evaluator,
            corr_evaluators=corr_evaluators,
            for_DRtoSR=for_DRtoSR,
            logger=f"ff_corrections.{process}",
        )

        corrections[process]["non_closure_" + closure_corr] = result

        correction_set = corrlib.generate_correction_corrlib(
            config=corr_config,
            corrections=corrections,
            for_DRtoSR=for_DRtoSR,
        )

    return corrections


def run_ff_calculation_for_DRtoSR(
    args: Tuple[
        str,
        Dict[str, Union[Dict, List, str]],
        Dict[str, Union[Dict, str]],
        List[str],
        str,
    ]
) -> Union[Dict, None]:
    """
    This function can be used for multiprocessing. It runs the fake factor calculation step for a specified process for the DR to SR correction region.

    Args:
        args: Tuple with a process name, a configuration for this process, a configration for the correction for this process, a list of all sample paths and a path for the output

    Return:
        If a DR to SR correction is defined for the "process" a dictionary with fake factor function expressions is returned, otherwise None is returned
    """
    process, config, corr_config, sample_paths, output_path = args
    log = logging.getLogger(f"ff_corrections.{process}")

    if "DR_SR" in corr_config["target_processes"][process]:
        temp_conf = copy.deepcopy(config)
        func.modify_config(
            config=temp_conf,
            corr_config=corr_config["target_processes"][process]["DR_SR"],
            process=process,
            to_AR_SR=False,
        )
        log.info(f"Calculating fake factors for the SR-DR correction for the {process} process.")
        log.info("-" * 50)
        result = FF_calculation(
            config=temp_conf,
            sample_paths=sample_paths,
            output_path=output_path,
            process=process,
            logger=f"ff_corrections.{process}",
        )
    else:
        log.info(f"Target process {process} has no DR to SR calculation defined.")
        result = None

    return args, result


def run_non_closure_correction_for_DRtoSR(
    args: Tuple[
        str,
        Dict[str, Union[Dict, List, str]],
        Dict[str, Union[Dict, str]],
        List[str],
        str,
    ]
) -> Union[Dict, None]:
    """
    This function can be used for multiprocessing. It prepares and runs the non closure correction calculation step for a specified process for the DR to SR correction region.

    Args:
        args: Tuple with a process name, a configuration for this process, a configration for the correction for this process, a list of all sample paths and a path for the output
            process: Name of the process
            config: Dictionary with information for fake factor calculation
            corr_config: Dictionary with information for the correction calculation
            sample_paths: List of file paths where samples are stored
            output_path: Path where generated plots are

    Return:
        If a non closure correction for the DR to SR correction is defined for the "process"
        a dictionary with correction values is returned, otherwise None is returned
    """

    process, config, corr_config, sample_paths, output_path = args
    log = logging.getLogger(f"ff_corrections.{process}")
    corrections = {process: dict()}

    process_config = deepcopy(corr_config["target_processes"][process])
    if ("DR_SR" in process_config and "non_closure" in process_config["DR_SR"]):
        assert "ttbar" not in process, "ttbar is not supported for DR to SR corrections"

        var_dependences = [config["target_processes"][process]["var_dependence"]] + list(config["target_processes"][process]["split_categories"].keys())
        evaluator = FakeFactorEvaluator.loading_from_file(
            config=config,
            process=process,
            var_dependences=var_dependences,
            for_DRtoSR=True,
            logger=f"ff_corrections.{process}",
        )

        corrections.update(
            run_non_closure_correction(
                config=config,
                corr_config=corr_config,
                process_config=process_config["DR_SR"],
                process=process,
                evaluator=evaluator,
                sample_paths=sample_paths,
                output_path=output_path,
                for_DRtoSR=True,
                logger=f"ff_corrections.{process}",
            )
        )
    else:
        log.info(f"Target process {process} has no closure correction for DR to SR calculation defined.")

    return args, corrections


def run_correction(
    args,
) -> Dict[str, Dict[str, Any]]:
    """
    Function to calculate the fake factor corrections for a given process.

    Args:
        args: Tuple containing all the necessary information for the calculation of the fake factors
            process: Name of the process
            config: Dictionary with information for fake factor calculation
            corr_config: Dictionary with information for the correction calculation
            sample_paths: List of file paths where samples are stored
            save_path_plots: Path where generated plots are stored

    Return:
        Dictionary with the process name as key and a dictionary with the corrections
    """
    (
        process,
        config,
        corr_config,
        sample_paths,
        save_path_plots,
    ) = args

    log = logging.getLogger(f"ff_corrections.{process}")
    corrections = {process: dict()}

    var_dependences = [config["target_processes"][process]["var_dependence"]] + list(config["target_processes"][process]["split_categories"].keys())

    if "non_closure" in corr_config["target_processes"][process]:
        evaluator = FakeFactorEvaluator.loading_from_file(
            config=config,
            process=process,
            var_dependences=var_dependences,
            for_DRtoSR=False,
            logger=f"ff_corrections.{process}",
        )

        corrections.update(
            run_non_closure_correction(
                config=config,
                corr_config=corr_config,
                process_config=corr_config["target_processes"][process],
                process=process,
                evaluator=evaluator,
                sample_paths=sample_paths,
                output_path=save_path_plots,
                for_DRtoSR=False,
                logger=f"ff_corrections.{process}",
            )
        )

    if "DR_SR" in corr_config["target_processes"][process]:
        evaluator = FakeFactorEvaluator.loading_from_file(
            config=config,
            process=process,
            var_dependences=var_dependences,
            for_DRtoSR=True,
            logger=f"ff_corrections.{process}",
        )

        corr_evaluators = []
        DR_SR_conf = corr_config["target_processes"][process]["DR_SR"]

        for corr_var in DR_SR_conf["non_closure"].keys():
            non_closure_corr_vars_DR_SR = corr_var
            if "split_categories" in DR_SR_conf["non_closure"][corr_var]:
                split_variables = list(DR_SR_conf["non_closure"][corr_var]["split_categories"].keys())
                assert len(split_variables) == 1, "Only one split variable is supported"
                non_closure_corr_vars_DR_SR = (corr_var, split_variables[0])

            corr_evaluators.append(
                FakeFactorCorrectionEvaluator.loading_from_file(
                    config=config,
                    process=process,
                    corr_variable=non_closure_corr_vars_DR_SR,
                    for_DRtoSR=True,
                    logger=f"ff_corrections.{process}",
                )
            )

        log.info(f"Calculating DR to SR correction for the {process} process.")
        log.info("-" * 50)
        temp_conf = copy.deepcopy(config)
        func.modify_config(
            config=temp_conf,
            corr_config=DR_SR_conf,
            process=process,
            to_AR_SR=True,
        )

        if "split_categories" in DR_SR_conf:
            split_variables, split_combinations, split_binnings = ff_func.get_split_combinations(
                categories=DR_SR_conf["split_categories"],
                binning=DR_SR_conf["var_bins"],
            )
        else:
            split_variables = []
            split_combinations = [None]
            split_binnings = [DR_SR_conf["var_bins"]]

        results = func.optional_process_pool(
            args_list=[
                (
                    split,
                    binning,
                    config,
                    DR_SR_conf,
                    process,
                    split_variables,
                    sample_paths,
                    save_path_plots,
                    f"ff_corrections.{process}",
                    evaluator,
                    corr_evaluators,
                )
                for split, binning in zip(split_combinations, split_binnings)
            ],
            function=DR_SR_CORRECTION_FUNCTIONS[process],
        )

        if len(split_combinations) == 1 and split_combinations[0] is None:
            corrections[process]["DR_SR"] = results[0]
        else:
            corrections[process]["DR_SR"] = ff_func.fill_corrlib_expression(results, split_variables)

    return corrections


if __name__ == "__main__":
    args = parser.parse_args()

    func.RuntimeVariables.USE_MULTIPROCESSING = not args.disable_multiprocessing

    # loading of the chosen config file
    corr_config = func.load_config(args.config_file)
    workdir_path = os.path.join("workdir", corr_config["workdir_name"], corr_config["era"])
    func.check_path(path=os.path.join(os.getcwd(), workdir_path))

    with open(os.path.join(workdir_path, "fake_factors", corr_config["channel"], "config.yaml"), "r") as file:
        config = yaml.load(file, yaml.FullLoader)

    save_path_plots = os.path.join(workdir_path, "corrections", config["channel"])
    func.check_path(path=os.path.join(os.getcwd(), save_path_plots))

    with open(save_path_plots + "/config.yaml", "w") as config_file:
        yaml.dump(corr_config, config_file, default_flow_style=False, sort_keys=False)

    # start output logging
    func.setup_logger(
        log_file=save_path_plots + "/ff_corrections.log",
        log_name="ff_corrections",
        subcategories=corr_config["target_processes"].keys(),
    )

    # getting all the input files
    sample_paths = func.get_samples(config=config)
    if len(sample_paths) == 0:
        raise Exception("No input files found!")

    func.RuntimeVariables.RDataFrameWrapper = PassThroughWrapper
    if config.get("use_center_of_mass_bins", True):
        func.RuntimeVariables.RDataFrameWrapper = Histo1DPatchedRDataFrame

    ########### needed precalculations for DR to SR corrections ###########

    if not args.only_main_corrections:
        # initializing the fake factor calculation for DR to SR corrections
        if not args.skip_DRtoSR_ffs:
            fake_factors = dict()

            if "target_processes" in corr_config:
                results = func.optional_process_pool(
                    args_list=[
                        (process, config, corr_config, sample_paths, save_path_plots)
                        for process in corr_config["target_processes"]
                    ],
                    function=run_ff_calculation_for_DRtoSR,
                )

                for args, result in results:
                    if result is not None:
                        fake_factors[args[0]] = result
            else:
                raise Exception("No target processes are defined!")

            corrlib.generate_ff_corrlib_json(
                config=config,
                ff_functions=fake_factors,
                fractions=None,
                fractions_subleading=None,
                output_path=save_path_plots,
                for_corrections=True,
            )
        else:
            assert os.path.exists(
                os.path.join(
                    workdir_path,
                    "corrections",
                    config["channel"],
                    f"fake_factors_{config['channel']}_for_corrections.json",
                ),
            ), "Fake factors for DR to SR corrections are missing!"

        DR_SR_corrections = {
            "QCD": {},
            "QCD_subleading": {},
            "Wjets": {},
        }

        if "target_processes" in corr_config:
            results = func.optional_process_pool(
                args_list=[
                    (process, config, corr_config, sample_paths, save_path_plots)
                    for process in corr_config["target_processes"]
                ],
                function=run_non_closure_correction_for_DRtoSR,
            )

            for args, result in results:
                if result is not None:
                    DR_SR_corrections.update(result)
        else:
            raise Exception("No target processes are defined!")

        corrlib.generate_correction_corrlib(
            config=corr_config,
            corrections=func.remove_empty_keys(DR_SR_corrections),
            output_path=save_path_plots,
            for_DRtoSR=True,
        )

    ########### real fake factor corrections ###########

    corrections = {
        "QCD": {},
        "QCD_subleading": {},
        "Wjets": {},
        "ttbar": {},
        "ttbar_subleading": {},
    }

    if "target_processes" in corr_config:
        results = func.optional_process_pool(
            args_list=[
                (process, config, corr_config, sample_paths, save_path_plots)
                for process in corr_config["target_processes"]
            ],
            function=run_correction,
        )

        for result in results:
            corrections.update(result)

    corrlib.generate_correction_corrlib(
        config=corr_config,
        corrections=func.remove_empty_keys(corrections),
        output_path=workdir_path,
        for_DRtoSR=False,
    )

    with open(os.path.join(save_path_plots, "done"), "w") as done_file:
        done_file.write("")
