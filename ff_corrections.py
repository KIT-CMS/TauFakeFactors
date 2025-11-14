"""
Main script initializing the fake factor correction calculation
"""

import argparse
import contextlib
import copy
import logging
import os
import pickle
from copy import deepcopy
from typing import Any, Dict, List, Tuple, Union

import numpy as np

import FF_calculation.FF_QCD as FF_QCD
import FF_calculation.FF_ttbar as FF_ttbar
import FF_calculation.FF_Wjets as FF_Wjets
import helper.correctionlib_json as corrlib
import helper.ff_functions as ff_func
import helper.functions as func
from ff_calculation import FF_calculation
from helper.ff_evaluators import FakeFactorCorrectionEvaluator, FakeFactorEvaluator, DRSRCorrectionEvaluator
from helper.hooks_and_patches import Histo1DPatchedRDataFrame, PassThroughWrapper

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config-file",
    default=None,
    help="Path to the config file which contains information for the fake factor corrections step.",
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
        Flag to use cached ROOT RDataFrames and full non closure, FF and DRtoSR corrections if
        available. Will check for existance assuming the same order of all previous
        correction calculations.
    """,
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

    split_collections = ff_func.SplitQuantities(correction_conf)

    results = func.optional_process_pool(
        args_list=[
            (
                split_collection,
                config,
                correction_conf,
                process,
                closure_variable,
                sample_paths,
                output_path,
                f"{logger}.{correction_conf['var_dependence']}{'.DR_SR' if for_DRtoSR else ''}",
                evaluator,
                corr_evaluators,
                for_DRtoSR,
            )
            for split_collection in split_collections
        ],
        function=NON_CLOSURE_CORRECTION_FUNCTIONS[process],
    )

    if len(split_collections) == 1 and split_collections.split[0] is None:
        return results[0]
    else:
        return ff_func.fill_corrlib_expression(results, split_collections.split_variables)


def run_non_closure_correction(
    config: Dict[str, Union[str, Dict, List]],
    corr_config: Dict[str, Union[str, Dict]],
    process: str,
    evaluator: FakeFactorEvaluator,
    sample_paths: List[str],
    output_path: str,
    for_DRtoSR: bool,
    DR_SR_evaluator: Union[FakeFactorCorrectionEvaluator, None],
    logger: str,
):
    """
    Function to calculate the non closure corrections for a given process.
    If DR_SR_evaluator is given and for_DRtoSR is False the non closure corrections will be
    calculated including the DR to SR correction applied first.

    Args:
        config: Dictionary with information for fake factor calculation
        corr_config: A dictionary with information for calculating corrections
        process_config: A dictionary with information for calculating corrections for the specific process
        process: Name of the process
        evaluator: Fake factor evaluator
        sample_paths: List of file paths where samples are stored
        output_path: Path where generated plots are stored
        for_DRtoSR: If True closure correction for the DR to SR correction will be calculated
        DR_SR_evaluator: Evaluator with the DR to SR correction if it is already calculated
        logger: Name of the logger that should be used

    Return:
        Dictionary with the process name as key and a dictionary with the corrections
    """

    log = logging.getLogger(logger)
    corrections = {process: {}}
    _chained_DR_SR_process_config = None
    if for_DRtoSR:
        process_config = deepcopy(corr_config["target_processes"][process]["DR_SR"])
    else:
        process_config = deepcopy(corr_config["target_processes"][process])
        if DR_SR_evaluator is None:
            _chained_DR_SR_process_config = deepcopy(corr_config["target_processes"][process].get("DR_SR"))

    all_non_closure_corr_vars, correction_set, is_valid_cache = [], None, True
    for idx, (closure_var, closure_var_config) in enumerate(
        process_config["non_closure"].items(),
    ):
        ff_config = copy.deepcopy(config)

        if for_DRtoSR:
            log.info(f"Calculating closure correction for the DR to SR correction of the {process} process dependent on {closure_var}.")
            log.info("-" * 50)
            func.modify_config(
                config=ff_config,
                corr_config=process_config,
                process=process,
                to_AR_SR=False,
            )

        if "split_categories" in closure_var_config:
            split_variables = list(closure_var_config["split_categories"].keys())
            assert len(split_variables) == 1, "Only one split variable is supported"
            all_non_closure_corr_vars.append((closure_var, split_variables[0]))
        else:
            all_non_closure_corr_vars.append(closure_var)

        func.modify_config(
            config=ff_config,
            corr_config=closure_var_config,
            process=process,
            to_AR_SR=False,
        )

        cached_non_closure = func.get_cached_file_path(
            output_path=output_path,
            process=process,
            variables=all_non_closure_corr_vars,
            for_DRtoSR=for_DRtoSR,
        )

        if os.path.exists(cached_non_closure):
            with open(cached_non_closure, "rb") as f:
                __corrections, __corr_config, __ff_config, __chained_DR_SR_process_config = pickle.load(f)
                is_valid_cache &= func.correction_config_comparison(
                    __corr_config,
                    corr_config,
                    process=process,
                    closure_corr=closure_var,
                    for_DRtoSR=for_DRtoSR,
                )
                # fake factor config comparison is done for the case that
                # the fake factors are rerun with changes with the same tag
                is_valid_cache &= func.nested_object_comparison(
                    __ff_config,
                    ff_config,
                )
                is_valid_cache &= func.nested_object_comparison(
                    __chained_DR_SR_process_config,
                    _chained_DR_SR_process_config,
                )
        else:
            is_valid_cache = False

        if func.RuntimeVariables.USE_CACHED_INTERMEDIATE_STEPS and is_valid_cache:
            corrections[process].update(__corrections[process])
            correction_set = corrlib.generate_correction_corrlib(
                config=corr_config,
                corrections=corrections,
                for_DRtoSR=for_DRtoSR,
            )
        else:
            if for_DRtoSR:
                log.info(f"Removing cached DR_SR corrections for {process} process due to changes in non closure corrections")
                with contextlib.suppress(FileNotFoundError):
                    os.remove(func.get_cached_file_path(output_path=output_path, process=process))

            corr_evaluators = [deepcopy(DR_SR_evaluator)] if DR_SR_evaluator is not None else []

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
                config=ff_config,
                corr_config=corr_config,
                sample_paths=sample_paths,
                output_path=output_path,
                process=process,
                closure_variable=closure_var,
                evaluator=evaluator,
                corr_evaluators=corr_evaluators,
                for_DRtoSR=for_DRtoSR,
                logger=f"ff_corrections.{process}",
            )

            corrections[process]["non_closure_" + closure_var] = result

            correction_set = corrlib.generate_correction_corrlib(
                config=corr_config,
                corrections=corrections,
                for_DRtoSR=for_DRtoSR,
            )

            with open(cached_non_closure, "wb") as f:
                pickle.dump(
                    tuple(
                        [
                            corrections,
                            corr_config,
                            ff_config,
                            _chained_DR_SR_process_config,
                        ]
                    ),
                    f,
                    protocol=pickle.HIGHEST_PROTOCOL,
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
        ff_config = copy.deepcopy(config)
        func.modify_config(
            config=ff_config,
            corr_config=corr_config["target_processes"][process]["DR_SR"],
            process=process,
            to_AR_SR=False,
        )
        log.info(f"Calculating fake factors for the DR to SR correction for the {process} process.")
        log.info("-" * 50)
        result = FF_calculation(
            config=ff_config,
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
                process=process,
                evaluator=evaluator,
                sample_paths=sample_paths,
                output_path=output_path,
                for_DRtoSR=True,
                DR_SR_evaluator=None,
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
            save_path: Path where generated results are stored

    Return:
        Dictionary with the process name as key and a dictionary with the corrections
    """
    (
        process,
        config,
        corr_config,
        sample_paths,
        save_path,
    ) = args

    log = logging.getLogger(f"ff_corrections.{process}")
    corrections = {process: dict()}

    var_dependences = [config["target_processes"][process]["var_dependence"]] + list(config["target_processes"][process]["split_categories"].keys())

    chain_DR_SR = all(
        [
            "DR_SR" in corr_config["target_processes"][process],
            corr_config["target_processes"][process].get("chain_DR_SR_to_non_closure", False)
        ]
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
        DR_SR_config = corr_config["target_processes"][process]["DR_SR"]

        for corr_var in DR_SR_config["non_closure"].keys():
            non_closure_corr_vars_DR_SR = corr_var
            if "split_categories" in DR_SR_config["non_closure"][corr_var]:
                split_variables = list(DR_SR_config["non_closure"][corr_var]["split_categories"].keys())
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
        ff_config = copy.deepcopy(config)
        func.modify_config(
            config=ff_config,
            corr_config=DR_SR_config,
            process=process,
            to_AR_SR=True,
        )

        split_collections = ff_func.SplitQuantities(DR_SR_config)

        is_valid_cache = True
        cached_DR_SR = func.get_cached_file_path(output_path=save_path, process=process)
        if os.path.exists(cached_DR_SR):
            with open(cached_DR_SR, "rb") as f:
                __corrections, __DR_SR_config = pickle.load(f)
                is_valid_cache = func.nested_object_comparison(
                    __DR_SR_config,
                    DR_SR_config,
                )
        else:
            is_valid_cache = False

        if func.RuntimeVariables.USE_CACHED_INTERMEDIATE_STEPS and is_valid_cache:
            log.info(f"Using cached DR to SR corrections for {process} process.")
            corrections[process]["DR_SR"] = __corrections[process]["DR_SR"]
        else:
            results = func.optional_process_pool(
                args_list=[
                    (
                        split_collection,
                        config,
                        DR_SR_config,
                        process,
                        sample_paths,
                        save_path,
                        f"ff_corrections.{process}",
                        evaluator,
                        corr_evaluators,
                    )
                    for split_collection in split_collections
                ],
                function=DR_SR_CORRECTION_FUNCTIONS[process],
            )

            if len(split_collections) == 1 and split_collections.split[0] is None:
                corrections[process]["DR_SR"] = results[0]
            else:
                corrections[process]["DR_SR"] = ff_func.fill_corrlib_expression(results, split_collections.split_variables)

            with open(cached_DR_SR, "wb") as f:
                pickle.dump(
                    (corrections, DR_SR_config),
                    f,
                    protocol=pickle.HIGHEST_PROTOCOL,
                )

        DR_SR_correction = DRSRCorrectionEvaluator.loading_from_CorrectionSet(
            correction=corrlib.generate_correction_corrlib(
                config=corr_config,
                corrections=func.remove_empty_keys(corrections),
                for_DRtoSR=True,
            ),
            process=process,
            corr_variable=tuple(
                [
                    corr_config["target_processes"][process]["DR_SR"]["var_dependence"],
                ] + list(corr_config["target_processes"][process]["DR_SR"].get("split_categories", {}).keys())
            ),
            logger=f"ff_corrections.{process}",
        )

    if "non_closure" in corr_config["target_processes"][process]:
        evaluator = FakeFactorEvaluator.loading_from_file(
            config=config,
            process=process,
            var_dependences=var_dependences,
            for_DRtoSR=False,
            logger=f"ff_corrections.{process}",
        )

        corrections[process].update(
            run_non_closure_correction(
                config=config,
                corr_config=corr_config,
                process=process,
                evaluator=evaluator,
                sample_paths=sample_paths,
                output_path=save_path,
                for_DRtoSR=False,
                DR_SR_evaluator=DR_SR_correction if chain_DR_SR else None,
                logger=f"ff_corrections.{process}",
            )[process]
        )

    return corrections


if __name__ == "__main__":
    args = parser.parse_args()

    func.RuntimeVariables.USE_MULTIPROCESSING = not args.disable_multiprocessing
    func.RuntimeVariables.USE_CACHED_INTERMEDIATE_STEPS = not args.ignore_cached_intermediary_steps

    # loading of the chosen config file
    corr_config = func.load_config(args.config_file)
    workdir_path = os.path.join("workdir", corr_config["workdir_name"], corr_config["era"])
    func.check_path(path=os.path.join(os.getcwd(), workdir_path))

    with open(os.path.join(workdir_path, "fake_factors", corr_config["channel"], "config.yaml"), "r") as file:
        config = func.configured_yaml.load(file)

    save_path = os.path.join(workdir_path, "corrections", config["channel"])
    func.check_path(path=os.path.join(os.getcwd(), save_path))

    with open(save_path + "/config.yaml", "w") as config_file:
        func.configured_yaml.dump(corr_config, config_file)

    # start output logging
    func.setup_logger(
        log_file=save_path + "/ff_corrections.log",
        log_name="ff_corrections",
        log_level=logging.INFO,
        subcategories=corr_config["target_processes"].keys(),
    )

    # getting all the input files
    sample_paths = func.get_samples(config=config)
    if len(sample_paths) == 0:
        raise Exception("No input files found!")

    func.RuntimeVariables.RDataFrameWrapper = PassThroughWrapper
    if config.get("use_center_of_mass_bins", True):
        func.RuntimeVariables.RDataFrameWrapper = Histo1DPatchedRDataFrame

    # setting default systematic variations if not present in the config
    if "correction_variations" not in corr_config:
        corr_config["correction_variations"] = ("Stat_1Sigma", "Syst_MCShift", "Syst_BandAsym")

    ########### needed precalculations for DR to SR corrections ###########

    # initializing the fake factor calculation for DR to SR corrections
    ff_for_DRtoSR_file = os.path.join(
        workdir_path,
        "corrections",
        config["channel"],
        f"fake_factors_{config['channel']}_for_corrections.json",
    )

    is_valid_cache = True
    cached_DR_SR_ffs = func.get_cached_file_path(
        output_path=save_path,
    )

    if os.path.exists(cached_DR_SR_ffs):
        with open(cached_DR_SR_ffs, "rb") as f:
            __ffs, __corr_config = pickle.load(f)
            for proc in corr_config["target_processes"]:
                if "DR_SR" in corr_config["target_processes"][proc]:
                    __test_config = __corr_config["target_processes"][proc].get("DR_SR", {})
                    test_config = corr_config["target_processes"][proc].get("DR_SR", {})

                    is_valid_cache = all(
                        func.nested_object_comparison(__test_config[k], test_config[k])
                        for k in ("SRlike_cuts", "ARlike_cuts", "AR_SR_cuts")
                    )
    else:
        is_valid_cache = False

    log = logging.getLogger("ff_corrections")
    if func.RuntimeVariables.USE_CACHED_INTERMEDIATE_STEPS and is_valid_cache:
        log.info("Using cached DR to SR fake factors.")
        fake_factors = __ffs
    else:
        fake_factors = dict()

        if "target_processes" in corr_config:
            results = func.optional_process_pool(
                args_list=[
                    (process, config, corr_config, sample_paths, save_path)
                    for process in corr_config["target_processes"]
                ],
                function=run_ff_calculation_for_DRtoSR,
            )

            for args, result in results:
                if result is not None:
                    fake_factors[args[0]] = result
        else:
            raise Exception("No target processes are defined!")

        if fake_factors:
            corrlib.generate_ff_corrlib_json(
                config=config,
                ff_functions=fake_factors,
                fractions=None,
                fractions_subleading=None,
                output_path=save_path,
                for_corrections=True,
            )

        with open(cached_DR_SR_ffs, "wb") as f:
            pickle.dump(
                (fake_factors, corr_config),
                f,
                protocol=pickle.HIGHEST_PROTOCOL,
            )

    DR_SR_corrections = {
        "QCD": {},
        "QCD_subleading": {},
        "Wjets": {},
    }

    if "target_processes" in corr_config:
        results = func.optional_process_pool(
            args_list=[
                (process, config, corr_config, sample_paths, save_path)
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
        output_path=save_path,
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
                (process, config, corr_config, sample_paths, save_path)
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

    with open(os.path.join(save_path, "done"), "w") as done_file:
        done_file.write("")
