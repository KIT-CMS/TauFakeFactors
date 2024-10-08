"""
Main script initializing the fake factor correction calculation
"""

import argparse
import yaml
import os
import copy
import logging
from typing import Tuple, Dict, List, Union
import concurrent.futures

import FF_calculation.FF_QCD as FF_QCD
import FF_calculation.FF_Wjets as FF_Wjets
import FF_calculation.FF_ttbar as FF_ttbar
import helper.correctionlib_json as corrlib
import helper.functions as func
from helper.ff_evaluators import FakeFactorEvaluator, FakeFactorCorrectionEvaluator

parser = argparse.ArgumentParser()

parser.add_argument(
    "--config-file",
    default=None,
    help="Path to the config file which contains information for the fake factor corrections step.",
)
parser.add_argument(
    "--common-config-file",
    default=None,
    help="""
        Path to the common config file which contains common information i.e. ntuple path, workdir etc.
        Settings loaded here are overwritten by the settings in --config-file if present there.
        'common_settings.yaml' in the same folder as --config-file is loaded by default if present.
    """,
)
parser.add_argument(
    "--only-main-corrections",
    action="store_true",
    help="Using this argument means to skip the calculation of the fake factors and corrections for the DR to SR correction and directly calculate the main corrections. This is useful if you already calculated the needed additional fake factors.",
)
parser.add_argument(
    "--skip-DRtoSR-ffs",
    action="store_true",
    help="Using this argument means to skip the calculation of the fake factors for the DR to SR correction and directly calculate the corrections. This is useful if you already calculated the needed additional fake factors.",
)


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
    This function can be used for multiprocessing. It runs the fake factor calculation step for a specified process.

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

        if process == "QCD":
            log.info(
                "Calculating fake factors for the SR-DR correction for the QCD process."
            )
            log.info("-" * 50)
            result = FF_QCD.calculation_QCD_FFs(
                config=temp_conf,
                sample_paths=sample_paths,
                output_path=output_path,
                process=process,
                logger=f"ff_corrections.{process}",
            )
        elif process == "QCD_subleading":
            log.info(
                "Calculating fake factors for the SR-DR correction for the QCD(subleading) process."
            )
            log.info("-" * 50)
            result = FF_QCD.calculation_QCD_FFs(
                config=temp_conf,
                sample_paths=sample_paths,
                output_path=output_path,
                process=process,
                logger=f"ff_corrections.{process}",
            )
        elif process == "Wjets":
            log.info(
                "Calculating fake factors for the SR-DR correction for the Wjets process."
            )
            log.info("-" * 50)
            result = FF_Wjets.calculation_Wjets_FFs(
                config=temp_conf,
                sample_paths=sample_paths,
                output_path=output_path,
                logger=f"ff_corrections.{process}",
            )
        else:
            log.info(f"No DR to SR calculation function is defined for {process}.")
            result = None

    else:
        log.info(f"Target process {process} has no DR to SR calculation defined.")
        result = None

    return result


def run_non_closure_for_DRtoSR(
    args: Tuple[
        str,
        Dict[str, Union[Dict, List, str]],
        Dict[str, Union[Dict, str]],
        List[str],
        str,
    ]
) -> Union[Dict, None]:
    """
    This function can be used for multiprocessing. It runs the non closure correction calculation step for a specified process.

    Args:
        args: Tuple with a process name, a configuration for this process, a configration for the correction for this process, a list of all sample paths and a path for the output

    Return:
        If a non closure correction for the DR to SR correction is defined for the "process"
        a dictionary with correction values is returned, otherwise None is returned
    """

    process, config, corr_config, sample_paths, output_path = args
    log = logging.getLogger(f"ff_corrections.{process}")

    if (
        "DR_SR" in corr_config["target_processes"][process]
        and "non_closure" in corr_config["target_processes"][process]["DR_SR"]
    ):
        evaluator = FakeFactorEvaluator(
            config=config,
            process=process,
            for_DRtoSR=True,
            logger=f"ff_corrections.{process}",
        )
        temp_conf = copy.deepcopy(config)
        func.modify_config(
            config=temp_conf,
            corr_config=corr_config["target_processes"][process]["DR_SR"],
            process=process,
            to_AR_SR=False,
        )
        
        closure_correction = list(
            corr_config["target_processes"][process]["DR_SR"]["non_closure"]
        )[0]
        func.modify_config(
            config=temp_conf,
            corr_config=corr_config["target_processes"][process]["DR_SR"]["non_closure"][closure_correction],
            process=process,
            to_AR_SR=False,
        )

        if process == "QCD":
            log.info(
                f"Calculating closure correction for the DR to SR correction of the QCD process dependent on {closure_correction}."
            )
            log.info("-" * 50)
            result = FF_QCD.non_closure_correction(
                config=temp_conf,
                corr_config=corr_config,
                sample_paths=sample_paths,
                output_path=output_path,
                process=process,
                closure_variable=closure_correction,
                evaluator=evaluator,
                corr_evaluator=None,
                for_DRtoSR=True,
                logger=f"ff_corrections.{process}",
            )

        elif process == "QCD_subleading":
            log.info(
                f"Calculating closure correction for the DR to SR correction of the QCD(subleading) process dependent on {closure_correction}."
            )
            log.info("-" * 50)
            result = FF_QCD.non_closure_correction(
                config=temp_conf,
                corr_config=corr_config,
                sample_paths=sample_paths,
                output_path=output_path,
                process=process,
                closure_variable=closure_correction,
                evaluator=evaluator,
                corr_evaluator=None,
                for_DRtoSR=True,
                logger=f"ff_corrections.{process}",
            )

        elif process == "Wjets":
            log.info(
                f"Calculating closure correction for the DR to SR correction of the Wjets process dependent on {closure_correction}."
            )
            log.info("-" * 50)
            result = FF_Wjets.non_closure_correction(
                config=temp_conf,
                corr_config=corr_config,
                sample_paths=sample_paths,
                output_path=output_path,
                closure_variable=closure_correction,
                evaluator=evaluator,
                corr_evaluator=None,
                for_DRtoSR=True,
                logger=f"ff_corrections.{process}",
            )
        else:
            log.info(
                f"No closure correction function for DR to SR calculation is defined for {process}."
            )
            result = None

    else:
        log.info(
            f"Target process {process} has no closure correction for DR to SR calculation defined."
        )
        result = None
        closure_correction = None

    return result, closure_correction


if __name__ == "__main__":
    args = parser.parse_args()

    # loading of the chosen config file
    config = func.load_config(args.config_file, args.common_config_file)

    with open(
        os.path.join(
            "workdir",
            corr_config["workdir_name"],
            corr_config["era"],
            "fake_factors",
            corr_config["channel"],
            "config.yaml",
        ),
        "r",
    ) as file:
        config = yaml.load(file, yaml.FullLoader)

    save_path_ffs = os.path.join("workdir", config["workdir_name"], config["era"])
    func.check_path(path=os.path.join(os.getcwd(), save_path_ffs))
    save_path_plots = os.path.join(
        "workdir",
        config["workdir_name"],
        config["era"],
        "corrections",
        config["channel"],
    )
    func.check_path(path=os.path.join(os.getcwd(), save_path_plots))

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

    ########### needed precalculations for DR to SR corrections ###########

    if not args.only_main_corrections:
        # initializing the fake factor calculation for DR to SR corrections
        if not args.skip_DRtoSR_ffs:
            fake_factors = dict()

            if "target_processes" in corr_config:
                args_list = [
                    (process, config, corr_config, sample_paths, save_path_plots)
                    for process in corr_config["target_processes"]
                ]

                with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
                    for args, result in zip(
                        args_list, executor.map(run_ff_calculation_for_DRtoSR, args_list)
                    ):
                        if result is not None:
                            fake_factors[args[0]] = result
            else:
                raise Exception("No target processes are defined!")

            corrlib.generate_ff_corrlib_json(
                config=config,
                ff_functions=fake_factors,
                fractions=None,
                output_path=save_path_plots,
                for_corrections=True,
            )

        DR_SR_corrections = {
            "QCD": dict(),
            "QCD_subleading": dict(),
            "Wjets": dict(),
            "ttbar": dict(),
        }

        if "target_processes" in corr_config:
            args_list = [
                (process, config, corr_config, sample_paths, save_path_plots)
                for process in corr_config["target_processes"]
            ]

            with concurrent.futures.ProcessPoolExecutor(max_workers=8) as executor:
                for args, result in zip(
                    args_list, executor.map(run_non_closure_for_DRtoSR, args_list)
                ):
                    if result[0] is not None:
                        DR_SR_corrections[args[0]]["non_closure_" + result[1]] = result[
                            0
                        ]
        else:
            raise Exception("No target processes are defined!")

        corrlib.generate_correction_corrlib_json(
            config=corr_config,
            corrections=DR_SR_corrections,
            output_path=save_path_plots,
            for_DRtoSR=True,
        )

    ########### real fake factor corrections ###########

    corrections = {
        "QCD": dict(),
        "QCD_subleading": dict(),
        "Wjets": dict(),
        "ttbar": dict(),
    }

    if "target_processes" in corr_config:
        for process in corr_config["target_processes"]:
            log = logging.getLogger(f"ff_corrections.{process}")

            if "non_closure" in corr_config["target_processes"][process]:
                evaluator = FakeFactorEvaluator(
                    config=config,
                    process=process,
                    for_DRtoSR=False,
                    logger=f"ff_corrections.{process}",
                )

                for idx, closure_corr in enumerate(
                    corr_config["target_processes"][process]["non_closure"]
                ):
                    log.info(
                        f"Calculating closure correction for the {process} process dependent on {closure_corr}."
                    )
                    log.info("-" * 50)
                    temp_conf = copy.deepcopy(config)
                    func.modify_config(
                        config=temp_conf,
                        corr_config=corr_config["target_processes"][process][
                            "non_closure"
                        ][closure_corr],
                        process=process,
                        to_AR_SR=False,
                    )

                    if idx == 0:
                        corr_evaluator = None
                    else:
                        corr_evaluator = FakeFactorCorrectionEvaluator(
                            config=config,
                            process=process,
                            for_DRtoSR=False,
                            logger=f"ff_corrections.{process}",
                        )

                    if process in ["QCD", "QCD_subleading"]:
                        corr = FF_QCD.non_closure_correction(
                            config=temp_conf,
                            corr_config=corr_config,
                            sample_paths=sample_paths,
                            output_path=save_path_plots,
                            process=process,
                            closure_variable=closure_corr,
                            evaluator=evaluator,
                            corr_evaluator=corr_evaluator,
                            for_DRtoSR=False,
                            logger=f"ff_corrections.{process}",
                        )
                    elif process == "Wjets":
                        corr = FF_Wjets.non_closure_correction(
                            config=temp_conf,
                            corr_config=corr_config,
                            sample_paths=sample_paths,
                            output_path=save_path_plots,
                            closure_variable=closure_corr,
                            evaluator=evaluator,
                            corr_evaluator=corr_evaluator,
                            for_DRtoSR=False,
                            logger=f"ff_corrections.{process}",
                        )
                    elif process == "ttbar":
                        corr = FF_ttbar.non_closure_correction(
                            config=temp_conf,
                            corr_config=corr_config,
                            sample_paths=sample_paths,
                            output_path=save_path_plots,
                            closure_variable=closure_corr,
                            evaluator=evaluator,
                            corr_evaluator=corr_evaluator,
                            logger=f"ff_corrections.{process}",
                        )
                    else:
                        raise ValueError(f"Process {process} not known!")

                    corrections[process]["non_closure_" + closure_corr] = corr

                    corrlib.generate_correction_corrlib_json(
                        config=corr_config,
                        corrections=corrections,
                        output_path=save_path_ffs,
                        for_DRtoSR=False,
                    )

            if "DR_SR" in corr_config["target_processes"][process]:
                evaluator = FakeFactorEvaluator(
                    config=config,
                    process=process,
                    for_DRtoSR=True,
                    logger=f"ff_corrections.{process}",
                )
                corr_evaluator = FakeFactorCorrectionEvaluator(
                    config=config,
                    process=process,
                    for_DRtoSR=True,
                    logger=f"ff_corrections.{process}",
                )

                log.info(f"Calculating DR to SR correction for the {process} process.")
                log.info("-" * 50)
                temp_conf = copy.deepcopy(config)
                func.modify_config(
                    config=temp_conf,
                    corr_config=corr_config["target_processes"][process]["DR_SR"],
                    process=process,
                    to_AR_SR=True,
                )

                if process in ["QCD", "QCD_subleading"]:
                    corr = FF_QCD.DR_SR_correction(
                        config=temp_conf,
                        corr_config=corr_config,
                        sample_paths=sample_paths,
                        output_path=save_path_plots,
                        process=process,
                        evaluator=evaluator,
                        corr_evaluator=corr_evaluator,
                        logger=f"ff_corrections.{process}",
                    )
                elif process == "Wjets":
                    corr = FF_Wjets.DR_SR_correction(
                        config=temp_conf,
                        corr_config=corr_config,
                        sample_paths=sample_paths,
                        output_path=save_path_plots,
                        evaluator=evaluator,
                        corr_evaluator=corr_evaluator,
                        logger=f"ff_corrections.{process}",
                    )
                else:
                    raise ValueError(
                        f"Process {process} not known or DR to SR correction not defined!"
                    )

                corrections[process]["DR_SR"] = corr

    corrlib.generate_correction_corrlib_json(
        config=corr_config,
        corrections=corrections,
        output_path=save_path_ffs,
        for_DRtoSR=False,
    )
