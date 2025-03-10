import logging
import os
from typing import Any, Dict, List, Union, Tuple

import correctionlib
import correctionlib.schemav2 as cs
import ROOT


class FakeFactorEvaluator:
    """
    Evaluator class to initiate a fake factor setup. The fake factors are loaded from an already produced correctionlib file.
    """

    @classmethod
    def loading_from_file(
        cls,
        config: Dict[str, Union[str, Dict, List]],
        process: str,
        var_dependences: List[str],
        for_DRtoSR: bool,
        logger: str,
    ) -> "FakeFactorEvaluator":
        log = logging.getLogger(logger)

        directories = ["workdir", config["workdir_name"], config["era"]]
        if not for_DRtoSR:
            path = os.path.join(*directories, f"fake_factors_{config['channel']}.json")
        else:
            path = os.path.join(
                *directories,
                "corrections",
                f"{config['channel']}",
                f"fake_factors_{config['channel']}_for_corrections.json",
            )

        correctionlib.register_pyroot_binding()

        log.info(f"Loading fake factor file {path} for process {process}")
        ROOT.gInterpreter.Declare(
            f'auto {process}_{"for_DRtoSR" if for_DRtoSR else ""} = '
            f'correction::CorrectionSet::from_file("{path}")->'
            f'at("{process}_fake_factors");'
        )

        return cls(process, var_dependences, for_DRtoSR)

    @classmethod
    def loading_from_CorrectionSet(
        cls,
        fake_factors: cs.CorrectionSet,
        process: str,
        var_dependences: List[str],
        for_DRtoSR: bool,
        logger: str,
    ) -> "FakeFactorEvaluator":
        log = logging.getLogger(logger)

        assert isinstance(fake_factors, cs.CorrectionSet), "face_factors must be of type correctionlib.schemav2.CorrectionSet"
        literal = fake_factors.json().replace('"', r'\"')

        log.info(f"Loading fake factor from string for process {process}")
        correctionlib.register_pyroot_binding()
        ROOT.gInterpreter.Declare(
            f'auto {process}_{"for_DRtoSR" if for_DRtoSR else ""} = '
            f'correction::CorrectionSet::from_string("{literal}")->'
            f'at("{process}_fake_factors");'
        )

        return cls(process, var_dependences, for_DRtoSR)

    def __init__(
        self,
        process: str,
        var_dependences: List[str],
        for_DRtoSR: bool,
    ):
        """
        Initiating a new evaluator for fake factors using correctionlib.

        Args:
            config: A dictionary with all the relevant information for the fake factor calculation
            process: Name of the process the fake factors were calculated for
            var_dependences: List of variable dependences of the fake factors
            for_DRtoSR: If True fake factors calculated specifically for the DR to SR correction will be loaded, if False the general fake factors will be used
            logger: Name of the logger that should be used
        """

        self.process = process
        self.for_DRtoSR = "for_DRtoSR" if for_DRtoSR else ""
        self.var_dependences = var_dependences

    @property
    def str_var_dependences(self) -> List[str]:
        return ", ".join([f'(float){var}' for var in self.var_dependences])

    def evaluate_fake_factor(self, rdf: Any) -> Any:
        """
        Evaluating the fake factors based on the variables it depends on.
        In this function it is the transverse momentum of the subleading lepton in the tau pair and the number of jets.

        Args:
            rdf: root DataFrame object

        Return:
            root DataFrame object with a new column with the evaluated fake factors
        """
        eval_str = self.str_var_dependences + ', "nominal"'
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{{eval_str}}})',
        )
        return rdf


class FakeFactorCorrectionEvaluator:
    """
    Evaluator class to initiate a fake factor correction setup. The fake factor corrections are loaded from an already produced correctionlib file.
    Currently only a readout for a correction dependent on the leading or subleading pt is implemented.
    """

    @classmethod
    def loading_from_file(
        cls,
        config: Dict[str, Union[str, Dict, List]],
        process: str,
        corr_variable: Union[str, Tuple[str, ...]],
        for_DRtoSR: bool,
        logger: str,
    ) -> "FakeFactorCorrectionEvaluator":
        log = logging.getLogger(logger)

        _for_DRtoSR = "for_DRtoSR" if for_DRtoSR else ""
        directories = ["workdir", config["workdir_name"], config["era"]]

        if not for_DRtoSR:
            path = os.path.join(*directories, f"FF_corrections_{config['channel']}.json")
        else:
            path = os.path.join(
                *directories,
                "corrections",
                f"{config['channel']}",
                f"FF_corrections_{config['channel']}_{_for_DRtoSR}.json",
            )

        variable = corr_variable if isinstance(corr_variable, str) else corr_variable[-1]

        correctionlib.register_pyroot_binding()
        log.info(f"Loading fake factor correction file {path} for process {process}")
        ROOT.gInterpreter.Declare(
            f'auto {process}_corr_{variable}_{_for_DRtoSR} = '
            f'correction::CorrectionSet::from_file("{path}")'
            f'->at("{process}_non_closure_{variable}_correction");'
        )

        return cls(process, corr_variable, for_DRtoSR, logger)

    @classmethod
    def loading_from_CorrectionSet(
        cls,
        correction: cs.CorrectionSet,
        process: str,
        corr_variable: Union[str, Tuple[str, ...]],
        for_DRtoSR: bool,
        logger: str,
    ) -> "FakeFactorCorrectionEvaluator":
        log = logging.getLogger(logger)

        assert isinstance(correction, cs.CorrectionSet), "Correction must be of type correctionlib.schemav2.CorrectionSet"

        variable = corr_variable if isinstance(corr_variable, str) else corr_variable[-1]
        literal = correction.json().replace('"', r'\"')

        log.info(f"Loading fake factor correction from string for process {process}")
        correctionlib.register_pyroot_binding()
        ROOT.gInterpreter.Declare(
            f'auto {process}_corr_{variable}_{"for_DRtoSR" if for_DRtoSR else ""} = '
            f'correction::CorrectionSet::from_string("{literal}")'
            f'->at("{process}_non_closure_{variable}_correction");'
        )

        return cls(process, corr_variable, for_DRtoSR, logger)

    def __init__(
        self,
        process: str,
        corr_variable: Union[str, Tuple[str, ...]],
        for_DRtoSR: bool,
    ):
        """
        Initiating a new evaluator for fake factor corrections using correctionlib.

        Args:
            config: A dictionary with all the relevant information for the fake factor correction calculation
            process: Name of the process the fake factor corrections were calculated for
            corr_variable: Name of the variable dependence of the correction
            for_DRtoSR: If True fake factor corrections calculated specifically for the DR to SR correction will be loaded, if False the general fake factor corrections will be used
            logger: Name of the logger that should be used
        """

        self.for_DRtoSR = "for_DRtoSR" if for_DRtoSR else ""

        self.process = process
        if not isinstance(corr_variable, str) and isinstance(corr_variable[-1], str):
            self.variable = corr_variable[-1]
            self.var_dependences = corr_variable
        else:
            self.variable = corr_variable
            self.var_dependences = [corr_variable]

    @property
    def str_var_dependences(self) -> List[str]:
        return ", ".join([f'(float){var}' for var in self.var_dependences])

    def evaluate_correction(self, rdf: Any) -> Any:
        """
        Evaluating the fake factor corrections based on the variables it depends on.

        Args:
            rdf: root DataFrame object

        Return:
            root DataFrame object with a new column with the evaluated fake factor corrections
        """
        eval_str = self.str_var_dependences + ', "nominal"'
        rdf = rdf.Define(
            f"{self.process}_ff_corr_{self.variable}",
            f"{self.process}_corr_{self.variable}_{self.for_DRtoSR}->evaluate({{{eval_str}}})",
        )
        return rdf
