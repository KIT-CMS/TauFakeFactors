import os
import ROOT
import correctionlib
import logging
from typing import Dict, Union, List, Any


class FakeFactorEvaluator:
    """
    Evaluator class to initiate a fake factor setup. The fake factors are loaded from an already produced correctionlib file.
    """

    def __init__(
        self,
        config: Dict[str, Union[str, Dict, List]],
        process: str,
        for_DRtoSR: bool,
        logger: str,
    ):
        """
        Initiating a new evaluator for fake factors using correctionlib.

        Args:
            config: A dictionary with all the relevant information for the fake factor calculation
            process: Name of the process the fake factors were calculated for
            for_DRtoSR: If True fake factors calculated specifically for the DR to SR correction will be loaded, if False the general fake factors will be used
            logger: Name of the logger that should be used
        """
        log = logging.getLogger(logger)

        if not for_DRtoSR:
            self.for_DRtoSR = ""
            self.ff_path = os.path.join(
                "workdir",
                config["workdir_name"],
                config["era"],
                f"fake_factors_{config['channel']}.json",
            )
        else:
            self.for_DRtoSR = "for_DRtoSR"
            self.ff_path = os.path.join(
                "workdir",
                config["workdir_name"],
                config["era"],
                f"corrections/{config['channel']}/fake_factors_{config['channel']}_for_corrections.json",
            )

        self.process = process

        correctionlib.register_pyroot_binding()
        log.info(f"Loading fake factor file {self.ff_path} for process {self.process}")
        ROOT.gInterpreter.Declare(
            f'auto {self.process}_{self.for_DRtoSR} = correction::CorrectionSet::from_file("{self.ff_path}")->at("{self.process}_fake_factors");'
        )

    def evaluate_subleading_lep_pt_njets(self, rdf: Any) -> Any:
        """
        Evaluating the fake factors based on the variables it depends on.
        In this function it is the transverse momentum of the subleading lepton in the tau pair and the number of jets.

        Args:
            rdf: root DataFrame object

        Return:
            root DataFrame object with a new column with the evaluated fake factors
        """
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{pt_2, (float)njets, "nominal"}})',
        )
        return rdf

    def evaluate_leading_lep_pt_njets(self, rdf: Any) -> Any:
        """
        Evaluating the fake factors based on the variables it depends on.
        In this function it is the transverse momentum of the leading lepton in the tau pair and the number of jets.

        Args:
            rdf: root DataFrame object

        Return:
            root DataFrame object with a new column with the evaluated fake factors
        """
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{pt_1, (float)njets, "nominal"}})',
        )
        return rdf

    def evaluate_subleading_lep_pt_njets_deltaR(self, rdf: Any) -> Any:
        """
        Evaluating the fake factors based on the variables it depends on.
        In this function it is the transverse momentum of the subleading lepton in the tau pair,
        the number of jets, and the deltaR distance between the leptons in the tau pair.

        Args:
            rdf: root DataFrame object

        Return:
            root DataFrame object with a new column with the evaluated fake factors
        """
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{pt_2, (float)njets, deltaR_ditaupair, "nominal"}})',
        )
        return rdf


class FakeFactorCorrectionEvaluator:
    """
    Evaluator class to initiate a fake factor correction setup. The fake factor corrections are loaded from an already produced correctionlib file.
    Currently only a readout for a correction dependent on the leading or subleading pt is implemented.
    """

    def __init__(
        self,
        config: Dict[str, Union[str, Dict, List]],
        process: str,
        for_DRtoSR: bool,
        logger: str,
    ):
        """
        Initiating a new evaluator for fake factor corrections using correctionlib.

        Args:
            config: A dictionary with all the relevant information for the fake factor correction calculation
            process: Name of the process the fake factor corrections were calculated for
            for_DRtoSR: If True fake factor corrections calculated specifically for the DR to SR correction will be loaded, if False the general fake factor corrections will be used
            logger: Name of the logger that should be used
        """
        log = logging.getLogger(logger)

        if not for_DRtoSR:
            self.for_DRtoSR = ""
            self.corr_path = os.path.join(
                "workdir",
                config["workdir_name"],
                config["era"],
                f"FF_corrections_{config['channel']}.json",
            )
        else:
            self.for_DRtoSR = "for_DRtoSR"
            self.corr_path = os.path.join(
                "workdir",
                config["workdir_name"],
                config["era"],
                f"corrections/{config['channel']}/FF_corrections_{config['channel']}_{self.for_DRtoSR}.json",
            )

        self.process = process

        correctionlib.register_pyroot_binding()
        log.info(
            f"Loading fake factor correction file {self.corr_path} for process {self.process}"
        )
        if self.process == "QCD" and config["channel"] == "tt":
            ROOT.gInterpreter.Declare(
                f'auto {self.process}_corr_{self.for_DRtoSR} = correction::CorrectionSet::from_file("{self.corr_path}")->at("{self.process}_non_closure_subleading_lep_pt_correction");'
            )
        else:
            ROOT.gInterpreter.Declare(
                f'auto {self.process}_corr_{self.for_DRtoSR} = correction::CorrectionSet::from_file("{self.corr_path}")->at("{self.process}_non_closure_leading_lep_pt_correction");'
            )

    def evaluate_leading_lep_pt(self, rdf: Any) -> Any:
        """
        Evaluating the fake factor corrections based on the variables it depends on.
        In this function it is the transverse momentum of the leading lepton in the tau pair.

        Args:
            rdf: root DataFrame object

        Return:
            root DataFrame object with a new column with the evaluated fake factor corrections
        """
        rdf = rdf.Define(
            f"{self.process}_ff_corr",
            f'{self.process}_corr_{self.for_DRtoSR}->evaluate({{pt_1, "nominal"}})',
        )
        return rdf

    def evaluate_subleading_lep_pt(self, rdf: Any) -> Any:
        """
        Evaluating the fake factor corrections based on the variables it depends on.
        In this function it is the transverse momentum of the subleading lepton in the tau pair.

        Args:
            rdf: root DataFrame object

        Return:
            root DataFrame object with a new column with the evaluated fake factor corrections
        """
        rdf = rdf.Define(
            f"{self.process}_ff_corr",
            f'{self.process}_corr_{self.for_DRtoSR}->evaluate({{pt_2, "nominal"}})',
        )
        return rdf
