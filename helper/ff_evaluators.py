import logging
import os
import ctypes
from typing import Any, Dict, List, Tuple, Union
import traceback

import correctionlib
import correctionlib.schemav2 as cs
import ROOT

try:
    import onnxruntime as ort
    import numpy as np
except ImportError:
    ort = None
    np = None

import helper.functions as func
import CustomLogging as logging_helper


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

        return cls(process, var_dependences, for_DRtoSR, logger=logger)

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

        return cls(process, var_dependences, for_DRtoSR, logger=logger)

    def __init__(
        self,
        process: str,
        var_dependences: List[str],
        for_DRtoSR: bool,
        logger: Union[str, logging.Logger, None] = None,
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

        self.log = logging.getLogger(logger) if logger else logging.getLogger(__name__)

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
        self.log.debug(f"Evaluating fake factor for process {self.process} with variables {self.var_dependences} and for_DRtoSR={self.for_DRtoSR}")
        eval_str = self.str_var_dependences + ', "nominal"'
        rdf = rdf.Define(
            f"{self.process}_fake_factor",
            f'{self.process}_{self.for_DRtoSR}->evaluate({{{eval_str}}})',
        )
        self.log.debug(f"Defined column '{self.process}_fake_factor' with evaluation string: {self.process}_{self.for_DRtoSR}->evaluate({{{eval_str}}})")
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

        variable = corr_variable if isinstance(corr_variable, str) else corr_variable[0]

        correctionlib.register_pyroot_binding()
        log.info(f"Loading fake factor correction file {path} for process {process}")
        ROOT.gInterpreter.Declare(
            f'auto {process}_corr_{variable}_{_for_DRtoSR} = '
            f'correction::CorrectionSet::from_file("{path}")'
            f'->at("{process}_non_closure_{variable}_correction");'
        )

        return cls(process, corr_variable, for_DRtoSR, logger=logger)

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

        variable = corr_variable if isinstance(corr_variable, str) else corr_variable[0]
        literal = correction.json().replace('"', r'\"')

        log.info(f"Loading fake factor correction from string for process {process}")
        correctionlib.register_pyroot_binding()
        ROOT.gInterpreter.Declare(
            f'auto {process}_corr_{variable}_{"for_DRtoSR" if for_DRtoSR else ""} = '
            f'correction::CorrectionSet::from_string("{literal}")'
            f'->at("{process}_non_closure_{variable}_correction");'
        )

        return cls(process, corr_variable, for_DRtoSR, logger=logger)

    def __init__(
        self,
        process: str,
        corr_variable: Union[str, Tuple[str, ...]],
        for_DRtoSR: bool,
        logger: Union[str, logging.Logger, None] = None,
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
        if not isinstance(corr_variable, str) and isinstance(corr_variable[0], str):
            self.variable = corr_variable[0]
            self.var_dependences = corr_variable
        else:
            self.variable = corr_variable
            self.var_dependences = [corr_variable]

        self.log = logging.getLogger(logger) if logger else logging.getLogger(__name__)

    @property
    def corr_str(self) -> str:
        return f"{self.process}_ff_corr_{self.variable}"

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
        self.log.debug(f"Evaluating fake factor correction for process {self.process} with variable {self.variable} and for_DRtoSR={self.for_DRtoSR}")
        eval_str = self.str_var_dependences + ', "nominal"'
        rdf = rdf.Define(
            self.corr_str,
            f"{self.process}_corr_{self.variable}_{self.for_DRtoSR}->evaluate({{{eval_str}}})",
        )
        self.log.debug(f"Defined column '{self.corr_str}' with evaluation string: {self.process}_corr_{self.variable}_{self.for_DRtoSR}->evaluate({{{eval_str}}})")
        return rdf


class DRSRCorrectionEvaluator:
    """
    Dedicated evaluator class for DR to SR corrections.
    Loads a correction named '{process}_DR_SR_correction' from a correctionlib file.
    """

    @classmethod
    def loading_from_file(
        cls,
        config: Dict[str, Union[str, Dict, List]],
        process: str,
        corr_variable: Union[str, Tuple[str, ...]],
        logger: str,
    ) -> "DRSRCorrectionEvaluator":
        """
        Loads a DR_SR correction from a correctionlib file.

        Args:
            config: A dictionary with all the relevant information for the fake factor correction calculation
            process: Name of the process the DR to SR correction was calculated for
            corr_variable: Name of the variable dependence of the correction
            logger: Name of the logger that should be used

        Returns:
            An instance of the DRSRCorrectionEvaluator class.
        """
        log = logging.getLogger(logger)

        directories = ["workdir", config["workdir_name"], config["era"]]
        path = os.path.join(*directories, f"FF_corrections_{config['channel']}.json")

        correctionlib.register_pyroot_binding()
        log.info(f"Loading DR_SR correction from file {path} for process {process}")

        ROOT.gInterpreter.Declare(
            f'auto {process}_corr_DR_SR = '
            f'correction::CorrectionSet::from_file("{path}")'
            f'->at("{process}_DR_SR_correction");'
        )

        return cls(process, corr_variable, logger=logger)

    @classmethod
    def loading_from_CorrectionSet(
        cls,
        correction: cs.CorrectionSet,
        process: str,
        corr_variable: Union[str, Tuple[str, ...]],
        logger: str,
    ) -> "DRSRCorrectionEvaluator":
        """
        Loads a DR_SR correction from a correctionlib.CorrectionSet object.

        Args:
            correction: A correctionlib.CorrectionSet object containing the DR to SR correction.
            process: Name of the process the DR to SR correction was calculated for.
            corr_variable: Name of the variable dependence of the correction.
            logger: Name of the logger that should be used.

        Returns:
            An instance of the DRSRCorrectionEvaluator class.
        """
        log = logging.getLogger(logger)

        assert isinstance(correction, cs.CorrectionSet), "Correction must be of type correctionlib.schemav2.CorrectionSet"

        literal = correction.json().replace('"', r'\"')

        log.info(f"Loading DR_SR correction from string for process {process}")
        correctionlib.register_pyroot_binding()

        ROOT.gInterpreter.Declare(
            f'auto {process}_corr_DR_SR = '
            f'correction::CorrectionSet::from_string("{literal}")'
            f'->at("{process}_DR_SR_correction");'
        )

        return cls(process, corr_variable, logger=logger)

    def __init__(
        self,
        process: str,
        corr_variable: Union[str, Tuple[str, ...]],
        logger: Union[str, logging.Logger, None] = None
    ):
        """
        Initializes the DR to SR correction evaluator.

        Args:
            process: Name of the process the DR to SR correction was calculated for.
            corr_variable: Name of the variable dependence of the correction.

        Returns:
            None
        """
        self.process = process
        if not isinstance(corr_variable, str) and isinstance(corr_variable[0], str):
            self.variable = corr_variable[0]
            self.var_dependences = corr_variable
        else:
            self.variable = corr_variable
            self.var_dependences = [corr_variable]

        self.log = logging.getLogger(logger) if logger else logging.getLogger(__name__)

    @property
    def corr_str(self) -> str:
        """
        Returns the unique column name for the applied correction.
        """
        return f"{self.process}_DR_SR_correction"

    @property
    def str_var_dependences(self) -> List[str]:
        """
        Returns a list of strings with the variable dependences for the correction.
        """
        return [f"(float){var}" for var in self.var_dependences]

    def evaluate_correction(self, rdf: Any) -> Any:
        """
        Applies the DR to SR correction to the RDataFrame.

        Args:
            rdf: root DataFrame object

        Returns:
            root DataFrame object with a new column with the evaluated DR to SR corrections
        """
        self.log.debug(f"Evaluating DR to SR correction for process {self.process} with variable {self.variable}")
        eval_str = ", ".join(self.str_var_dependences) + ', "nominal"'
        root_corr_name = f"{self.process}_corr_DR_SR"
        rdf = rdf.Define(
            self.corr_str,
            f"{root_corr_name}->evaluate({{{eval_str}}})",
        )
        self.log.debug(f"Defined column '{self.corr_str}' with evaluation string: {root_corr_name}->evaluate({{{eval_str}}})")
        return rdf


def _onnx_ctypes_callback(proc_id: int, features_ptr) -> float:
    """
    Pure Python function called natively from C++.
    Intercepts raw memory pointer of the RDataFrame array, zero-copy casts  it to numpy, and
    runs ONNX inference.
    """
    try:
        proc_info = _PROCESS_ID_MAP[proc_id]
        sess = proc_info["sess"]
        ort_input_names = proc_info["ort_input_names"]
        ort_input_shapes = proc_info["ort_input_shapes"]
        n_features = proc_info["n_features"]

        features_array = np.ctypeslib.as_array(features_ptr, shape=(n_features,))

        if len(ort_input_names) == 1:
            expected_shape = ort_input_shapes[0]

            target_shape = [1 if not isinstance(dim, int) else dim for dim in expected_shape]  # Potential dynamic batch dim -> 1.

            tensor = features_array.astype(np.float32).reshape(target_shape)
            inputs = {ort_input_names[0]: tensor}
        else:
            inputs = {}
            for name, shape, val in zip(ort_input_names, ort_input_shapes, features_array):
                target_shape = [1 if not isinstance(dim, int) else dim for dim in shape]
                inputs[name] = np.array([val], dtype=np.float32).reshape(target_shape)

        out = sess.run(None, inputs)

        return float(np.asarray(out[0]).item())
    except Exception:
        traceback.print_exc()
        return -999.0


# Global reference to the ctypes callback to avoid Python's Garbage Collector
_PROCESS_ID_MAP = {}
_C_CALLBACK_TYPE = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_int, ctypes.POINTER(ctypes.c_double))

_GLOBAL_C_CALLBACK = _C_CALLBACK_TYPE(_onnx_ctypes_callback)
_GLOBAL_C_CALLBACK_PTR = ctypes.cast(_GLOBAL_C_CALLBACK, ctypes.c_void_p).value


class ONNXFakeFactorEvaluator:
    """
    Evaluator class to initiate a fake factor setup utilizing ONNX models.
    The fake factors are evaluated dynamically inside the RDataFrame via ONNXRuntime.
    """

    @classmethod
    def loading_from_config(
        cls,
        nn_config: Dict[str, Any],
        process: str,
        logger: str,
    ) -> "ONNXFakeFactorEvaluator":
        log = logging.getLogger(logger)
        if process not in nn_config.get("target_processes", {}):
            raise KeyError(f"Process {process} not found in NN config.")

        proc_config = nn_config["target_processes"][process]
        model_path = proc_config["model_path"]
        model_inputs = proc_config["model_input"]
        define_cols = proc_config.get("define_columns", {})

        log.info(f"Loading ONNX model from {model_path} for process {process}")
        return cls(process, model_path, model_inputs, define_cols, logger)

    def __init__(
        self,
        process: str,
        model_path: str,
        model_inputs: List[str],
        define_columns: Dict[str, str],
        logger: Union[str, logging.Logger, None] = None
    ):
        if ort is None:
            raise ImportError("onnxruntime and numpy are required for ONNXFakeFactorEvaluator")

        self.process = process
        self.model_path = model_path
        self.model_inputs = model_inputs
        self.define_columns = define_columns

        self.log = logging.getLogger(logger) if logger else logging.getLogger(__name__)
        self._is_initialized = False

    def __getstate__(self):
        """
        Called when passing this object to a Multiprocessing Worker.
        """
        state = self.__dict__.copy()
        state['_is_initialized'] = False
        if 'log' in state:
            del state['log']
        return state

    def __setstate__(self, state):
        """
        Called inside the Multiprocessing Worker to rebuild the object.
        """
        self.__dict__.update(state)
        self.log = logging.getLogger(self.log_name)

    def _initialize_worker_state(self):
        """
        Lazy Initialization: Sets up ONNX and ROOT C++ injection the first
        time it is needed in whatever process calls it.
        """
        if self._is_initialized:
            return

        self.log.debug(f"Initializing ONNX Session and C++ context for {self.process} in PID {os.getpid()}")

        if self.process not in [info.get("process") for info in _PROCESS_ID_MAP.values()]:
            self.proc_id = len(_PROCESS_ID_MAP)
            sess = ort.InferenceSession(self.model_path)
            _PROCESS_ID_MAP[self.proc_id] = {
                "process": self.process,
                "sess": sess,
                "ort_input_names": [inp.name for inp in sess.get_inputs()],
                "ort_input_shapes": [inp.shape for inp in sess.get_inputs()],
                "n_features": len(self.model_inputs)
            }
        else:
            self.proc_id = next(k for k, v in _PROCESS_ID_MAP.items() if v["process"] == self.process)

        if not hasattr(ROOT, self.cxx_func_name):  # Inject C++ Code into ROOT's Cling compiler just once
            cpp_code = f"""
            #ifndef ONNX_EVAL_DEFINED_{self.process}
            #define ONNX_EVAL_DEFINED_{self.process}
            double {self.cxx_func_name}(const ROOT::VecOps::RVec<double>& features) {{
                typedef double (*CallbackType)(int, const double*);
                CallbackType cb = (CallbackType){_GLOBAL_C_CALLBACK_PTR}ULL;
                return cb({self.proc_id}, features.data());
            }}
            #endif
            """
            self.log.debug(f"Registering C++ execution wrapper for {self.process}")
            ROOT.gInterpreter.Declare(cpp_code)

        self._is_initialized = True

    def evaluate_fake_factor(self, rdf: Any) -> Any:
        """
        Evaluating the fake factors using the loaded ONNX model.
        Missing columns are defined on the fly before inference.
        """
        self._initialize_worker_state()

        existing_cols = [str(_column) for _column in rdf.GetColumnNames()]

        for column, expression in self.define_columns.items():
            if column not in existing_cols:
                self.log.debug(f"Defining column '{column}' as '{expression}'")
                rdf = rdf.Define(column, expression)

        rvec_col = f"{self.process}_onnx_features"
        inputs_csv = ", ".join([f"(double){inp}" for inp in self.model_inputs])
        rvec_expr = f"ROOT::VecOps::RVec<double>{{{inputs_csv}}}"
        self.log.debug(f"Packaging {len(self.model_inputs)} inputs into RVec column '{rvec_col}'")
        rdf = rdf.Define(rvec_col, rvec_expr)

        self.log.info(f"Applying ONNX fake factor evaluation for {self.process}")
        rdf = rdf.Define(f"{self.process}_fake_factor", f"{self.cxx_func_name}({rvec_col})")

        return rdf


def get_fake_factor_evaluator(
    config: Dict[str, Any],
    process: str,
    var_dependences: List[str],
    for_DRtoSR: bool,
    logger: str,
) -> Union[FakeFactorEvaluator, ONNXFakeFactorEvaluator]:
    """
    Factory function to decide whether to load the traditional correctionlib JSON
    evaluator or the new ONNX YAML-based NN evaluator.
    """

    if process in config.get("target_processes", {}) and "model_path" in config["target_processes"][process]:
        return ONNXFakeFactorEvaluator.loading_from_config(config, process, logger)

    directories = ["workdir", config["workdir_name"], config["era"]]
    if not for_DRtoSR:
        yaml_path = os.path.join(*directories, f"fake_factors_models_{config['channel']}.yaml")
    else:
        yaml_path = os.path.join(*directories, f"fake_factors_models_DR_SR_{config['channel']}.yaml")

    # with open(yaml_path, "r") as f:
    #     nn_config = func.configured_yaml.load(f)

    # return ONNXFakeFactorEvaluator.loading_from_config(nn_config, process, logger)

    return FakeFactorEvaluator.loading_from_file(config, process, var_dependences, for_DRtoSR, logger)

    # if os.path.exists(yaml_path):
    #     with open(yaml_path, "r") as f:
    #         nn_config = func.configured_yaml.load(f)
    #     if process in nn_config.get("target_processes", {}) and "model_path" in nn_config["target_processes"][process]:
    #         return ONNXFakeFactorEvaluator.loading_from_config(nn_config, process, logger)

    # return FakeFactorEvaluator.loading_from_file(config, process, var_dependences, for_DRtoSR, logger)
