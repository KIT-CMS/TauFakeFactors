"""
Collection of helpful functions for other scripts
"""

import concurrent.futures
import glob
import hashlib
import inspect
import json
import logging
import os
import re
import sys
from decimal import Decimal
from typing import Any, Callable, Dict, List, Tuple, Union

import numpy as np
import ROOT
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.scalarstring import DoubleQuotedScalarString
from XRootD import client


THIS_DIR = os.path.dirname(os.path.abspath(__file__))
TAU_FAKE_FACTORS_DIR = os.path.dirname(THIS_DIR)


class CachingKeyHelper:
    @staticmethod
    def make_hashable(obj: Union[Dict, List, Tuple, Any]) -> Union[Dict, Tuple, bytes, Any]:
        """
        Recursively convert an object to a hashable representation.

        Args:
            obj: The object to convert.

        Returns:
            A hashable representation of the object.
        """
        if isinstance(obj, dict):
            return {k: CachingKeyHelper.make_hashable(v) for k, v in sorted(obj.items())}
        elif isinstance(obj, (list, tuple)):
            return tuple(CachingKeyHelper.make_hashable(x) for x in obj)
        elif isinstance(obj, np.ndarray):
            return obj.tobytes()
        else:
            return obj

    @staticmethod
    def generate_key(args: Tuple[Any, ...], kwargs: Dict[str, Any]) -> str:
        """
        Generate a stable hash key for a function's arguments and keyword arguments.

        This key can be used to identify cached results uniquely.

        Args:
            args: Positional arguments passed to the function.
            kwargs: Keyword arguments passed to the function.

        Returns:
            A hash key as a string.
        """
        stable_representation = {"args": args, "kwargs": kwargs}
        key_string = json.dumps(stable_representation, sort_keys=True, separators=(',', ':'))
        return hashlib.md5(key_string.encode()).hexdigest()

    @staticmethod
    def get_assignment_variable(stack_depth: int = 2) -> str:
        """
        Try to extract the variable name on the left hand sign of the assignment where this function was
        called. Works for single- and multi-line assignments, useful to determine caching placement.

        Args:
            stack_depth (int): The depth of the stack to inspect. Defaults to 2.

        Returns:
            str: The name of the variable on the left hand side of the assignment, or "unknown" if it cannot be determined.
        """
        try:
            frame = inspect.stack()[stack_depth]
            lines, start = inspect.getsourcelines(frame.frame)
            call_lineno = frame.lineno - start  # Find line with function call
            src_lines = []
            for line in lines[call_lineno:]:
                src_lines.append(line)
                if ")" in line:   # look ahead until the closing parenthesis
                    break

            src = "".join(src_lines)
            src = re.sub(r"\s+", " ", src.strip())

            match = re.match(r"(\w+)\s*=\s*.*apply_region_filters", src)  # Regex for variable assignment
            if match:
                return match.group(1)
        except Exception:
            pass
        return "unknown"


def to_commented_map(items: dict) -> CommentedMap:
    """
    Recursively converts a Python dictionary to a ruamel.yaml CommentedMap,
    ensuring that list leaves are in flow style and keys with operators are quoted.
    This function is now centralized here.
    """
    if not isinstance(items, dict):
        return items
    commented_map = CommentedMap()
    for key, value in items.items():
        if any(op in str(key) for op in ["<=", ">=", "<", ">", "==", "#"]):
            key = DoubleQuotedScalarString(key)

        if isinstance(value, dict):
            commented_map[key] = to_commented_map(value)
        elif isinstance(value, list) and not any(isinstance(i, (dict, list)) for i in value):
            cs = CommentedSeq(value)
            cs.fa.set_flow_style()
            commented_map[key] = cs
        elif isinstance(value, list):
            commented_map[key] = [to_commented_map(item) for item in value]
        else:
            commented_map[key] = value
    return commented_map


class ConfiguredYAML(YAML):
    """
    Pre-configured YAML for handling specific formatting setting:
    - Better float representer to avoid scientific notation (those are stored as strings).
    - Improved dump() method for automatic apply of to_commented_map for better readable list structures.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(typ="rt", *args, **kwargs)

        self.preserve_quotes = True
        self.indent(mapping=2, sequence=4, offset=2)
        self.width = 4096  # No line break in the output, never!

        def smart_float_representer(representer, data):
            # Handle special float values inf and -inf before parsing the Decimal, which would throw
            # an error message.
            if data == float('inf'):
                return representer.represent_scalar('tag:yaml.org,2002:float', '.inf')
            if data == float('-inf'):
                return representer.represent_scalar('tag:yaml.org,2002:float', '-.inf')

            # At this point, data is expected to be a regular float, here we can go on with the
            # cleaning.
            cleaned_decimal = Decimal(str(data)).quantize(Decimal('1e-12'))
            formatted_str = f"{cleaned_decimal:.15f}".rstrip('0')
            if formatted_str.endswith('.'):
                formatted_str += '0'
            return representer.represent_scalar('tag:yaml.org,2002:float', formatted_str)

        self.representer.add_representer(float, smart_float_representer)
        self.representer.add_representer(np.floating, smart_float_representer)

        merge_tag = 'tag:yaml.org,2002:merge'
        if merge_tag in self.constructor.yaml_constructors:
            self.constructor.yaml_constructors[merge_tag].deep = True


configured_yaml = ConfiguredYAML()


class RuntimeVariables(object):
    """
    A singleton-like container class holding variables that can be adjusted at runtime.

    Attributes:
        USE_MULTIPROCESSING (bool): Flag to enable or disable multiprocessing globally
    """
    USE_MULTIPROCESSING = True
    USE_CACHED_INTERMEDIATE_STEPS = False
    RDataFrameWrapper = None

    def __new__(cls) -> "RuntimeVariables":
        if not hasattr(cls, "instance"):
            cls.instance = super(RuntimeVariables, cls).__new__(cls)
            return cls.instance


def get_cached_file_path(
    output_path: str,
    process: Union[str, None] = None,
    variables: Union[List[Union[Tuple[str, ...], str]], None] = None,
    for_DRtoSR: bool = False,
) -> str:
    """
    Function to get the path of a cached file.
    The cached file is stored in the ".cache" folder in the output path, the file name
    is generated based on the process and the variables if it is given.
    If the ".cache" folder does not exist, it is created.

    In case of DR_SR fake factors the name is generated without the process.
    In case of DR_SR correction the name is generated without the variables.
    The file name is generated in the following format:
    <corr_type>_<for_DRtoSR>_<process>_<variables>.pickle
    where <corr_type> is either "_DR_SR_" or "_non_closure" depending on the variables.
    The <for_DRtoSR> is only added if the for_DRtoSR argument is set to True.
    The <process> is the name of the process the file correspond to.
    The <variables> are the variables which are used for the non-closure correction, joined by an underscore if given.

    Args:
        output_path: Path to the folder where the file should be stored at
        process: Name of the process the file correspond to
        variables: List of variables which are used for the non-closure correction
        for_DRtoSR: If True, the cached path is generated for the DR to SR correction
    Return:
        String with the file name
    """
    cache_path = os.path.join(output_path, ".cache")
    os.makedirs(cache_path, exist_ok=True)

    corr_type_str = "DR_SR" if variables is None else "non_closure"
    for_DRtoSR_str = "for_DRtoSR" if for_DRtoSR else ""
    process_str = "" if process is None else process
    if variables is None:
        variables_str = ""
    elif isinstance(variables[0], str):
        variables_str = "_".join(variables)
    elif isinstance(variables[0], Tuple):
        variables_str = "_".join("_".join(it) for it in variables)
    else:
        raise TypeError(
            f"""
                Unsupported type: {type(variables)}.
                You sure you are doing the right thing?
            """
        )
    file_name = hashlib.md5(
        "_".join(
            [
                corr_type_str,
                for_DRtoSR_str,
                process_str,
                variables_str
            ]
        ).encode()
    ).hexdigest() + ".pickle"
    return os.path.join(cache_path, file_name)


def correction_config_comparison(
    test_config: dict,
    config: dict,
    *,
    process: str,
    closure_corr: str,
    for_DRtoSR: bool = False,
) -> bool:
    """
    Function to compare two configurations for a specific process non-closure correction.
    It will compare the non-closure correction configuration for the specified process,
    in case of non-closure for DRtoSR correction it will also compare the SRlike and
    ARlike cuts.

    Args:
        test_config: The configuration to be tested for equality
        config: The reference configuration
        process: The process to be compared
        closure_corr: The non-closure correction to be compared
    Returns:
        bool: True if the configurations are equal, False otherwise
    """

    _test_config = test_config["target_processes"][process]
    _config = config["target_processes"][process]

    is_same = True
    if for_DRtoSR:
        if "DR_SR" not in _test_config or "DR_SR" not in _config:
            return False

        _test_config, _config = _config["DR_SR"], _test_config["DR_SR"]

        is_same &= nested_object_comparison(_test_config["SRlike_cuts"], _config["SRlike_cuts"])
        is_same &= nested_object_comparison(_test_config["ARlike_cuts"], _config["ARlike_cuts"])

    _test_config = _test_config["non_closure"][closure_corr]
    _config = _config["non_closure"][closure_corr]

    return is_same and nested_object_comparison(_test_config, _config)


def nested_object_comparison(obj_1: Any, obj_2: Any) -> bool:
    """
    Function to compare two objects of potentially different types.
    It will return True if they are equal, and False otherwise.

    The function will compare:
    - int, float, str: by value
    - list, tuple: by length and value (recursively)
    - dict: by length and value (recursively)
    - np.ndarray: by shape and value (using np.array_equal)
    - Any other type: raises TypeError

    Args:
        obj_1: The first object to compare.
        obj_2: The second object to compare.

    Returns:
        bool: True if the objects are equal, False otherwise.
    """
    if type(obj_1) != type(obj_2):
        return False
    if obj_1 is None and obj_2 is None:
        return True
    elif isinstance(obj_1, (int, float, str)):
        return obj_1 == obj_2
    elif isinstance(obj_1, (list, tuple)):
        return len(obj_1) == len(obj_2) and all(nested_object_comparison(a, b) for a, b in zip(obj_1, obj_2))
    elif isinstance(obj_1, dict):
        return set(obj_1.keys()) == set(obj_2.keys()) and all(nested_object_comparison(obj_1[k], obj_2[k]) for k in obj_1)
    elif isinstance(obj_1, np.ndarray):
        return obj_1.shape == obj_2.shape and np.array_equal(obj_1, obj_2)
    else:
        raise TypeError(
            f"""
                Unsupported type: {type(obj_1)}.
                Add more logic to handle this type or do you really want to compare it?
            """
        )


def remove_empty_keys(data: Union[Dict, List, Any]) -> Union[Dict, List, Any]:
    """
    Recursively remove keys from a dictionary whose values are:
        * empty lists
        * empty dictionaries
        * empty numpy arrays.
    arrays/lists with at least one element are retained.

    Args:
        data (dict, list, or any): Input data to filter.

    Returns:
        Filtered data without empty lists, dicts, numpy arrays.
    """
    if isinstance(data, dict):
        new_dict = {}
        for key, value in data.items():
            cleaned_value = remove_empty_keys(value)
            if isinstance(cleaned_value, np.ndarray) and cleaned_value.size > 0:
                new_dict[key] = cleaned_value
            elif isinstance(cleaned_value, (list, dict)) and cleaned_value:
                new_dict[key] = cleaned_value
        return new_dict
    elif isinstance(data, list):
        new_list = [remove_empty_keys(item) for item in data]
        filtered = []
        for item in new_list:
            if isinstance(item, np.ndarray) and item.size > 0:
                filtered.append(item)
            elif isinstance(item, (list, dict)) and item:
                filtered.append(item)
        return filtered
    else:
        return data


def optional_process_pool(
    args_list: List[Tuple[Any, ...]],
    function: Callable,
    max_workers: Union[int, None] = None,
) -> List[Any]:
    """
    Running a function with a list of arguments in parallel using multiprocessing if
    the list of arguments is longer than one and multiprocessing is enabled.

    Args:
        args_list: List of tuples with arguments for the function
        function: Function to be executed
        max_workers: Number of workers to be used in the multiprocessing pool (default: None)

    Return:
        List of results of the function

    """

    if len(args_list) == 1 or not RuntimeVariables.USE_MULTIPROCESSING:
        results = [function(args) for args in args_list]
    else:
        n = max_workers if max_workers is not None else len(args_list)
        with concurrent.futures.ProcessPoolExecutor(max_workers=n) as executor:
            results = list(executor.map(function, args_list))

    return results


def load_config(config_file: str) -> Dict:
    """
    This function loads the configuration file. It first loads the common settings which
    should be present in the same folder as the config file.
    The common settings are overwritten if they are also specified in the config file.

    Args:
        config_file: Path to the specific config file

    Return:
        Configuration as a dictionary
    """
    common_config_file_name = "common_settings.yaml"

    common_config_file = os.path.join(
        os.path.split(config_file)[0],
        common_config_file_name,
    )
    if os.path.exists(common_config_file):
        print(
            f"Using {common_config_file_name} file found in {os.path.split(config_file)[0]}"
        )
        print("-" * 50)
    else:
        print("No common config file found!")

    # Container of the loaded configuration
    # 
    # Some default values are pre-defined in the config dict that is going to contain the loaded
    # configuration. These values are overwritten if they are explicitly set in the common config file.
    #
    # The variables, for which defaults are set, are:
    #
    # - 'sample_database`: Path to the sample database directory. Usuallly, this path is set to the
    #   `datasets` submodule of the `TauFakeFactors` module. Users can set a custom path, e.g.,
    #   to an external path to a working version of their sample database.
    config = {
        "sample_database": os.path.join(TAU_FAKE_FACTORS_DIR, "datasets"),
    }

    # Update the config with common settings, applying to all steps
    with open(common_config_file, "r") as file:
        config.update(configured_yaml.load(file))

    # loading of the chosen config file
    try:
        with open(config_file, "r") as file:
            _config = configured_yaml.load(file)
            repeating_keys = set(config.keys()).intersection(set(_config.keys()))
            if repeating_keys:
                _overwriting = "\n\t".join(f"{k}: {_config[k]}" for k in repeating_keys)
                print(
                    f"""
    Warning: The following keys are present in both the common and the specific config file:
        {repeating_keys}
    and are overwritten by the ones in {config_file} to
        {_overwriting}
                    """
                )
            config.update(_config)

    except FileNotFoundError:
        print(f"Error: Config file {config_file} not found.")
        sys.exit(1)

    return config


def check_path(path: str) -> None:
    """
    This function checks if a given path exist. If not, this path is created.

    Args:
        path: path to check, as a string

    Return:
        None
    """
    if not os.path.exists(path):
        os.makedirs(path, exist_ok=True)


def setup_logger(
    log_file: str, log_name: str, log_level: int, subcategories: Union[List[str], None] = None
) -> None:
    """
    Setting up all relevant loggers and handlers.

    Args:
        log_file: Name of the file the logging information will be stored in
        log_name: General name of the logger
        log_level: Level of the logger, e.g. logging.INFO, logging.DEBUG, etc.
        subcategories: List of different sub logger names e.g. can be used to differentiate between processes (default: None)

    Return:
        None
    """
    # create file handler
    fh = logging.FileHandler(log_file)
    fh.setLevel(log_level)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(log_level)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)

    if subcategories is not None:
        for cat in subcategories:
            log = logging.getLogger(f"{log_name}.{cat}")
            log.setLevel(log_level)
            # add the handlers to logger
            log.addHandler(ch)
            log.addHandler(fh)
    else:
        log = logging.getLogger(f"{log_name}")
        log.setLevel(log_level)
        # add the handlers to logger
        log.addHandler(ch)
        log.addHandler(fh)


def get_ntuples(config: Dict, process: str, sample: str) -> List[str]:
    """
    This function generates a list of paths of all ntuples for a specific sample of a process.

    Args:
        config: Dictionary with the configuration information for the preselection
        process: General name of a process e.g. "ttbar"
        sample: Exact name of the folder where the ntuples are stored e.g. "SingleMuon_Run2018A-UL2018"

    Return:
        List of file paths
    """
    log = logging.getLogger(f"preselection.{process}")
    sample_path = os.path.join(
        config["ntuple_path"], config["era"], sample, config["channel"]
    )
    log.info(
        f"The following files are loaded for era: {config['era']}, channel: {config['channel']}, sample {sample}"
    )
    # now check, if the files exist
    selected_files = check_inputfiles(
        path=sample_path, process=process, tree=config["tree"]
    )

    log.info("-" * 50)
    for file in selected_files:
        log.debug(file)
    log.debug("-" * 50)

    return selected_files


def check_inputfiles(path: str, process: str, tree: str) -> List[str]:
    """
    Additional function to check if the input files are empty. If yes, they are skipped.

    Args:
        path: Path where files should be checked
        process: General name of a process e.g. "ttbar"
        tree: Name of the tree in the input files

    Return:
        List of file paths with not empty files
    """
    log = logging.getLogger(f"preselection.{process}")

    fsname = "root://cmsdcache-kit-disk.gridka.de/"
    xrdclient = client.FileSystem(fsname)
    status, listing = xrdclient.dirlist(path.replace(fsname, ""))

    if not status.ok:
        log.info(f"Error: {status.message}")
        sys.exit(1)

    selected_files = []
    for f in listing:
        if f.name.endswith(".root"):
            # check if file is empty
            if check_for_empty_tree(file_path=os.path.join(path, f.name), tree=tree):
                log.info(f"File {f.name} is empty. Skipping.")
                continue
            selected_files.append(os.path.join(path, f.name))

    return selected_files


def check_for_empty_tree(file_path: str, tree: str) -> bool:
    """
    Checking if the tree in a file is empty.

    Args:
        file_path: Path of the input file
        tree: Name of the tree in the input file

    Return:
        Boolean which is true if the file is empty, false otherwise.
    """
    f = ROOT.TFile.Open(file_path)
    t = f.Get(tree)
    if not isinstance(t, ROOT.TTree):
        return True

    return bool(t.GetEntries() == 0)


def rdf_is_empty(rdf: ROOT.RDataFrame) -> bool:
    """
    Function to check if a root DataFrame is empty.

    Args:
        rdf: root DataFrame object

    Return:
        Boolean which is true if the root DataFrame is empty, false otherwise.
    """
    try:
        cols = rdf.GetColumnNames()
        if len(cols) == 0:
            return True
    except:
        return True

    return False


def get_output_name(
    path: str, process: str, tau_gen_mode: str
) -> str:
    """
    Function to generate a file output path where the file will be stored.

    Args:
        path: Path to the folder where the file should be stored at
        process: Name of the process the file correspond to
        tau_gen_mode: Specifying the applied tau pair origin selection

    Return:
        String with the file path
    """
    if tau_gen_mode == "all":
        tau_gen_mode = ""
    else:
        tau_gen_mode = "_" + tau_gen_mode

    return os.path.join(path, f"{process}{tau_gen_mode}.root")


def rename_boosted_variables(rdf: Any, channel: str) -> Any:
    """
    Function to redefine variables to the boosted tau pair information. Redefining only variables
    which are written out for the fake factor measurement. Due to the hardcoded naming and redifinitions
    this function needs to be adjusted if something changes in the list of output variables.

    Args:
        rdf: root DataFrame
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"

    Return:
        root DataFrame with redefined variables
    """
    rdf = rdf.Redefine("njets", "njets_boosted")
    rdf = rdf.Redefine("nbtag", "nbtag_boosted")
    rdf = rdf.Redefine("metphi", "metphi_boosted")
    rdf = rdf.Redefine("met", "met_boosted")
    rdf = rdf.Redefine("pt_1", "boosted_pt_1")
    rdf = rdf.Redefine("q_1", "boosted_q_1")
    rdf = rdf.Redefine("pt_2", "boosted_pt_2")
    rdf = rdf.Redefine("q_2", "boosted_q_2")
    rdf = rdf.Redefine("mt_1", "boosted_mt_1")
    rdf = rdf.Redefine("iso_1", "boosted_iso_1")
    rdf = rdf.Redefine("mass_2", "boosted_mass_2")
    rdf = rdf.Redefine("tau_decaymode_2", "boosted_tau_decaymode_2")
    rdf = rdf.Redefine("deltaR_ditaupair", "boosted_deltaR_ditaupair")
    rdf = rdf.Redefine("m_vis", "boosted_m_vis")
    rdf = rdf.Redefine("fj_Xbb_pt", "fj_Xbb_pt_boosted")
    rdf = rdf.Redefine("fj_Xbb_eta", "fj_Xbb_eta_boosted")
    rdf = rdf.Redefine(
        "fj_Xbb_particleNet_XbbvsQCD", "fj_Xbb_particleNet_XbbvsQCD_boosted"
    )
    rdf = rdf.Redefine("bpair_pt_1", "bpair_pt_1_boosted")
    rdf = rdf.Redefine("bpair_pt_2", "bpair_pt_2_boosted")
    rdf = rdf.Redefine("bpair_btag_value_2", "bpair_btag_value_2_boosted")
    rdf = rdf.Redefine("bpair_eta_2", "bpair_eta_2_boosted")

    if "boosted_gen_match_2" in rdf.GetColumnNames():
        rdf = rdf.Redefine("gen_match_2", "boosted_gen_match_2")
    else:
        rdf = rdf.Define("boosted_gen_match_2", "-1.")
        rdf = rdf.Redefine("gen_match_2", "boosted_gen_match_2")

    if "btag_weight_boosted" in rdf.GetColumnNames():
        rdf = rdf.Redefine("btag_weight", "btag_weight_boosted")
    else:
        rdf = rdf.Define("btag_weight_boosted", "1.")
        rdf = rdf.Redefine("btag_weight", "btag_weight_boosted")

    if "pNet_Xbb_weight_boosted" in rdf.GetColumnNames():
        rdf = rdf.Redefine("pNet_Xbb_weight", "pNet_Xbb_weight_boosted")
    else:
        rdf = rdf.Define("pNet_Xbb_weight_boosted", "1.")
        rdf = rdf.Redefine("pNet_Xbb_weight", "pNet_Xbb_weight_boosted")

    if channel == "tt":
        rdf = rdf.Redefine("mass_1", "boosted_mass_1")
        rdf = rdf.Redefine("tau_decaymode_1", "boosted_tau_decaymode_1")
        if "boosted_gen_match_1" in rdf.GetColumnNames():
            rdf = rdf.Redefine("gen_match_1", "boosted_gen_match_1")
        else:
            rdf = rdf.Define("boosted_gen_match_1", "-1.")
            rdf = rdf.Redefine("gen_match_1", "boosted_gen_match_1")

    if channel == "et":
        rdf = rdf.Redefine("extraelec_veto", "boosted_extraelec_veto")
    if channel == "mt":
        rdf = rdf.Redefine("extramuon_veto", "boosted_extramuon_veto")

    return rdf

def define_columns(rdf: Any, column_definitions: dict, process: str) -> Any:
    """
    Customizer function to define additional columns in the ntuples.
     
    The `column_definitions` dictionary is usually provided with the preselection configuration
    file. The keys of the dictionary correspond to the columns to be created. The values are
    dictionaries which contain the information for the column information. The keys of these inner
    dictionaries have the following meaning:

    - `expression`: The expression string which is used to define the new column.

    - `exclude_processes` (_optional_): A list of process names for which the definition should be
      skipped. If `process` is in this list, the definition is not applied.

    Note that the new column names must not exist in the ntuples, otherwise an error is raised.

    Args:
        rdf: root DataFrame
        column_definitions: Dictionary mapping new column names (keys) to expressions (values)
        process: Name of the current process

    Return:
        root DataFrame with redefined variables
    """

    # Ensure that the new column names are not already present in the ntuple
    rdf_columns = set(rdf.GetColumnNames())
    new_columns = set(column_definitions.keys())
    intersection = rdf_columns.intersection(new_columns)
    if intersection:
        raise ValueError(
            f"The following new column names already exist in the ntuple: {intersection}"
        )

    # Perform the define declarations on the RDataFrame object
    for new_column, define_dict in column_definitions.items():
        expression = define_dict["expression"]
        exclude_processes = define_dict.get("exclude_processes", [])
        if process in exclude_processes:
            continue
        rdf = rdf.Define(new_column, expression)

    return rdf

def get_samples(config: Dict[str, Union[str, Dict, List]]) -> List[str]:
    """
    Function to get a list of all sample paths which will be used for the fake factor calculation.
    This function assumes that the preselection step was already finished and takes into account
    if embedded events or MC events should be used.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation

    Return:
        List of all paths to the relevant samples
    """
    log = logging.getLogger("ff_calc")

    general_sample_path = os.path.join(
        config["file_path"], "preselection", config["era"], config["channel"], "*.root"
    )
    log.info(
        f"The following files are loaded for era: {config['era']}, channel: {config['channel']} from {general_sample_path}"
    )
    log.info("-" * 50)
    sample_paths = glob.glob(general_sample_path)
    tmp_list = glob.glob(general_sample_path)

    for f in tmp_list:
        sample = f.rsplit("/")[-1].rsplit(".")[0]
        if config["use_embedding"] and "_T" in sample:
            sample_paths.remove(f)
        elif not config["use_embedding"] and sample == "embedding":
            sample_paths.remove(f)

    for f in sample_paths:
        log.info(f)
    log.info("-" * 50)

    return sample_paths


def check_categories(config: Dict[str, Union[str, Dict, List]]) -> None:
    """
    A function to check if the categories in the configuration file are defined properly.
    Categories are defined by to parameters, the first one are the orthogonal cuts which split the data
    into the categories and the second one is a binning corresponding to this categories.
    The binning is needed to write a correct correctionlib file. The number of categories and bins have to match.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation

    Return:
        None
    """
    def _recursive_check(categories_config, binedges_config, process_name, var_name):
        if isinstance(binedges_config, list):  # binedge_config is a list -> leaf -> break
            if not isinstance(categories_config, list):
                raise TypeError(f"Category definition for {process_name}/{var_name} is not a list, but binedges are.")
            if len(categories_config) != len(binedges_config) - 1:
                raise Exception(
                    f"Categories and binning for the categories does not match up for {process_name}, {var_name}.\n"
                    f"Found {len(categories_config)} categories and {len(binedges_config)} bin edges."
                )
            return None

        if not isinstance(categories_config, dict) or not isinstance(binedges_config, dict):
            return None    # old style config compatible (not nested 1evel dicts)

        if categories_config.keys() != binedges_config.keys():
            raise Exception(
                f"Keys for nested categories and binnings do not match for {process_name}, {var_name}.\n"
                f"Category keys: {categories_config.keys()}\nBinedge keys: {binedges_config.keys()}"
            )

        for key in categories_config:
            _recursive_check(categories_config[key], binedges_config[key], process_name, f"{var_name}/{key}")

    if "target_processes" in config:
        for process, process_conf in config["target_processes"].items():
            if "split_categories" in process_conf and "split_categories_binedges" in process_conf:
                split_cats = process_conf["split_categories"]
                split_bins = process_conf["split_categories_binedges"]

                if split_cats.keys() != split_bins.keys():
                    raise Exception(f"Split variables in split_categories and split_categories_binedges do not match for {process}.")

                for var in split_cats:
                    _recursive_check(split_cats[var], split_bins[var], process, var)

    if "process_fractions" in config:
        fraction_categories = config["process_fractions"]["split_categories"]
        fraction_categories_edges = config["process_fractions"][
            "split_categories_binedges"
        ]
        for cat in fraction_categories:
            if len(fraction_categories[cat]) != (
                len(fraction_categories_edges[cat]) - 1
            ):
                raise Exception("Categories and binning for the categories does not match up for {cat} for fractions.")


def modify_config(
    config: Dict[str, Union[str, Dict, List]],
    corr_config: Dict[str, Union[str, Dict, List]],
    process: str,
    to_AR_SR: bool = False,
) -> None:
    """
    This function modifies or adds cuts in the general configuration based on the cut changes
    defined in the corresponding correction configuration for a specific process.

    Args:
        config: A dictionary with all the relevant information for the fake factor calculation
        corr_config: A dictionary with all the relevant information for calculating corrections to the measured fake factors
        process: Name of the process for which the configuration should be modified
        to_AR_SR: If True the cut configuration is modified to the signal/application region and not the determination region anymore, this change is relevant for the DR to SR correction

    Return:
        None
    """
    if "SRlike_cuts" in corr_config:
        for mod in corr_config["SRlike_cuts"]:
            config["target_processes"][process]["SRlike_cuts"][mod] = corr_config[
                "SRlike_cuts"
            ][mod]
    if "ARlike_cuts" in corr_config:
        for mod in corr_config["ARlike_cuts"]:
            config["target_processes"][process]["ARlike_cuts"][mod] = corr_config[
                "ARlike_cuts"
            ][mod]
    if "AR_SR_cuts" in corr_config and to_AR_SR:
        for mod in corr_config["AR_SR_cuts"]:
            config["target_processes"][process]["ARlike_cuts"][mod] = corr_config[
                "AR_SR_cuts"
            ][mod]
            config["target_processes"][process]["SRlike_cuts"][mod] = corr_config[
                "AR_SR_cuts"
            ][mod]
