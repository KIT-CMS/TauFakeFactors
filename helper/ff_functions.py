"""
Collection of helpful functions for the fake factor calculation scripts
"""

import array
import functools
import inspect
import itertools as itt
import logging
import os
import random
from contextlib import contextmanager
from dataclasses import dataclass
from io import StringIO
from typing import Any, Callable, Dict, Generator, Iterator, List, Tuple, Union

import numpy as np
import ROOT
from wurlitzer import STDOUT, pipes

import configs.general_definitions as gd
import helper.fitting_helper as fitting_helper
import helper.functions as func
import helper.weights as weights
from configs.general_definitions import random_seed
from helper.hooks_and_patches import (_EXTRA_PARAM_COUNTS, _EXTRA_PARAM_FLAG,
                                      _EXTRA_PARAM_MEANS)


def cache_rdf_snapshot(cache_dir: str = "./.RDF_CACHE") -> Callable:
    """
    A decorator to cache the result of a ROOT RDataFrame operation to a file. If the cached file
    already exists, it loads the RDataFrame from the file instead of recomputing it.

    The wrapped function must accept an RDataFrame as its first positional argument or as a keyword
    argument named 'rdf'. The decorator computes a hash based on the function's arguments (excluding
    the RDataFrame) to uniquely identify the cache file. The tree name used in the ROOT file is
    hardcoded as "ntuple". The function assumes that the wrapped function returns a ROOT RDataFrame.
    The cache file is a ROOT file with a name derived from a hash of the function's arguments.

    The decorator performs the following steps:
        1. Checks if the cache file for the given arguments exists in the cache directory.
        2. If the cache file exists, loads the RDataFrame from the file.
        3. If the cache file does not exist, calls the wrapped function to process the RDataFrame,
           saves the result to the cache file, and then reloads the RDataFrame from the file.

    Args:
        cache_dir (str): The directory where cached ROOT files will be stored. Defaults to "./.RDF_CACHE".

    Returns:
        Callable: A decorator that wraps a function which processes an RDataFrame.
    """
    def decorator(function: Callable) -> Callable:
        @functools.wraps(function)
        def wrapper(*args: Any, **kwargs: Any) -> ROOT.RDataFrame:
            log = logging.getLogger(kwargs.get("logger") or function.__module__ + '.' + function.__name__)
            tree_name = "ntuple"

            if 'rdf' in kwargs:
                base_rdf = kwargs['rdf']
                key_args = {k: v for k, v in kwargs.items() if k != 'rdf'}
            elif args:
                base_rdf = args[0]
                key_args = {'__args__': args[1:], **kwargs}
            else:
                raise ValueError("Could not find RDataFrame argument ('rdf').")

            key_args['_var_'] = func.CachingKeyHelper.get_assignment_variable()
            caller_frame = inspect.stack()[1]
            caller_info = f"{caller_frame.function}@{caller_frame.filename}:{caller_frame.lineno}"
            key_args["_caller_"] = caller_info

            if "logger" in key_args:
                del key_args["logger"]

            os.makedirs(cache_dir, exist_ok=True)

            cache_hash = func.CachingKeyHelper.generate_key((), func.CachingKeyHelper.make_hashable(key_args))
            cache_filepath = os.path.join(cache_dir, f"{cache_hash}.root")

            if os.path.exists(cache_filepath) and func.RuntimeVariables.USE_CACHED_INTERMEDIATE_STEPS:
                log.info(f"Using existent filtered Rdf: {cache_filepath}")
                return ROOT.RDataFrame(tree_name, cache_filepath)

            log.info(f"Creating filtered Rdf under: {cache_filepath}")

            cols = [str(c) for c in base_rdf.GetColumnNames()]
            filtered_rdf = function(*args, **kwargs)

            if filtered_rdf.Count().GetValue() == 0:
                log.warning("Filter resulted in zero events. Creating an empty snapshot with the correct schema.")
                f = ROOT.TFile(cache_filepath, "RECREATE")
                tree = ROOT.TTree(tree_name, tree_name)
                for c in cols:
                    arr = ROOT.std.vector('float')()
                    tree.Branch(c, arr)
                tree.Write()
                f.Close()
            else:
                snapshot_result = filtered_rdf.Snapshot(tree_name, cache_filepath, cols)
                snapshot_result.GetValue()  # force execution

            if not os.path.exists(cache_filepath):
                log.error(f"Snapshot failed, file {cache_filepath} not created")
                raise RuntimeError(f"Snapshot failed, file {cache_filepath} not created")

            return ROOT.RDataFrame(tree_name, cache_filepath)

        return wrapper
    return decorator


@dataclass
class SplitQuantitiesContainer:
    """
    Container for split quantities.

    Stores quantities associated with a given split combination, including variables being
    split, category definitions, current split, and the binning. It also caches derived
    options such as the correction_option, fit_option, bandwidth for smoothing, and
    additional keyword arguments used in fitting (limit_kwargs).

    Attributes:
        variables (list): List of variable names used for splitting.
        categories (dict): Dictionary defining available category cuts.
        split (dict): Dictionary holding  current split combination.
        var_bins (list): Edges of the bins for histograms creation
        _correction_option (Union[str, None]): Correction option (set explicitly or from defaults in general definitions).
        _fit_option (Union[list, None]): Fitting option(s) (set explicitly or from defaults).
        _bandwidth (Union[float, int, None]): Bandwidth for smoothing corrections, derived from var_bins.
        _limit_kwargs (Union[dict, None]): Additional configuration (limits) for fit functions, derived from var_bins and an optional histogram.
    """

    variables: list
    categories: dict
    split: dict
    var_bins: list
    # --------------------------
    _correction_option: Union[str, None] = None  # set or from general_definitions
    _fit_option: Union[list, None] = None  # set or from general_definitions
    _bandwidth: Union[float, int, None] = None  # derived from var_bins
    _limit_kwargs: Union[dict, None] = None  # derived from var_bins and optional hist

    def limit_kwargs(self, hist: Union[None, ROOT.TH1] = None) -> dict:
        """
        Returns the keyword arguments for fit function limits.

        If already set this value is returned from the instance. Otherwise, a default value is
        obtained from the binning configuration and optional histogram.

        Args:
            hist (Union[None, ROOT.TH1], optional): ROOT histogram from which to derive limits.

        Returns:
            dict: Dictionary of keyword arguments for the fit limits.
        """

        if self._limit_kwargs is not None:
            return self._limit_kwargs

        self._limit_kwargs = gd.get_default_fit_function_limit_kwargs(
            binning=self.var_bins,
            hist=hist,
        )

        return self._limit_kwargs

    @property
    def bandwidth(self) -> float:
        """
        Retrieves the bandwidth parameter used for smoothing fits.

        If already set, the value is returned from the instance. Otherwise, a default value is
        obtained from the binning configuration.

        Returns:
            float: The computed bandwidth value.
        """

        if self._bandwidth is not None:
            return self._bandwidth

        self._bandwidth = gd.get_default_bandwidth(self.var_bins)

        return self._bandwidth

    @property
    def fit_option(self) -> Union[str, list]:
        """
        Retrieves the fitting option(s) for applying the fit.

        If already set, the value is returned from the instance. Otherwise, a default value is
        obtained from the binning configuration.

        Returns:
            Union[str, list]: The fit option(s) used for fitting.
        """

        if self._fit_option is not None:
            return self._fit_option

        self._fit_option = gd.default_fit_option

        return self._fit_option

    @property
    def correction_option(self) -> Union[str, None]:
        """
        Retrieves the correction option to be used with correctionlib.

        If already set, the value is returned from the instance. Otherwise, a default value is
        obtained from the general definitions.

        Returns:
            Union[str, None]: The correction option.
        """

        if self._correction_option is not None:
            return self._correction_option

        self._correction_option = gd.default_correction_option

        return self._correction_option


class SplitQuantities:
    """
    Handles splitting of quantities based on a sprovided split from configuration dictionary.

    This class processes a configuration (from YAML) that defines how the data are to be split
    into categories based on split category cuts.
    Supports iteration over all split combinations (via __iter__), returning a
    SplitQuantitiesContainer for each split.

    New options can be added by:
      - Extending the configuration dictionary with a new key (e.g. "new_option").
      - Implementing a corresponding property accessor here (similar to fit_option, bandwidth, etc.).
      - Adding new option in the __iter__ method.
      - Adding the new option to the SplitQuantitiesContainer class with the appropriate logic.

    Attributes:
        config (dict): Dictionary containing the overall configuration.
        categories (dict): Dictionary with category definitions (split_categories).
        split_variables (list): List of variable names dictating the available split dimensions.
    """
    def __init__(
        self,
        config: dict,
    ) -> None:
        """
        Initialization using a configuration dictionary.

        Args:
            config (dict): Configuration with at least "split_categories".
        """

        self.config = config
        self.categories = config.get("split_categories", None)
        self.split_variables = list(self.categories.keys()) if self.categories else []

    def to_dict(self, item: Union[list, str]) -> dict:
        """
        Converts a list or a property name into a dictionary of split options.

        Creates a dictionary with keys formatted as "variable#value". For a single split
        variable, the resulting dictionary maps "variable#value" to corresponding value in
        item; for two or more, nested dictionaries are created.

        Args:
            item (Union[list, str]): List of values corresponding to each split combination or the name of an attribute in the config.

        Returns:
            dict: Dictionary with keys based on split variable names and their values.
        """

        if isinstance(item, str):
            try:
                item = getattr(self, item)
            except AttributeError:
                raise KeyError(f"Item {item} not found in SplitQuantities")

        collection = {}
        for split, it in zip(self.split, item):
            keys = [f"{var}#{split[var]}" for var in self.split_variables]
            if len(keys) == 1:
                collection.update({keys[0]: it})
            elif len(keys) == 2:
                collection.setdefault(keys[0], {})[keys[1]] = it

        return collection

    def _fill_dict_based(self, key: str) -> list:
        """
        Fill and return a list from a nested dictionary configuration for a given key.

        For each split combination, checks if first value is present in the configuration
        under 'key', then returns direct object or the corresponding nested value. Asserts
        the number of split combinations matching the length of the filled list.

        Args:
            key (str): The configuration key for which to fill the list.

        Returns:
            list: List with each element corresponds to the configuration for each split combination.
        """

        collection = []
        for split_combination in self.split:
            if split_combination is None:
                collection.append(self.config[key])
                continue

            temp_config = self.config[key]
            for variable_name in self.split_variables:  # transverse the split variables
                category_value = split_combination[variable_name]
                if isinstance(temp_config, dict) and category_value in temp_config:
                    temp_config = temp_config[category_value]
                else:  # old style or missing category, might apply to all sub-categories
                    break

            collection.append(temp_config)

        assert len(self) == len(collection), f"Length of split combinations and {key} do not match"

        return collection

    def _fill(self, key: str, variable_condition: bool) -> Any:
        """
        Generic helper to fill a configuration value for a given key from the config.

        Based on whether the item is a list, a dict, or missing, returns either a cycle over
        a single element, the dictionary-based filled list, or raises an exception if the
        type is invalid.

        Args:
            key (str): Configuration key.
            variable_condition (bool): Flag indicating whether the config value should be treated as a simple list-based option.

        Returns:
            Any: Iterator (itertools.cycle) over configuration values or a list of filled values.
        """

        if key not in self.config:
            return itt.cycle([None])
        elif variable_condition:
            return itt.cycle([self.config[key]])
        elif isinstance(self.config[key], dict):
            return self._fill_dict_based(key)
        elif isinstance(self.config[key], list):
            if len(self.config[key]) != len(self):
                raise ValueError(f"Length of {key} list does not match number of split combinations")
            return self.config[key]
        else:
            raise Exception(f"Invalid type for {key}")

    @property
    def split(self) -> list:
        """
        Generate all combinations of split options based on the split_variables and categories.
        if split_variables are defined else returns [None].

        This method handles both flat lists for categories and nested dictionaries.

        Returns:
            list: List of dictionaries; each dictionary maps split variable names to category values.
        """
        if hasattr(self, "_split") and self._split is not None:
            return self._split

        if not self.split_variables:
            self._split = [None]
            return self._split

        variable1_name = self.split_variables[0]
        variable1_categories = self.categories[variable1_name]

        combinations = [{variable1_name: cat} for cat in variable1_categories]

        for i in range(1, len(self.split_variables)):  # iteratively build combinations
            next_variable_name = self.split_variables[i]
            next_variable_categories_config = self.categories[next_variable_name]

            new_combinations = []
            for combo in combinations:
                if isinstance(next_variable_categories_config, dict):
                    parent_key = combo[self.split_variables[i - 1]]
                    if parent_key not in next_variable_categories_config:
                        raise KeyError(f"Category '{parent_key}' not found in '{next_variable_name}' split definition.")
                    sub_categories = next_variable_categories_config[parent_key]
                else:  # old style: same list of sub-categories -> all parents
                    sub_categories = next_variable_categories_config

                for sub_category in sub_categories:  # new combination
                    new_combination = combo.copy()
                    new_combination[next_variable_name] = sub_category
                    new_combinations.append(new_combination)

            combinations = new_combinations

        self._split = combinations
        return self._split

    @property
    def var_bins(self) -> List[float]:
        """
        Retrieves binning configuration ("var_bins") from the configuration.

        Returns:
            List[float]: Binning configuration as defined in the configuration (could be a list or a nested configuration).
        """

        if hasattr(self, "_var_bins") and self._var_bins is not None:
            return self._var_bins

        key = "var_bins"
        assert key in self.config, f"{key} must always be defined"
        binnings = self._fill(
            key=key,
            variable_condition=isinstance(self.config[key], list),
        )

        self._var_bins = binnings
        return self._var_bins

    @property
    def fit_option(self) -> Union[str, list, None]:
        """
        Retrieves the fit option(s) from the configuration.

        Returns:
            Union[str, list]: The fit option(s) to be used.
        """

        if hasattr(self, "_fit_option") and self._fit_option is not None:
            return self._fit_option

        key = "fit_option"
        self._fit_option = self._fill(
            key=key,
            variable_condition=(
                key in self.config
                and isinstance(self.config[key], (str, list))
            ),
        )

        return self._fit_option

    @property
    def limit_kwargs(self) -> Union[dict, None]:
        """
        Retrieves the limit keyword arguments for fitting from the configuration.

        Returns:
            Union[dict, None]: The limit keyword arguments for fitting.
        """

        if hasattr(self, "_limit_kwargs") and self._limit_kwargs is not None:
            return self._limit_kwargs

        key = "limit_kwargs"
        self._limit_kwargs = self._fill(
            key=key,
            variable_condition=(
                key in self.config
                and isinstance(self.config[key], dict) 
                and self.config[key].get("limit_x") is not None
            ),
        )

        return self._limit_kwargs

    @property
    def bandwidth(self) -> Union[float, int, None]:
        """
        Retrieves the bandwidth for smoothing from the configuration.

        Returns:
            Union[float, int, None]: The bandwidth value or None (if not set).
        """

        if hasattr(self, "_bandwidth") and self._bandwidth is not None:
            return self._bandwidth

        key = "bandwidth"
        self._bandwidth = self._fill(
            key=key,
            variable_condition=(
                key in self.config 
                and isinstance(self.config[key], (float, int))
            ),
        )

        return self._bandwidth

    @property
    def correction_option(self) -> Union[str, None]:
        """
        Retrieves the correction option for correctionlib from the configuration.

        Returns:
            Union[str, None]: The correction option (if set) or None.
        """

        if hasattr(self, "_correction_option") and self._correction_option is not None:
            return self._correction_option

        key = "correction_option"
        self._correction_option = self._fill(
            key=key,
            variable_condition=(
                key in self.config 
                and isinstance(self.config[key], str)
            ),
        )

        return self._correction_option

    def __iter__(self) -> Iterator[Any]:
        """
        Iterates over all split combinations and yields a SplitQuantitiesContainer for each.

        For each combination of:
            - split (unique category combination),
            - var_bins,
            - fit_option,
            - limit_kwargs,
            - bandwidth, and
            - correction_option,
        a new SplitQuantitiesContainer is constructed and yielded.

        Yields:
            SplitQuantitiesContainer: One container for each split combination.
        """

        for (
            split,
            var_bins,
            fit_option,
            limit_kwargs,
            bandwidth,
            correction_option,
        ) in zip(
            self.split,
            self.var_bins,
            self.fit_option,
            self.limit_kwargs,
            self.bandwidth,
            self.correction_option,
        ):
            yield SplitQuantitiesContainer(
                variables=self.split_variables,
                categories=self.categories,
                split=split,
                var_bins=var_bins,
                _fit_option=fit_option,
                _limit_kwargs=limit_kwargs,
                _bandwidth=bandwidth,
                _correction_option=correction_option,
            )

    def __len__(self) -> int:
        """
        Returns the number of split combinations.

        Returns:
            int: Total number of unique splits based on the defined split variables.
        """

        return len(self.split)

    def __getitem__(
        self,
        index: Union[int, slice],
    ) -> Union[SplitQuantitiesContainer, List[SplitQuantitiesContainer]]:
        """
        Retrieves one or more SplitQuantitiesContainer objects by index.

        Supports both integer indexing and slicing.

        Args:
            index (Union[int, slice]): Index or slice for the requested container(s).

        Returns:
            Union[SplitQuantitiesContainer, List[SplitQuantitiesContainer]]: The split container or list of containers.
        """

        if isinstance(index, int):
            return list(self)[index]
        elif isinstance(index, slice):
            return [list(self)[i] for i in range(len(self))[index]]
        else:
            raise TypeError("Index must be an integer or a slice.")


@contextmanager
def rng_seed(seed: int) -> Generator[None, None, None]:
    """
    Context manager to set the random seed for numpy and python random.
    This is used to ensure reproducibility of the results (and help with debugging).

    Args:
        seed (int, optional): Seed to be used. Defaults to 19.

    Yields:
        None
    """
    np_rng_state, py_rng_state = np.random.get_state(), random.getstate()

    np.random.seed(seed)
    random.seed(seed)
    try:
        yield
    finally:
        np.random.set_state(np_rng_state)
        random.setstate(py_rng_state)


def controlplot_samples(
    use_embedding: bool,
    add_qcd: bool = True,
) -> List[str]:
    """
    Returns the list of samples that should be used for the control plots.

    Args:
        use_embedding: Boolean to use embedding or MC for genuine tau processes
        add_qcd: Add QCD samples to the collection of samples to be plotted.

    Returns:
        List of samples used for controlplots
    """
    samples = [
        "diboson_J",
        "diboson_L",
        "Wjets",
        "ttbar_J",
        "ttbar_L",
        "DYjets_J",
        "DYjets_L",
        "ST_J",
        "ST_L",
    ]
    if add_qcd:
        samples.append("QCD")

    if use_embedding:
        samples.append("embedding")
    else:
        samples.extend(["diboson_T", "ttbar_T", "DYjets_T", "ST_T"])

    return samples


def fill_corrlib_expression(
    item: Union[List[dict], dict],
    split_variables: List[str],
    split: Union[None, dict] = None,
):
    """
    This function fills the correctionlib expressions with the results from the fake factor
    calculation, process fraction or non-closure correction calculation (if a split is
    provided). If a Dictionary is provided, the function will fill the correctionlib
    expressions for each category cut combination. If a List is provided, the function will
    fill a Dictionary with the individual results.

    Args:
        item (Union[List[dict], dict]): Results from the fake factor calculation, process
                                        fraction or non-closure correction calculation
        split_variables (List[str]): List of variables the categories are defined in
        split (Union[None, dict], optional): Dictionary with the category cut combinations.
                                             Defaults to None.

    Returns:
        Dict: Dictionary with the filled correctionlib expressions
    """
    results = {}
    if split is not None and not isinstance(item, list) and isinstance(item, dict):  # Single result from multiprocessing
        keys = [f"{var}#{split[var]}" for var in split_variables]
        if len(keys) == 1:
            results[keys[0]] = item
        elif len(keys) == 2:
            results.setdefault(keys[0], {})[keys[1]] = item
        else:
            raise Exception("Something went wrong with the category splitting.")

    elif split is None and isinstance(item, list) and all(isinstance(it, dict) for it in item):  # Multiple results from multiprocessing
        for it in item:
            if len(split_variables) == 1:
                results.update(it)
            elif len(split_variables) == 2:
                key = list(it.keys())[0]
                results.setdefault(key, {}).update(it[key])
    else:
        raise ValueError("Item can only be a list (of dictionaries) or a dictionary")

    return results


@cache_rdf_snapshot(cache_dir="./.RDF_CACHE")
def apply_region_filters(
    rdf: Any,
    channel: str,
    sample: str,
    category_cuts: Union[Dict[str, str], None],
    region_cuts: Dict[str, str],
    logger: str,
) -> Any:
    """
    Function which applies filters to a root DataFrame for the fake factor calculation. This includes the region cuts and the category splitting.
    Additionally some weights are applied in this step because they depend on the full event selection.

    Args:
        rdf: root DataFrame object
        channel: Analysis channel of the tau analysis e.g. "et", "mt" or "tt"
        sample: Name of the sample/process of the "rdf", needed to prevent weight application to data
        category_cuts: Dictionary of cuts for the category splitting, for inclusive category this is None
        region_cuts: Dictionary of cuts for a fake factor calculation region
        logger: Logger name for logging purposes

    Return:
        Filtered root DataFrame with adjusted weights
    """
    tmp = dict()
    if category_cuts is not None:
        for cut in category_cuts:
            if "#" not in category_cuts[cut]:
                tmp[cut] = f"{cut} {category_cuts[cut]}"
            else:
                a, op, b = category_cuts[cut].split("#")
                tmp[cut] = f"(({cut} {a}) {op} ({cut} {b}))"
    sum_cuts = {**tmp, **region_cuts}

    for cut in sum_cuts:
        if cut not in ["nbtag", "bb_selection"]:
            if "had_tau_id_vs_jet" in cut:
                wps = get_wps(cut_string=sum_cuts[cut])
                try:
                    idx = cut.rsplit("_")[5]
                except:
                    idx = None
                if sample not in ["data"]:
                    rdf = weights.apply_tau_id_vsJet_weight(
                        rdf=rdf, channel=channel, wps=wps, idx=idx
                    )
                rdf = rdf.Filter(f"({sum_cuts[cut]})", f"cut on {cut}")

            elif (
                "had_boostedtau_id_iso" in cut
            ):  # this is only relevant for an analysis with boosted tau pairs
                wp = get_wps(cut_string=sum_cuts[cut])
                try:
                    idx = cut.rsplit("_")[4]
                except:
                    idx = None
                if sample not in ["data"]:
                    rdf = weights.apply_boostedtau_id_iso_weight(
                        rdf=rdf,
                        channel=channel,
                        cut_string=sum_cuts[cut],
                        wp=wp,
                        idx=idx,
                    )
                rdf = rdf.Filter(f"({sum_cuts[cut]})", f"cut on {cut}")

            else:
                rdf = rdf.Filter(f"({sum_cuts[cut]})", f"cut on {cut}")
    # cut on number of b-tagged jets needs to be the last cut to do an on-the-fly calculation of the b-tagger weight
    if "nbtag" in sum_cuts.keys():
        if sample not in ["data", "embedding"]:
            rdf = weights.apply_btag_weight(rdf=rdf)
        rdf = rdf.Filter(f"({sum_cuts['nbtag']})", "cut on nbtag")
    if "bb_selection" in sum_cuts.keys():
        if sample not in ["data", "embedding"] and "fj_Xbb" in sum_cuts["bb_selection"]:
            rdf = weights.apply_pNet_weight(rdf=rdf)
        if sample not in ["data", "embedding"] and "nbtag" not in sum_cuts.keys():
            rdf = weights.apply_btag_weight(rdf=rdf)
        rdf = rdf.Filter(f"({sum_cuts['bb_selection']})", "cut on bb pair")

    log = logging.getLogger(logger)
    # redirecting C++ stdout for Report() to python stdout
    out = StringIO()
    with pipes(stdout=out, stderr=STDOUT):
        rdf.Report().Print()
    log.info(out.getvalue())
    log.info("-" * 50)

    return rdf


def get_wps(cut_string: str) -> Union[List[str], str]:
    """
    This function reads out the tau id vs jet working points based on the cut string.
    The names of the tau id variable is expected to be "id_tau_vsJet_{WP}_*".

    Args:
        cut_string: String with the tau id vs jet cuts

    Return:
        One working point or if the cut is between two working points a list of two working points
    """
    if "&&" in cut_string:
        wp_1 = cut_string.rsplit("&&")[0].rsplit("_")[3]
        wp_2 = cut_string.rsplit("&&")[1].rsplit("_")[3]
        wps = [wp_1, wp_2]
    else:
        wps = cut_string.rsplit("_")[3]

    return wps


def QCD_SS_estimate(hists: Dict[str, Any]) -> Any:
    """
    Function which estimates QCD as data minus all other background processes (from MC) exept QCD.
    For Higgs to tautau analysis this is normally done within the same sign selection to estimate the opposite sign QCD contribution.

    Args:
        hists: Dictionary with histograms for all relevant processes

    Return:
        Histogram with the estimated QCD yield
    """
    qcd = hists["data"].Clone()
    for sample in hists:
        if sample not in ["data", "data_subtracted", "data_ff", "QCD"]:
            qcd.Add(hists[sample], -1)

    # check for negative bins
    for i in range(1, qcd.GetNbinsX() + 1):
        if qcd.GetBinContent(i) < 0.0:
            qcd.SetBinContent(i, 0.0)

    return qcd


def calculate_QCD_FF(
    SRlike: Dict[str, Any], ARlike: Dict[str, Any]
) -> Tuple[Any, Any, Any]:
    """
    Function which calculates fake factors based on the histograms from the determination regions.
    Additionally up and down variations are calculated for the MC subtraction.
    This function is specified for QCD.

    Args:
        SRlike: Dictionary with histograms from the signal-like determination region for all relevant processes
        ARlike: Dictionary with histograms from the application-like determination region for all relevant processes

    Return:
        Ratio histograms for nominal, up and down variations (MC subtraction)
    """
    ratio = SRlike["data_subtracted"].Clone()
    ratio.Divide(ARlike["data_subtracted"])
    ratio_up = SRlike["data_subtracted_up"].Clone()
    ratio_up.Divide(ARlike["data_subtracted_up"])
    ratio_down = SRlike["data_subtracted_down"].Clone()
    ratio_down.Divide(ARlike["data_subtracted_down"])

    return ratio, ratio_up, ratio_down


def calculate_Wjets_FF(
    SRlike: Dict[str, Any], ARlike: Dict[str, Any]
) -> Tuple[Any, Any, Any]:
    """
    Function which calculates fake factors based on the histograms from the determination regions.
    Additionally up and down variations are calculated for the MC subtraction.
    This function is specified for Wjets.

    Args:
        SRlike: Dictionary with histograms from the signal-like determination region for all relevant processes
        ARlike: Dictionary with histograms from the application-like determination region for all relevant processes

    Return:
        Ratio histograms for nominal, up and down variations (MC subtraction)
    """
    ratio = SRlike["data_subtracted"].Clone()
    ratio.Divide(ARlike["data_subtracted"])
    ratio_up = SRlike["data_subtracted_up"].Clone()
    ratio_up.Divide(ARlike["data_subtracted_up"])
    ratio_down = SRlike["data_subtracted_down"].Clone()
    ratio_down.Divide(ARlike["data_subtracted_down"])

    return ratio, ratio_up, ratio_down


def calculate_ttbar_FF(
    SR: Dict[str, Any],
    AR: Dict[str, Any],
    SRlike: Dict[str, Any],
    ARlike: Dict[str, Any],
) -> Any:
    """
    Function which calculates fake factors based on the histograms from the signal and application regions.
    No additional up and down variations are calculated for the MC subtraction because the ttbar fake factors
    are calculated based on MC only. The measured fake factors are scaled with a ratio between inclusive fake factors
    from data and MC in the determination regions.

    Args:
        SR: Dictionary with histograms from the signal region for all relevant processes
        AR: Dictionary with histograms from the application region for all relevant processes
        SRlike: Dictionary with histograms from the signal-like determination region for all relevant processes
        ARlike: Dictionary with histograms from the application-like determination region for all relevant processes

    Return:
        Ratio histogram for nominal case
    """
    ratio_mc = SR["ttbar_J"].Clone()
    ratio_mc.Divide(AR["ttbar_J"])

    ratio_DR_data = SRlike["data_subtracted"].Clone()
    ratio_DR_data.Divide(ARlike["data_subtracted"])

    ratio_DR_mc = SRlike["ttbar_J"].Clone()
    ratio_DR_mc.Divide(ARlike["ttbar_J"])

    sf = ratio_DR_data.GetMaximum() / ratio_DR_mc.GetMaximum()
    ratio_mc.Scale(sf)

    return ratio_mc


def calc_fraction(hists: Dict[str, Any], target: str, processes: List[str]) -> Any:
    """
    This function calculates a fraction histogram for a target process out of a list of all considered processes.
    Later the sum of all fraction histograms results to 1 (in each bin).

    Args:
        hists: Dictionary with histograms for all considered processes
        target: Target process for which the fraction should be calculated
        processes: List of all cosidered processes for the fraction calculation

    Return:
        Fraction histogram
    """
    mc = hists[target].Clone()
    for p in processes:
        if p != target:
            mc.Add(hists[p])

    frac = hists[target].Clone()
    frac.Divide(mc)

    return frac


def add_fraction_variations(
    hists: Dict[str, Any], processes: List[str]
) -> Dict[str, Dict[str, Any]]:
    """
    Function which calculates variations of fraction histograms. The fraction of each process
    is once varied up by 7% and all other processes are varied down by 7%. The same is done
    the other way around to get the down variations.

    Args:
        hists: Dictionary with fraction histograms for all considered processes
        processes: List of all cosidered processes for the fraction calculation

    Return:
        Dictionary with nominal, up and down variations of all considered processes
    """
    variations = dict()
    variations["nominal"] = dict(hists)

    for p in processes:
        hists_up = dict(hists)
        for hist in hists_up:
            hists_up[hist] = hists_up[hist].Clone()
            if p == hist:
                hists_up[hist].Scale(1.07)
            else:
                hists_up[hist].Scale(0.93)
        variations["frac_{}_up".format(p)] = hists_up

    for p in processes:
        hists_down = dict(hists)
        for hist in hists_down:
            hists_down[hist] = hists_down[hist].Clone()
            if p == hist:
                hists_down[hist].Scale(0.93)
            else:
                hists_down[hist].Scale(1.07)
        variations["frac_{}_down".format(p)] = hists_down

    return variations


def get_yields_from_hists(
    hists: Dict[str, Dict[str, Any]], processes: List[str]
) -> Dict[str, Dict[str, Dict[str, List[float]]]]:
    """
    This function transforms fraction histograms (with variations) obtained from calc_fraction() and add_fraction_variations()
    into lists of fraction values which are later used to produce in the production of the correctionlib file.

    Args:
        hists: Dictionary with all categories which include nominal, up and down fraction histograms of all considered processes
        processes: List of all cosidered processes for the fraction calculation
    Return:
        Dictionary with all categories which include nominal, up and down fraction values of all considered processes
    """
    fracs = dict()
    categories = hists.keys()

    for category in categories:
        v_fracs = dict()
        for variation in hists[category]:
            p_fracs = dict()
            for p in processes:
                h = hists[category][variation][p]
                l = list()
                for b in range(1, h.GetNbinsX() + 1):
                    l.append(h.GetBinContent(b))
                p_fracs[p] = l
            v_fracs[variation] = p_fracs
        fracs[category] = v_fracs

    return fracs


def build_TGraph(
    hist: Union[ROOT.TH1, None] = None,
    return_components: bool = False,
    add_xerrors_in_graph: bool = False,
) -> Union[ROOT.TGraphAsymmErrors, Tuple[ROOT.TGraphAsymmErrors, List[Any]], None]:
    """
    Function which builds a TGraph from a histogram. The TGraph is built from the bin
    centers and the bin contents of the histogram. If the histogram is None, None is
    returned, if the histogram has the _EXTRA_PARAM_FLAG set to True, the bin centers are
    taken from the _EXTRA_PARAM_MEANS attribute of the histogram.

    Args:
        hist: Histogram to be converted to TGraph
        return_components: Boolean to return the components of the TGraph as well
        add_xerrors_in_graph: Boolean to add x errors to the TGraph, if not set, the x errors are set to 0
    Return:
        1. TGraph of the histogram,
        2. List of components of the TGraph if return_components is True containing:
            (x, y, x_err_down, x_err_up, y_err_down, y_err_up) 
    """
    if hist is None:
        return None

    nbins = hist.GetNbinsX()
    x, y, y_err_up, y_err_down, x_err_up, x_err_down = [], [], [], [], [], []

    for nbin in range(nbins):
        if hasattr(hist, _EXTRA_PARAM_FLAG) and getattr(hist, _EXTRA_PARAM_FLAG):
            nominal = getattr(hist, _EXTRA_PARAM_MEANS)[nbin]
        else:
            nominal = hist.GetBinCenter(nbin + 1)
        x.append(nominal)
        x_err_down.append(nominal - hist.GetBinLowEdge(nbin + 1))
        x_err_up.append(hist.GetBinLowEdge(nbin + 2) - nominal)

        y.append(hist.GetBinContent(nbin + 1))
        y_err_down.append(hist.GetBinErrorLow(nbin + 1))
        y_err_up.append(hist.GetBinErrorUp(nbin + 1))

    x, y = array.array("d", x), array.array("d", y)
    y_err_up, y_err_down = array.array("d", y_err_up), array.array("d", y_err_down)
    x_err_up, x_err_down = array.array("d", x_err_up), array.array("d", x_err_down)

    if add_xerrors_in_graph:
        args = (x, y, x_err_down, x_err_up, y_err_down, y_err_up)
    else:
        args = (x, y, 0, 0, y_err_down, y_err_up)

    if return_components:
        return (ROOT.TGraphAsymmErrors(nbins, *args), *args)
    else:
        return ROOT.TGraphAsymmErrors(nbins, *args)


def fit_function(
    ff_hists: Union[List[Any], Any],
    bin_edges: List[int],
    logger: str,
    fit_option: Union[str, List[str]],
    limit_kwargs: Dict[str, Any],
) -> Tuple[
    Union[ROOT.TH1, ROOT.TGraphAsymmErrors],
    Dict[str, Any],
    Dict[str, str],
    str,
]:
    """
    This function performs fits of the ratio histogram. The fitted function is then used
    to produce expressions for correctionlib (including variations). Additionally graphs of
    the fitted function variations are generated which are later used for plotting purposes.

    Args:
        ff_hists: Either a list of nominal and MC varied ratio histograms or only the nominal ratio histogram
        bin_edges: Bins edges of the fitted variable, needed for the graphs for plotting
        logger: Name of the logger that should be used
        fit_option: List[str] correspond to a list of poly_n fits to be performed
                    with best fit being the one with the lowest chi2/ndf value and nominal and
                    variations being > 0 in fit range.
                    str: "poly_n" or "binwise", where n is the order of the polynomial fit

    Return:
        1. Dictionary with the fitted function for each variation,
        2. Dictionary with graphs for each variation,
        3. Dictionary of function expressions for correctionlib (nominal and variations)
    """
    if not isinstance(fit_option, list) and fit_option != "binwise":
        fit_option = [fit_option]

    do_mc_subtr_unc, ff_hist_up, ff_hist_down = False, None, None
    if isinstance(ff_hists, list):
        ff_hist, ff_hist_up, ff_hist_down = ff_hists
        do_mc_subtr_unc = True
    else:
        ff_hist = ff_hists

    retrival_function, convert = fitting_helper.get_wrapped_functions_from_fits, True
    if fit_option == "binwise":
        retrival_function, convert = fitting_helper.get_wrapped_hists, False

    callable_expression, correctionlib_expression, used_fit = retrival_function(
        bounds=(bin_edges[0], bin_edges[-1]),
        ff_hist=build_TGraph(ff_hist) if convert else ff_hist,
        ff_hist_up=build_TGraph(ff_hist_up) if convert else ff_hist_up,
        ff_hist_down=build_TGraph(ff_hist_down) if convert else ff_hist_down,
        do_mc_subtr_unc=do_mc_subtr_unc,
        logger=logger,
        function_collection=fit_option,
        verbose=True,
        limit_x=limit_kwargs["limit_x"],
    )

    return (
        build_TGraph(ff_hist, add_xerrors_in_graph=True) if convert else ff_hist,
        callable_expression,
        correctionlib_expression,
        used_fit,
    )


def calculate_non_closure_correction(
    SRlike: Dict[str, Any], ARlike: Dict[str, Any]
) -> Tuple[Any, Any]:
    """
    Function which calculates non closure corrections based on the histograms from the determination regions.

    Args:
        SRlike: Dictionary with histograms from the signal-like determination region for all relevant processes
        ARlike: Dictionary with histograms from the application-like determination region for all relevant processes

    Return:
        1. Ratio histogram of data (MC subtracted) in a signal-like region and data (scaled to MC subtracted) with applied fake factors in an application-like region,
        2. Process fraction in the application-like region
    """
    frac = ARlike["data_subtracted"].GetMaximum() / ARlike["data"].GetMaximum()
    predicted = ARlike["data_ff"].Clone()
    predicted.Scale(frac)

    corr = SRlike["data_subtracted"].Clone()
    corr.Divide(predicted)
    hist_bins = corr.GetNbinsX()

    # check for empty bins, if empty the bin value is set to 1 +/- 1
    for x in range(1, hist_bins + 1):
        bincontent = corr.GetBinContent(x)
        if bincontent <= 0.0:
            if hasattr(corr, _EXTRA_PARAM_FLAG) and getattr(corr, _EXTRA_PARAM_FLAG):
                getattr(
                    corr,
                    _EXTRA_PARAM_MEANS,
                )[x - 1] = (corr.GetBinLowEdge(x) + corr.GetBinLowEdge(x + 1)) / 2
                getattr(corr, _EXTRA_PARAM_COUNTS)[x - 1] = 1.0
            corr.SetBinContent(x, 1.0)
            corr.SetBinError(x, 1.0)

    return corr, frac


def calculate_non_closure_correction_Wjets_fromMC(
    SRlike: Dict[str, Any], ARlike: Dict[str, Any]
) -> Any:
    """
    Function which calculates non closure corrections with MC based on the histograms from the determination regions for the ttbar process.

    Args:
        SRlike: Dictionary with histograms from the signal-like determination region for all relevant processes
        ARlike: Dictionary with histograms from the application-like determination region for all relevant processes

    Return:
        Ratio histogram of MC in a signal-like region MC and with applied fake factors in an application-like region
    """
    corr = SRlike["Wjets"].Clone()
    predicted = ARlike["Wjets_ff"].Clone()
    corr.Divide(predicted)

    return corr


def calculate_non_closure_correction_ttbar_fromMC(
    SR: Dict[str, Any], AR: Dict[str, Any]
) -> Any:
    """
    Function which calculates non closure corrections with MC based on the histograms from the signal/application regions for the ttbar process.

    Args:
        SR: Dictionary with histograms from the signal region for all relevant processes
        AR: Dictionary with histograms from the application region for all relevant processes

    Return:
        Ratio histogram of MC in a signal region and MC with applied fake factors in an application region
    """
    corr = SR["ttbar_J"].Clone()
    predicted = AR["ttbar_ff"].Clone()
    corr.Divide(predicted)

    return corr


def _get_index_and_slices(
    correction_option: str,
) -> Tuple[int, int, slice, slice, slice, bool]:
    """
    This function derives the bin index and slices for the hybrid correction,
    in case of no hybrid correction, default values are returned.

    Args:
        correction_option: Correction option string

    Return:
        1. Start index of the left part of the hybrid correction
        2. End index of the right part of the hybrid correction
        3. Slice for the left part of the hybrid correction
        4. Slice for the right part of the hybrid correction
        5. Slice for the overlap part of the hybrid correction
        6. Boolean if the hybrid correction
    """

    # default values
    start_idx, end_idx, has_left_part, has_right_part = 0, -1, False, False
    left_slice, right_slice, overlap_slice = slice(0, 0), slice(-1, -1), slice(None, None)

    if "binwise" in correction_option and "smoothed" in correction_option:
        binwise_part = next(it for it in correction_option.split("+") if "binwise" in it)
        bin_idx_groups = [eval(expr) for expr in binwise_part.split("#")[1:]]

        if len(bin_idx_groups) == 2:
            start_idx, end_idx = bin_idx_groups[0][-1] + 1, bin_idx_groups[1][-1] - 1
            has_left_part, has_right_part = True, True
        elif len(bin_idx_groups) == 1:
            group = bin_idx_groups[0]
            if all(it >= 0 for it in group):
                start_idx, has_left_part = group[-1] + 1, True
            elif all(it < 0 for it in group):
                end_idx, has_right_part = group[-1] - 1, True
            else:
                raise ValueError("Invalid bin index for the binned part of the correction")

        left_slice = slice(None, start_idx) if has_left_part else slice(0, 0)
        right_slice = slice(end_idx + 1, None) if has_right_part else slice(-1, -1)
        overlap_slice = slice(None, -1) if has_right_part else slice(None, None)

    return (
        start_idx,
        end_idx,
        left_slice,
        right_slice,
        overlap_slice,
        has_left_part or has_right_part,
    )


def smooth_function(
    hist: Any,
    bin_edges: List[float],
    correction_option: str,
    bandwidth: float,
) -> Tuple[Any, Dict[str, np.ndarray]]:
    """
    This function performs a smoothing fit of a histogram. Smoothing is mainly used for the corrections of the fake factors.
    For smoothing the Nadaraya-Watson kernel regression estimator is used.

    Args:
        hist: Histogram of a correction
        bin_edges: Bin edges of "hist"

    Return:
        1. root TGraph of the smoothed function,
        2. Dictionary of arrays with information about the smoothed function values to be stored with correctionlib (nominal and variations)
    """
    hist_bins = hist.GetNbinsX()

    # transforming bin information to arrays
    nominal_graph, x, y, _, _, error_y_down, error_y_up = build_TGraph(
        hist, return_components=True, add_xerrors_in_graph=True,
    )

    if "skip" in correction_option:
        return (
            nominal_graph,
            nominal_graph,
            {
                "edges": np.array(bin_edges),
                "nominal": np.array([1.0]),
                "up": np.array([1.0]) + 1e-6,
                "down": np.array([1.0]) - 1e-6,
            },
        )

    # values for slicing and range determination in case of hybrid correction
    (
        start_idx,
        end_idx,
        left_slice,
        right_slice,
        overlap_slice,
        is_hybrid_correction,
    ) = _get_index_and_slices(correction_option)

    # sampling values for y based on a normal distribution with the bin yield as mean
    # value and with the measured statistical uncertainty
    n_samples = 20

    with rng_seed(seed=random_seed):
        sampled_y = np.array([np.random.normal(*it, n_samples) for it in zip(y, error_y_up)])
        sampled_y[sampled_y < 0.0] = 0.0

    n_bins = 100 * hist_bins
    smooth_x = array.array("d", np.linspace(x[start_idx], x[end_idx], n_bins + 1))
    fit_y_binned = [[] for _ in range(n_bins)]

    for sample in sampled_y.T:
        graph = ROOT.TGraphAsymmErrors(len(x), x, array.array("d", sample), 0, 0, error_y_down, error_y_up)
        gs = ROOT.TGraphSmooth("normal")
        grout = gs.SmoothKern(graph, "normal", bandwidth, n_bins, smooth_x)
        for i in range(n_bins):
            fit_y_binned[i].append(grout.GetPointY(i))

    smooth_y, smooth_y_up, smooth_y_down = [], [], []

    for b in fit_y_binned:
        smooth_y.append(np.mean(b))
        smooth_y_up.append(np.std(b))
        smooth_y_down.append(np.std(b))

    smooth_y = array.array("d", smooth_y)
    smooth_y_up = array.array("d", smooth_y_up)
    smooth_y_down = array.array("d", smooth_y_down)

    if correction_option == "binwise":
        _bins = np.array(bin_edges)
        _nom = np.array(y)
        _up = _nom + np.array(error_y_up)
        _down = _nom - np.array(error_y_down)
    elif correction_option == "smoothed":
        _bins = np.array(smooth_x)
        _nom = np.array(smooth_y)
        _up = _nom + np.array(smooth_y_up)
        _down = _nom - np.array(smooth_y_down)
    elif is_hybrid_correction:
        _bins = np.concatenate(
            (
                bin_edges[left_slice],
                [bin_edges[start_idx]] if start_idx != 0 else [],
                np.array(smooth_x)[overlap_slice],
                [bin_edges[end_idx]] if end_idx != -1 else [],
                bin_edges[right_slice],
            ),
        )
        _nom = np.concatenate(
            (
                y[left_slice],
                [smooth_y[0]] if start_idx != 0 else [],
                np.array(smooth_y)[overlap_slice],
                [smooth_y[-1]] if end_idx != -1 else [],
                y[right_slice],
            ),
        )
        _up = _nom + np.concatenate(
            (
                error_y_up[left_slice],
                [smooth_y_up[0]] if start_idx != 0 else [],
                np.array(smooth_y_up)[overlap_slice],
                [smooth_y_up[-1]] if end_idx != -1 else [],
                error_y_up[right_slice],
            ),
        )
        _down = _nom - np.concatenate(
            (
                error_y_down[left_slice],
                [smooth_y_down[0]] if start_idx != 0 else [],
                np.array(smooth_y_down)[overlap_slice],
                [smooth_y_down[-1]] if end_idx != -1 else [],
                error_y_down[right_slice],
            ),
        )
    else:
        raise ValueError("Invalid correction option")

    _nom[_nom < 0] = 0
    _up[_up < 0] = 0
    _down[_down < 0] = 0

    # adjustment to bin edges
    if _bins[0] != bin_edges[0]:
        _prepend = lambda a, b = None: np.concatenate(([a[0] if b is None else b[0]], a))
        _bins = _prepend(_bins, bin_edges)
        _nom, _up, _down = _prepend(_nom), _prepend(_up), _prepend(_down)
    if _bins[-1] != bin_edges[-1]:
        _append = lambda a, b = None: np.concatenate((a, [a[-1] if b is None else b[-1]]))
        _bins = _append(_bins, bin_edges)
        _nom, _up, _down = _append(_nom), _append(_up), _append(_down)

    corr_dict = {"edges": _bins, "nominal": _nom, "up": _up, "down": _down}

    smooth_graph = ROOT.TGraphAsymmErrors(
        len(smooth_x), smooth_x, smooth_y, 0, 0, smooth_y_down, smooth_y_up
    )
    return nominal_graph, smooth_graph, corr_dict
