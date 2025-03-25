"""
Collection of helpful functions for the fake factor calculation scripts
"""

import array
import itertools as itt
import logging
import random
from contextlib import contextmanager
from typing import Any, Dict, Generator, List, Tuple, Union

import numpy as np
import ROOT

import helper.fitting_helper as fitting_helper
import helper.weights as weights
from configs.general_definitions import random_seed
from helper.hooks_and_patches import _EXTRA_PARAM_FLAG, _EXTRA_PARAM_MEANS, _EXTRA_PARAM_COUNTS


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


def get_split_combinations(
    categories: Dict[str, List[str]],
    binning: Union[
        List[float],
        Dict[str, List[float]],
        Dict[str, Dict[str, List[float]]],
        None,
    ],
    convert_binning_to_dict: bool = False,
    bandwidth: Union[
        float,
        Dict[str, float],
        None,
    ] = None,
) -> Tuple[List[str], List[Dict[str, str]], List[List[float]]]:
    """
    This function generates a dictionary for all category cut combinations.
    Categories can be defined based on one or two variables (more are not supported).
    Each variable has a list of cuts it should be split into.

    Further, if binning for individual categories is provided, the function will return the
    binning for each category cut combination.
    If binning is a list, the same binning is used for all category cut combinations.
    If binning is a dictionary, the binning is defined for each variable separately, where
    the first variable is the key of the first dictionary and the second variable is the key of the second dictionary.

    In case of two variables, the binning needs to be defined at least for the first variable,
    the second splitting variable can be defined in a nested dictionary, if not then the default binning of first
    variable is used for the second variable.

    Analogous to the binning, the bandwidth for the kernel regression can be defined for each category cut combination.
    (Supported up to one dimensional category splitting)

    Args:
        categories: Dictionary with the category definitions
        binning: Binning for the dependent variable
        convert_binning_to_dict: Boolean to convert the binning to a dictionary of binnings at the end
        bandwidth: Bandwidth for the kernel regression

    Return:
        1. List of variables the categories are defined in,
        2. List of all combinations of variable splits
        3. List of binning for each category cut combination
    """
    results, combinations, binnings = [], [], []
    split_variables = list(categories.keys())

    assert len(split_variables) <= 2, "Category splitting is only defined up to 2 dimensions."
    results.append(split_variables)

    combinations = [
        dict(zip(split_variables, v))
        for v in itt.product(*(categories[_v] for _v in split_variables))
    ]

    assert len(combinations) > 0, "No category combinations defined"
    results.append(combinations)

    if isinstance(binning, list):
        binnings = [binning] * len(combinations)
    elif isinstance(binning, dict):
        for values in (c.values() for c in combinations):
            values = list(values)
            _binning = binning.get(values[0])
            if len(values) == 1:
                binnings.append(_binning)
            elif len(values) == 2 and isinstance(_binning, dict):
                binnings.append(_binning.get(values[1]))
            elif len(values) == 2 and isinstance(_binning, list):
                logging.warning(f"Using default binning for {values} of {_binning}")
                binnings.append(_binning)
            else:
                raise Exception("Invalid type for binning")
    else:
        raise Exception("Invalid type for binning")

    assert len(combinations) == len(binnings), "Length of combinations and binnings do not match"

    if convert_binning_to_dict:
        _binnings = {}
        for split, _binning in zip(combinations, binnings):
            keys = [f"{var}#{split[var]}" for var in split_variables]
            if len(keys) == 1:
                _binnings.update({keys[0]: _binning})
            elif len(keys) == 2:
                _binnings.setdefault(keys[0], {})[keys[1]] = _binning
        binnings = _binnings
    results.append(binnings)

    if bandwidth is not None:
        if isinstance(bandwidth, dict):
            results.append([bandwidth.get(list(v)[0]) for v in (c.values() for c in combinations)])
        else:
            results.append([bandwidth] * len(combinations))

    return tuple(results)


def apply_region_filters(
    rdf: Any,
    channel: str,
    sample: str,
    category_cuts: Union[Dict[str, str], None],
    region_cuts: Dict[str, str],
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
                tmp[cut] = f"({cut} {a}) {op} ({cut} {b})"
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
        if sample not in ["data", "embedding"]:
            rdf = weights.apply_pNet_weight(rdf=rdf)
        if sample not in ["data", "embedding"] and "nbtag" not in sum_cuts.keys():
            rdf = weights.apply_btag_weight(rdf=rdf)
        rdf = rdf.Filter(f"({sum_cuts['bb_selection']})", "cut on bb pair")
        if sample not in ["data", "embedding"]:
            rdf = weights.apply_pNet_weight(rdf=rdf)

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
            (x, y, y_err_up, y_err_down, x_err_up, x_err_down) 
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

    y_fit, y_fit_up, y_fit_down = [], [], []
    if do_mc_subtr_unc:
        y_fit_mc_up, y_fit_mc_down = [], []

    x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000 * ff_hist.GetNbinsX())
    for value in x_fit:
        nominal, up, down = (
            callable_expression["nominal"](value),
            callable_expression["unc_up"](value),
            callable_expression["unc_down"](value),
        )
        y_fit.append(nominal)
        y_fit_up.append(up - nominal)
        y_fit_down.append(nominal - down)
        if do_mc_subtr_unc:
            mc_substr_up, mc_substr_down = (
                callable_expression["mc_subtraction_unc_up"](value),
                callable_expression["mc_subtraction_unc_down"](value),
            )
            y_fit_mc_up.append(abs(mc_substr_up - nominal))
            y_fit_mc_down.append(abs(nominal - mc_substr_down))

    x_fit = array.array("d", x_fit)
    y_fit = array.array("d", y_fit)
    y_fit_up = array.array("d", y_fit_up)
    y_fit_down = array.array("d", y_fit_down)
    if do_mc_subtr_unc:
        y_fit_mc_up = array.array("d", y_fit_mc_up)
        y_fit_mc_down = array.array("d", y_fit_mc_down)

    results, _args = {}, (len(x_fit), x_fit, y_fit, 0, 0)
    results["fit_graph_unc"] = ROOT.TGraphAsymmErrors(*_args, y_fit_down, y_fit_up)
    if do_mc_subtr_unc:
        results["fit_graph_mc_sub"] = ROOT.TGraphAsymmErrors(
            *_args, y_fit_mc_down, y_fit_mc_up
        )

    return (
        build_TGraph(ff_hist, add_xerrors_in_graph=True) if convert else ff_hist,
        results,
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


def smooth_function(
    hist: Any,
    bin_edges: List[float],
    write_corrections: str,
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

    # sampling values for y based on a normal distribution with the bin yield as mean value and with the measured statistical uncertainty
    with rng_seed(seed=random_seed):
        n_samples = 20
        sampled_y = list()
        for idx in range(len(x)):
            sampled_y.append(np.random.normal(y[idx], error_y_up[idx], n_samples))
        sampled_y = np.array(sampled_y)
        sampled_y[sampled_y < 0.0] = 0.0

    # calculate widths
    fit_y_binned = list()

    n_bins = 100 * hist_bins
    for i in range(n_bins):
        fit_y_binned.append(list())

    eval_bin_edges, bin_step = np.linspace(
        bin_edges[0], bin_edges[-1], (n_bins + 1), retstep=True,
    )
    bin_half = bin_step / 2.0
    smooth_x = (eval_bin_edges + bin_half)[:-1]
    smooth_x = array.array("d", smooth_x)

    for sample in range(n_samples):
        y_arr = array.array("d", sampled_y[:, sample])
        graph = ROOT.TGraphAsymmErrors(len(x), x, y_arr, 0, 0, error_y_down, error_y_up)
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

    if write_corrections == "binwise":
        _bins = np.array(bin_edges)
        _nom = np.array(y)
        _up = _nom + np.array(error_y_up)
        _down = _nom - np.array(error_y_down)
    else:
        _bins = np.array(eval_bin_edges)
        _nom = np.array(smooth_y)
        _up = _nom + np.array(smooth_y_up)
        _down = _nom - np.array(smooth_y_down)

    _nom[_nom < 0] = 0
    _up[_up < 0] = 0
    _down[_down < 0] = 0

    corr_dict = {"edges": _bins, "nominal": _nom, "up": _up, "down": _down}

    smooth_graph = ROOT.TGraphAsymmErrors(
        len(smooth_x), smooth_x, smooth_y, 0, 0, smooth_y_down, smooth_y_up
    )
    return nominal_graph, smooth_graph, corr_dict
