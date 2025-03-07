"""
Collection of helpful functions for the fake factor calculation scripts
"""

import logging
import array
from typing import Any, Dict, List, Tuple, Union

import numpy as np
import ROOT

import itertools as itt
from copy import deepcopy
import helper.fitting_helper as fitting_helper
import helper.weights as weights


class Defaults:
    @staticmethod
    def fit_function_limit_kwargs(binning):
        return {
            "limit_x": {
                "nominal": (binning[0], binning[-1]),
                "up": (-float("inf"), float("inf")),
                "down": (-float("inf"), float("inf")),
            },
        }


def get_split_combinations(
    categories: Dict[str, List[str]],
    binning: Union[
        List[float],
        Dict[str, List[float]],
        Dict[str, Dict[str, List[float]]],
        None,
    ],
    convert_binning_to_dict: bool = False,
) -> Tuple[List[str], List[Dict[str, str]], List[List[float]]]:
    """
    This function generates a dictionary for all category cut combinations.
    Categories can be defined based on one or two variables (more are not supported).
    Each variable has a list of cuts it should be split into.

    Further, if binning is provided, the function will return the binning for each category cut combination.
    If binning is a list, the same binning is used for all category cut combinations.
    If binning is a dictionary, the binning is defined for each variable separately, where
    the first variable is the key of the first dictionary and the second variable is the key of the second dictionary.

    In case of two variables, the binning needs to be defined at least for the first variable,
    the second splitting variable can be defined in a nested dictionary, if not then the default binning of first
    variable is used for the second variable.

    Args:
        categories: Dictionary with the category definitions
        binning: Binning for the dependent variable

    Return:
        1. List of variables the categories are defined in,
        2. List of all combinations of variable splits
        3. List of binning for each category cut combination
    """
    combinations, binnings = [], []
    split_variables = list(categories.keys())

    assert len(split_variables) <= 2, "Category splitting is only defined up to 2 dimensions."

    combinations = [
        dict(zip(split_variables, v))
        for v in itt.product(*(categories[_v] for _v in split_variables))
    ]

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
                print(f"Using default binning for {values} of {_binning}")
                logging.warning(f"Using default binning for {values} of {_binning}")
                binnings.append(_binning)
            else:
                raise Exception("Invalid type for binning")
    else:
        raise Exception("Invalid type for binning")

    assert len(combinations) > 0, "No category combinations defined"
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

    return split_variables, combinations, binnings


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
            tmp[cut] = f"{cut} {category_cuts[cut]}"
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
    for i in range(qcd.GetNbinsX()):
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
                for b in range(h.GetNbinsX()):
                    l.append(h.GetBinContent(b + 1))
                p_fracs[p] = l
            v_fracs[variation] = p_fracs
        fracs[category] = v_fracs

    return fracs


def fit_function(
    ff_hists: Union[List[Any], Any],
    bin_edges: List[int],
    logger: str,
    fit_option: Union[str, List[str]],
    limit_kwargs: Dict[str, Any],
) -> Tuple[Dict[str, Any], Dict[str, str], str]:
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
        1. Dictionary with graphs for each variation,
        2. Dictionary of function expressions for correctionlib (nominal and variations)
    """
    if not isinstance(fit_option, list) and fit_option != "binwise":
        fit_option = [fit_option]

    do_mc_subtr_unc, ff_hist_up, ff_hist_down = False, None, None
    if isinstance(ff_hists, list):
        ff_hist, ff_hist_up, ff_hist_down = ff_hists
        do_mc_subtr_unc = True
    else:
        ff_hist = ff_hists

    # Producing a graph based on the provided root histograms
    nbins = ff_hist.GetNbinsX()
    x, y, error_y_up, error_y_down = [], [], [], []

    for nbin in range(nbins):
        x.append(ff_hist.GetBinCenter(nbin + 1))
        y.append(ff_hist.GetBinContent(nbin + 1))
        error_y_up.append(ff_hist.GetBinErrorUp(nbin + 1))
        error_y_down.append(ff_hist.GetBinErrorLow(nbin + 1))

    x = array.array("d", x)
    y = array.array("d", y)
    error_y_up = array.array("d", error_y_up)
    error_y_down = array.array("d", error_y_down)

    graph = ROOT.TGraphAsymmErrors(nbins, x, y, 0, 0, error_y_down, error_y_up)

    retrival_function = fitting_helper.get_wrapped_functions_from_fits
    if fit_option == "binwise":
        retrival_function = fitting_helper.get_wrapped_hists

    callable_expression, correctionlib_expression, used_fit = retrival_function(
        graph=graph,
        bounds=(bin_edges[0], bin_edges[-1]),
        ff_hist=ff_hist,
        ff_hist_up=ff_hist_up,
        ff_hist_down=ff_hist_down,
        do_mc_subtr_unc=do_mc_subtr_unc,
        logger=logger,
        function_collection=fit_option,
        verbose=True,
        limit_x=limit_kwargs["limit_x"],
    )

    y_fit, y_fit_up, y_fit_down = [], [], []
    if do_mc_subtr_unc:
        y_fit_mc_up, y_fit_mc_down = [], []

    x_fit = np.linspace(bin_edges[0], bin_edges[-1], 1000 * nbins)
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

    return results, correctionlib_expression, used_fit


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

    # check for empty bins, if empty the bin value is set to 1 +/- 0.6
    for x in range(hist_bins):
        bincontent = corr.GetBinContent(x)
        if bincontent <= 0.0:
            corr.SetBinContent(x, 1.0)
            corr.SetBinError(x, 0.6)

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
    hist: Any, bin_edges: List[float], write_corrections: str
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
    x = list()
    y = list()
    error_y_up = list()
    error_y_down = list()
    for nbin in range(hist_bins):
        x.append(hist.GetBinCenter(nbin + 1))
        y.append(hist.GetBinContent(nbin + 1))
        error_y_up.append(hist.GetBinErrorUp(nbin + 1))
        error_y_down.append(hist.GetBinErrorLow(nbin + 1))

    x = array.array("d", x)
    y = array.array("d", y)
    error_y_up = array.array("d", error_y_up)
    error_y_down = array.array("d", error_y_down)

    # sampling values for y based on a normal distribution with the bin yield as mean value and with the measured statistical uncertainty
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
        bin_edges[0], bin_edges[-1], (n_bins + 1), retstep=True
    )
    bin_range = bin_edges[-1] - bin_edges[0]
    bin_half = bin_step / 2.0
    smooth_x = (eval_bin_edges + bin_half)[:-1]
    smooth_x = array.array("d", smooth_x)

    for sample in range(n_samples):
        y_arr = array.array("d", sampled_y[:, sample])
        graph = ROOT.TGraphAsymmErrors(len(x), x, y_arr, 0, 0, error_y_down, error_y_up)
        gs = ROOT.TGraphSmooth("normal")
        grout = gs.SmoothKern(graph, "normal", (bin_range / 5.0), n_bins, smooth_x)
        for i in range(n_bins):
            fit_y_binned[i].append(grout.GetPointY(i))

    smooth_y = list()
    smooth_y_up = list()
    smooth_y_down = list()

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
    
    corr_dict = {
            "edges": _bins,
            "nominal": _nom,
            "up": _up,
            "down": _down,
        }

    smooth_graph = ROOT.TGraphAsymmErrors(
        len(smooth_x), smooth_x, smooth_y, 0, 0, smooth_y_down, smooth_y_up
    )
    return smooth_graph, corr_dict
