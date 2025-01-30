import itertools as itt
import logging
from io import StringIO
from typing import Any, Callable, Dict, List, Tuple, Union

import numpy as np
import ROOT
from wurlitzer import STDOUT, pipes


def func_yerr(x, param, n, func):
    grad, eps, func_param = [], 1e-6, [param[i] for i in range(n)]
    cov = [[param[n + i * n + j] for j in range(n)] for i in range(n)]
    for i in range(n):
        func_param1 = [func_param[i] for i in range(n)]
        func_param2 = [func_param[i] for i in range(n)]
        func_param1[i] += eps
        func_param2[i] -= eps
        grad.append((func(x, func_param1) - func(x, func_param2)) / (2 * eps))
    variance = sum(sum(grad[i] * cov[i][j] * grad[j] for j in range(n)) for i in range(n))
    variance = (variance ** (2)) ** 0.5
    return (variance) ** 0.5


def str_func_yerr(x, param, n, func):
    grad, eps, func_param = [], 1e-6, [param[i] for i in range(n)]
    cov = [[param[n + i * n + j] for j in range(n)] for i in range(n)]
    for i in range(n):
        func_param1 = [func_param[i] for i in range(n)]
        func_param2 = [func_param[i] for i in range(n)]
        func_param1[i] += eps
        func_param2[i] -= eps
        grad.append(f"(({func(x, func_param1)} - {func(x, func_param2)}) / (2.0 * {eps}))")
    variance = " + ".join(" + ".join(f"(({grad[i]}) * ({cov[i][j]}) * ({grad[j]}))" for j in range(n)) for i in range(n))
    variance = f"pow(pow({variance}, 2), 0.5)"
    return f"pow({variance}, 0.5)"


def poly_n_func(n: int) -> Tuple[Callable, Callable, Callable]:
    """
    Definition of a polynomial function with degree n and an up and down varied version.

    Args:
        n (int): degree of the polynomial

    Returns:
        Callables of the defined polynomial function
    """
    def nominal(x, param):
        return sum(param[i] * x[0] ** i for i in range(n + 1))

    nominal.__name__ = f"poly_{n}"
    nominal.__n_param__ = n + 1

    def up(x, param):
        nom = nominal(x, [param[i] for i in range(n + 1)])
        unc = func_yerr(x, param, n + 1, nominal)
        return nom + unc / 2.0

    up.__name__ = f"poly_{n}_up"
    up.__n_param__ = (n + 1) + (n + 1) * (n + 1)

    def down(x, param):
        nom = nominal(x, [param[i] for i in range(n + 1)])
        unc = func_yerr(x, param, n + 1, nominal)
        return nom - unc / 2.0

    down.__name__ = f"poly_{n}_down"
    down.__n_param__ = (n + 1) + (n + 1) * (n + 1)

    return nominal, up, down


def poly_n_str_func(n: int) -> Tuple[Callable, Callable, Callable]:
    """
    Definition of a polynomial function as a string with degree n and an up and down varied version.

    Args:
        n (int): degree of the polynomial

    Returns:
        Callables of the defined polynomial string function
    """
    def nominal(x, param):
        expr = " + ".join((f"{param[i]} * pow({x}, {i})" for i in range(n + 1)))
        return f"( {expr} )"

    nominal.__name__ = f"poly_{n}"
    nominal.__n_param__ = n + 1

    def up(x, param):
        nom = nominal(x, [param[i] for i in range(nominal.__n_param__)])
        unc = str_func_yerr(x, param, nominal.__n_param__, nominal)
        return f"( ( {nom} ) + ( {unc} ) / 2.0 )"

    up.__name__ = f"poly_{n}_up"
    up.__n_param__ = (n + 1) + (n + 1) * (n + 1)

    def down(x, param):
        nom = nominal(x, [param[i] for i in range(nominal.__n_param__)])
        unc = str_func_yerr(x, param, nominal.__n_param__, nominal)
        return f"( ( {nom} ) - ( {unc} ) / 2.0 )"

    down.__name__ = f"poly_{n}_down"
    down.__n_param__ = (n + 1) + (n + 1) * (n + 1)

    return nominal, up, down


def extract_param_and_cov(fit: Any) -> Tuple[np.ndarray, np.ndarray]:
    param = np.array(
        [
            fit.Get().Parameter(i)
            for i in range(fit.Get().Parameters().size())
        ],
    )

    n_params = fit.Get().Parameters().size()
    cov, _cov = np.zeros((n_params, n_params)), fit.Get().GetCovarianceMatrix()
    for i, j in itt.product(range(n_params), range(n_params)):
        cov[i, j] = _cov(i, j)

    return param, cov


def check_function_validity_within_bounds(
    bounds: Tuple[float, float],
    funcs: Dict[str, Callable[..., List[float]]],
    parameters: Dict[str, List[float]],
) -> bool:
    x = np.linspace(*bounds, 100)
    keys = [
        "nominal",
        "up",  # TODO: is this needed?
        "down",  # TODO: is this needed?
    ]
    return all(all(funcs[k]([it], parameters[k]) > 0 for it in x) for k in keys)


def get_wrapped_functions_from_fits(
    graph: ROOT.TGraphAsymmErrors,
    bounds: Tuple[float, float],
    ff_hist_up: ROOT.TH1,
    ff_hist_down: ROOT.TH1,
    do_mc_subtr_unc: bool,
    logger: str,
    function_collection: Union[List[Callable], Tuple[Callable, ...]] = (
        "poly_1",
    ),
    verbose: bool = False,
    convert_to_callable: bool = True,
    **kwargs: Dict[str, Any],
) -> Tuple[Dict[str, Union[str, Callable]], str]:
    """
    Fits a set of functions to a given graph and returns functions
    (nominal, nominal + error, nominal - error) for the best fit.
    
    Args:
        graph (ROOT.TGraphAsymmErrors): The graph to fit.
        bounds (Tuple[float, float]): The bounds for the fit.
        ff_hist_up (ROOT.TH1): Histogram for the upward variation of the fit.
        ff_hist_down (ROOT.TH1): Histogram for the downward variation of the fit.
        do_mc_subtr_unc (bool): Whether to consider MC subtraction.
        logger (str): Logger name for logging fit information.
        function_collection (Tuple[str, ...], optional):
            Collection of functions to fit. Defaults to polynomial functions of degree 1 to 5.
            If more than one function is given, the one with the best chi2/ndf is chosen.
            function must be defined `generated_functions` and `generated_str_functions`.
        verbose (bool, optional): Whether to print detailed fit information. Defaults to True.
        convert_to_callable (bool): Boolean to define the return format string or Callable
        
    
    Returns: Dict[str, Callable]: Dictionary with the best fit function containing:
        - "nominal": The best fit function.
        - "up": The best fit function including the upper error.
        - "down": The best fit function including the lower error.
        - "mc_up": The best fit function including the upper error from MC subtraction.
        - "mc_down": The best fit function including the lower error from MC subtraction
    """
    log = logging.getLogger(logger)
    
    a, b = bounds
    TF1s, Fits = dict(), dict()
    best_chi2_over_ndf, name = float("inf"), ""
    TF1s_up, TF1s_down, Fits_up, Fits_down = dict(), dict(), dict(), dict()
    
    generated_functions, generated_str_functions = {}, {}
    generated_functions.update({f"{i}": poly_n_func(int(i.split("_")[1])) for i in function_collection})
    generated_str_functions.update({f"{i}": poly_n_str_func(int(i.split("_")[1])) for i in function_collection})

    for func_collection_name in function_collection:
        func, _, _ = generated_functions[func_collection_name]
        func_name, n_params = func.__name__, func.__n_param__
        TF1s[func_name] = ROOT.TF1(func_name, func, a, b, n_params)
        Fits[func_name] = graph.Fit(TF1s[func_name], "SFN")

        try:
            chi2_over_ndf = Fits[func_name].Chi2() / Fits[func_name].Ndf()
        except ZeroDivisionError:
            chi2_over_ndf = float("inf")

        __param, __cov = extract_param_and_cov(Fits[func_name])
        if not check_function_validity_within_bounds(
            bounds=bounds,
            funcs=dict(zip(["nominal", "up", "down"], generated_functions[func_name])),
            parameters=dict(
                zip(
                    ["nominal", "up", "down"],
                    [__param, [*__param, *__cov.flatten()], [*__param, *__cov.flatten()]],
                )
            ),
        ) and len(function_collection) > 1:
            chi2_over_ndf = float("inf")

        if abs(chi2_over_ndf - 1) < abs(best_chi2_over_ndf - 1):
            best_chi2_over_ndf, name = chi2_over_ndf, func_name

            if do_mc_subtr_unc:
                TF1s_up[func_name] = ROOT.TF1(func_name, func, a, b, n_params)
                TF1s_down[func_name] = ROOT.TF1(func_name, func, a, b, n_params)

                Fits_up[func_name] = ff_hist_up.Fit(TF1s_up[func_name], "SFN")
                Fits_down[func_name] = ff_hist_down.Fit(TF1s_down[func_name], "SFN")

    if best_chi2_over_ndf == float("inf"):
        raise ValueError("No valid function found for the given bounds.")

    if verbose:
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            log.info(f"Using: {name}")
            Fits[name].Print()
            Fits[name].Get().GetCovarianceMatrix().Print()
            log.info("-" * 50)
            if do_mc_subtr_unc:
                log.info("Up")
                Fits_up[name].Print()
                Fits_up[name].Get().GetCovarianceMatrix().Print()
                log.info("-" * 50)
                log.info("Down")
                Fits_down[name].Print()
                Fits_down[name].Get().GetCovarianceMatrix().Print()
                log.info("-" * 50)            
        log.info(out.getvalue())
        log.info("-" * 50)

    results = {}

    param, cov = extract_param_and_cov(Fits[name])
    for _k, _func, _vars in zip(
        ["nominal", "up", "down"],
        (generated_functions if convert_to_callable else generated_str_functions)[name],
        [param, [*param, *cov.flatten()], [*param, *cov.flatten()]],
    ):
        if convert_to_callable:
            results[_k] = ROOT.TF1(f"{_func.__name__}_{_k}", _func, a, b, _func.__n_param__)
            for i, p in enumerate(_vars):
                results[_k].SetParameter(i, p)
        else:
            results[_k] = _func(" ( x ) ", _vars)

    if do_mc_subtr_unc:
        param_up, _ = extract_param_and_cov(Fits_up[name])
        param_down, _ = extract_param_and_cov(Fits_down[name])
        _func, _, _ = (generated_functions if convert_to_callable else generated_str_functions)[name]
        for _k, _vars in zip(
            ("mc_up", "mc_down"),
            (param_up, param_down)
        ):
            if convert_to_callable:
                results[_k] = ROOT.TF1(f"{_func.__name__}_{_k}", _func, a, b, _func.__n_param__)
                for i, p in enumerate(_vars):
                    results[_k].SetParameter(i, p)
            else:
                results[_k] = _func(" ( x ) ", _vars)

    return results, name


def hist_func(
    hist: ROOT.TH1,
    return_callable: bool = True,
) -> Tuple[Callable, List[float], List[float]]:

    _nominal, _up, _down, _edges = [], [], [], []
    for i in range(1, hist.GetNbinsX() + 1):
        _nominal.append(hist.GetBinContent(i))
        _up.append(hist.GetBinContent(i) + hist.GetBinErrorUp(i))
        _down.append(hist.GetBinContent(i) - hist.GetBinErrorLow(i))
        _edges.append((hist.GetBinLowEdge(i), hist.GetBinLowEdge(i + 1)))

    def nominal(x):
        if x < _edges[0][0]:
            return _nominal[0]
        elif x >= _edges[-1][1]:
            return _nominal[-1]
        else:
            for value, window in zip(_nominal, _edges):
                if x >= window[0] and x < window[1]:
                    return value

    def up(x):
        if x < _edges[0][0]:
            return _up[0]
        elif x >= _edges[-1][1]:
            return _up[-1]
        else:
            for value, window in zip(_up, _edges):
                if x >= window[0] and x < window[1]:
                    return value

    def down(x):
        if x < _edges[0][0]:
            return _down[0]
        elif x >= _edges[-1][1]:
            return _down[-1]
        else:
            for value, window in zip(_down, _edges):
                if x >= window[0] and x < window[1]:
                    return value

    return (nominal, up, down) if return_callable else (_nominal, _up, _down)


def get_wrapped_hists(
    ff_hist: ROOT.TH1,
    ff_hist_up: ROOT.TH1,
    ff_hist_down: ROOT.TH1,
    do_mc_subtr_unc: bool,
    logger: str,
    verbose: bool = False, 
    convert_to_callable: bool = True,
    **kwargs: Dict[str, Any],
) -> Tuple[Dict[str, Union[List, Callable]], str]:
    """
    Returns histograms (nominal, nominal + error, nominal - error).
    
    Args:
        ff_hist (ROOT.TH1): Histogram for the nominal variation.
        ff_hist_up (ROOT.TH1): Histogram for the upward variation.
        ff_hist_down (ROOT.TH1): Histogram for the downward variation.
        do_mc_subtr_unc (bool): Whether to consider MC subtraction.
        logger (str): Logger name for logging fit information.
        convert_to_callable (bool): Boolean to define the return format List or Callable
    
    Returns: Dict[str, Union[List, Callable]]: Dictionary with the histograms containing:
        - "nominal": measured histogram.
        - "up": measured histogram including the upper error.
        - "down": measured histogram including the lower error.
        - "mc_up": measured histogram including the upper error from MC subtraction.
        - "mc_down": measured histogram including the lower error from MC subtraction
    """
    log = logging.getLogger(logger)
    if verbose:
        log.info("Measured histograms directly used instead of a fit.")
        log.info("-" * 50)
    
    results = dict(
        zip(
            ["nominal", "up", "down"],
            hist_func(ff_hist, return_callable=convert_to_callable),
        )
    )
    if do_mc_subtr_unc:
        results.update(
            {
                "mc_up": hist_func(ff_hist_up, return_callable=convert_to_callable)[0],
                "mc_down": hist_func(ff_hist_down, return_callable=convert_to_callable)[0],
            }
        )

    return results, "bin_wise"
