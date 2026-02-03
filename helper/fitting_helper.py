import itertools as itt
import logging
from copy import deepcopy
from functools import partial
from io import StringIO
from typing import Any, Callable, Dict, List, Tuple, Union

import numpy as np
import ROOT
from wurlitzer import STDOUT, pipes

import CustomLogging as logging_helper


def _get_gradients(
    x: Union[List[float], str],
    param: List[float],
    n: int,
    func: Callable,
    is_string: bool = False,
) -> Union[List[float], List[str]]:
    grad, eps = [], 1e-6
    for i in range(n):
        param1, param2 = deepcopy(param), deepcopy(param)
        param1[i] += eps
        param2[i] -= eps
        if is_string:
            grad.append(f"(({func(x, param1)} - {func(x, param2)}) / (2.0 * {eps}))")
        else:
            grad.append((func(x, param1) - func(x, param2)) / (2 * eps))
    return grad


def func_yerr(
    x: Union[List[float], str], par: List[float], n: int, func: Callable
) -> float:
    # par contains func parameters and the flattened covariance matrix
    grad = _get_gradients(x, [par[i] for i in range(n)], n, func, is_string=False)
    cov = [[par[n + i * n + j] for j in range(n)] for i in range(n)]
    variance = sum(
        sum(grad[i] * cov[i][j] * grad[j] for j in range(n)) for i in range(n)
    )
    variance = (variance ** (2)) ** 0.5
    return (variance) ** 0.5


def str_func_yerr(
    x: Union[List[float], List[str]], param: List[float], n: int, func: Callable
) -> str:
    # par contains func parameters and the flattened covariance matrix
    grad = _get_gradients(x, [param[i] for i in range(n)], n, func, is_string=True)
    cov = [[param[n + i * n + j] for j in range(n)] for i in range(n)]
    variance = " + ".join(
        " + ".join(f"(({grad[i]}) * ({cov[i][j]}) * ({grad[j]}))" for j in range(n))
        for i in range(n)
    )
    variance = f"pow(pow({variance}, 2), 0.5)"
    return f"pow({variance}, 0.5)"


def _limit(
    value: Union[float, str],
    *boundaries: float,
    is_string: bool = False,
) -> Union[float, str]:
    assert (
        len(boundaries) == 2
    ), "Boundaries must be a tuple of two (lower, upper) float values."
    minimum, maximum = boundaries
    if is_string:
        minimum = "-9999.0" if minimum == -float("inf") else str(float(minimum))
        maximum = "9999.0" if maximum == float("inf") else str(float(maximum))
        return f"(max(min({value}, {maximum}), {minimum}))"
    else:
        return max(min(value, maximum), minimum)


_no_limit_default = (-float("inf"), float("inf"))


def poly_n_func(
    n: int,
    #
    limit_x_nominal: Tuple[float, float] = _no_limit_default,
    limit_x_up: Tuple[float, float] = _no_limit_default,
    limit_x_down: Tuple[float, float] = _no_limit_default,
    #
    limit_y_nominal: Tuple[float, float] = _no_limit_default,
    limit_y_up: Tuple[float, float] = _no_limit_default,
    limit_y_down: Tuple[float, float] = _no_limit_default,
) -> Tuple[Callable, Callable, Callable]:
    """
    Definition of a polynomial function with degree n and an up and down varied version.

    Args:
        n (int): degree of the polynomial

    Returns:
        Callables of the defined polynomial function
    """

    def nominal(x, param):
        x = [_limit(x[0], *limit_x_nominal)]
        nom = sum(param[i] * x[0] ** i for i in range(n + 1))
        return _limit(nom, *limit_y_nominal)

    nominal.__name__ = f"poly_{n}"
    nominal.__n_param__ = n + 1

    def up(x, param):
        # due to ROOT it cannot be defined in outer scope
        def _nominal(x, param):
            return sum(param[i] * x[0] ** i for i in range(n + 1))

        x = [_limit(x[0], *limit_x_up)]
        nom = nominal(x, [param[i] for i in range(nominal.__n_param__)])
        unc = func_yerr(x, param, nominal.__n_param__, _nominal)
        return _limit(nom + unc / 2.0, *limit_y_up)

    up.__name__ = f"poly_{n}_up"
    up.__n_param__ = (n + 1) + (n + 1) * (n + 1)

    def down(x, param):
        # due to ROOT it cannot be defined in outer scope
        def _nominal(x, param):
            return sum(param[i] * x[0] ** i for i in range(n + 1))

        x = [_limit(x[0], *limit_x_down)]
        nom = nominal(x, [param[i] for i in range(nominal.__n_param__)])
        unc = func_yerr(x, param, nominal.__n_param__, _nominal)
        return _limit(nom - unc / 2.0, *limit_y_down)

    down.__name__ = f"poly_{n}_down"
    down.__n_param__ = (n + 1) + (n + 1) * (n + 1)

    return nominal, up, down


def poly_n_str_func(
    n: int,
    #
    limit_x_nominal: Tuple[float, float] = _no_limit_default,
    limit_x_up: Tuple[float, float] = _no_limit_default,
    limit_x_down: Tuple[float, float] = _no_limit_default,
    #
    limit_y_nominal: Tuple[float, float] = _no_limit_default,
    limit_y_up: Tuple[float, float] = _no_limit_default,
    limit_y_down: Tuple[float, float] = _no_limit_default,
) -> Tuple[Callable, Callable, Callable]:
    """
    Definition of a polynomial function as a string with degree n and an up and down varied version.

    Args:
        n (int): degree of the polynomial

    Returns:
        Callables of the defined polynomial string function
    """

    def _nominal(x, param):
        nom = " + ".join((f"{param[i]} * pow({x}, {i})" for i in range(n + 1)))
        return f"({nom})"

    def nominal(x, param):
        x = _limit(x, *limit_x_nominal, is_string=True)
        nom = " + ".join((f"{param[i]} * pow({x}, {i})" for i in range(n + 1)))
        return _limit(f"({nom})", *limit_y_nominal, is_string=True)

    nominal.__name__ = f"poly_{n}"
    nominal.__n_param__ = n + 1

    def up(x, param):
        x = _limit(x, *limit_x_up, is_string=True)
        nom = nominal(x, [param[i] for i in range(nominal.__n_param__)])
        unc = str_func_yerr(x, param, nominal.__n_param__, _nominal)
        return _limit(f"(({nom}) + ({unc}) / 2.0)", *limit_y_up, is_string=True)

    up.__name__ = f"poly_{n}_up"
    up.__n_param__ = (n + 1) + (n + 1) * (n + 1)

    def down(x, param):
        x = _limit(x, *limit_x_down, is_string=True)
        nom = nominal(x, [param[i] for i in range(nominal.__n_param__)])
        unc = str_func_yerr(x, param, nominal.__n_param__, _nominal)
        return _limit(f"(({nom}) - ({unc}) / 2.0)", *limit_y_down, is_string=True)

    down.__name__ = f"poly_{n}_down"
    down.__n_param__ = (n + 1) + (n + 1) * (n + 1)

    return nominal, up, down


def extract_param_and_cov(fit: Any) -> Tuple[np.ndarray, np.ndarray]:
    param = np.array(
        [fit.Get().Parameter(i) for i in range(fit.Get().Parameters().size())],
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


def get_default_limit_kwargs(
    logger: str,
    bounds: Tuple[float, float],
    limit_x: Union[
        Tuple[float, float], Dict[str, Tuple[float, float]]
    ] = _no_limit_default,
    limit_y: Union[
        Tuple[float, float], Dict[str, Tuple[float, float]]
    ] = _no_limit_default,
    verbose: bool = True,
) -> Dict[str, Union[Tuple[float, float], Tuple[None, ...]]]:
    kwargs_dict = dict(
        limit_x_nominal=(
            limit_x
            if isinstance(limit_x, tuple)
            else limit_x.get("nominal", _no_limit_default)
        ),
        limit_x_up=(
            limit_x
            if isinstance(limit_x, tuple)
            else limit_x.get("up", _no_limit_default)
        ),
        limit_x_down=(
            limit_x
            if isinstance(limit_x, tuple)
            else limit_x.get("down", _no_limit_default)
        ),
        limit_y_nominal=(
            limit_y
            if isinstance(limit_y, tuple)
            else limit_y.get("nominal", _no_limit_default)
        ),
        limit_y_up=(
            limit_y
            if isinstance(limit_y, tuple)
            else limit_y.get("up", _no_limit_default)
        ),
        limit_y_down=(
            limit_y
            if isinstance(limit_y, tuple)
            else limit_y.get("down", _no_limit_default)
        ),
    )

    if isinstance(limit_x, tuple) and limit_x == _no_limit_default:
        msg = " ".join(
            [
                f"No custom limit_x was provided, using provided boundaries {bounds} for ",
                "nominal, up and down functions",
            ],
        )
        if verbose:
            logger.info(msg)

        kwargs_dict.update(
            dict(
                limit_x_nominal=bounds,
                limit_x_up=bounds,
                limit_x_down=bounds,
            ),
        )
    if isinstance(limit_y, tuple) and limit_y == _no_limit_default:
        msg = " ".join(
            [
                "No custom limit_y was provided, using default of (0.0, float('inf')) for ",
                "nominal, up and down functions",
            ],
        )
        if verbose:
            logger.info(msg)

        kwargs_dict.update(
            dict(
                limit_y_nominal=(0.0, float("inf")),
                limit_y_up=(0.0, float("inf")),
                limit_y_down=(0.0, float("inf")),
            ),
        )
    return kwargs_dict


def get_wrapped_functions_from_fits(
    bounds: Tuple[float, float],
    ff_hist: ROOT.TGraphAsymmErrors,
    ff_hist_up: ROOT.TGraphAsymmErrors,
    ff_hist_down: ROOT.TGraphAsymmErrors,
    do_mc_subtr_unc: bool,
    logger: str,
    function_collection: Union[List[Callable], Tuple[Callable, ...]] = ("poly_1",),
    verbose: bool = True,
    #
    limit_x: Union[
        Tuple[float, float], Dict[str, Tuple[float, float]]
    ] = _no_limit_default,
    limit_y: Union[
        Tuple[float, float], Dict[str, Tuple[float, float]]
    ] = _no_limit_default,
    #
    **kwargs: Dict[str, Any],
) -> Tuple[Dict[str, Callable], Dict[str, str]]:
    """
    Fits a set of functions to a given graph and returns functions
    (nominal, nominal + error, nominal - error) for the best fit.

    Args:
        bounds (Tuple[float, float]): The bounds for the fit.
        ff_hist (ROOT.TGraphAsymmErrors): Histogram for the nominal fit.
        ff_hist_up (ROOT.TGraphAsymmErrors): Histogram for the upward variation of the fit.
        ff_hist_down (ROOT.TGraphAsymmErrors): Histogram for the downward variation of the fit.
        do_mc_subtr_unc (bool): Whether to consider MC subtraction.
        logger (str): Logger name for logging fit information.
        function_collection (Tuple[str, ...], optional):
            Collection of functions to fit. Defaults to polynomial functions of degree 1 to 5.
            If more than one function is given, the one with the best chi2/ndf is chosen.
        verbose (bool, optional): Whether to print detailed fit information. Defaults to True.

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

    limit_kwargs = get_default_limit_kwargs(
        logger=log,
        bounds=bounds,
        limit_x=limit_x,
        limit_y=limit_y,
        verbose=verbose,
    )

    poly_n_func_limited = partial(poly_n_func, **limit_kwargs)
    poly_n_str_func_limited = partial(poly_n_str_func, **limit_kwargs)

    _functions, _str_functions = dict(), dict()
    _functions_limited, _str_functions_limited = dict(), dict()
    for func in function_collection:
        if "poly_" in func:
            n = int(func.split("_")[-1])
            _functions[func] = poly_n_func(n)
            _str_functions[func] = poly_n_str_func(n)
            _functions_limited[func] = poly_n_func_limited(n)
            _str_functions_limited[func] = poly_n_str_func_limited(n)
        else:
            raise ValueError(f"Unknown/unsupported function: {func}")

    for func_collection_name in function_collection:
        func, _, _ = _functions[func_collection_name]
        func_name, n_params = func.__name__, func.__n_param__
        TF1s[func_name] = ROOT.TF1(func_name, func, a, b, n_params)
        Fits[func_name] = ff_hist.Fit(TF1s[func_name], "SFN")

        try:
            chi2_over_ndf = Fits[func_name].Chi2() / Fits[func_name].Ndf()
        except ZeroDivisionError:
            chi2_over_ndf = float("inf")

        __param, __cov = extract_param_and_cov(Fits[func_name])
        if (
            not check_function_validity_within_bounds(
                funcs=dict(zip(["nominal", "up", "down"], _functions[func_name])),
                bounds=bounds,
                parameters=dict(
                    zip(
                        ["nominal", "up", "down"],
                        [
                            __param,
                            [*__param, *__cov.flatten()],
                            [*__param, *__cov.flatten()],
                        ],
                    )
                ),
            ) and len(function_collection) > 1
        ):
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

    callable_results, str_results = {}, {}

    nominal_param, cov = extract_param_and_cov(Fits[name])
    variation_param = [*nominal_param, *cov.flatten()]

    callable_results["nominal"] = lambda x: _functions_limited[name][0](
        [x], nominal_param
    )
    callable_results["unc_up"] = lambda x: _functions_limited[name][1](
        [x], variation_param
    )
    callable_results["unc_down"] = lambda x: _functions_limited[name][2](
        [x], variation_param
    )

    str_results["nominal"] = _str_functions_limited[name][0](" ( x ) ", nominal_param)
    str_results["unc_up"] = _str_functions_limited[name][1](" ( x ) ", variation_param)
    str_results["unc_down"] = _str_functions_limited[name][2](
        " ( x ) ", variation_param
    )

    if do_mc_subtr_unc:
        param_up, _ = extract_param_and_cov(Fits_up[name])
        param_down, _ = extract_param_and_cov(Fits_down[name])
        callable_results["mc_subtraction_unc_up"] = lambda x: _functions_limited[name][0](
            [x],
            param_up,
        )
        callable_results["mc_subtraction_unc_down"] = lambda x: _functions_limited[name][0](
            [x],
            param_down,
        )
        str_results["mc_subtraction_unc_up"] = _str_functions_limited[name][0](
            " ( x ) ",
            param_up,
        )
        str_results["mc_subtraction_unc_down"] = _str_functions_limited[name][0](
            " ( x ) ",
            param_down,
        )

    return callable_results, str_results, name


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
    verbose: bool = True,
    **kwargs: Dict[str, Any],
) -> Dict[str, Union[ROOT.TF1, str, Callable]]:
    """
    Returns histograms (nominal, nominal + error, nominal - error).

    Args:
        ff_hist (ROOT.TH1): Histogram for the nominal variation.
        ff_hist_up (ROOT.TH1): Histogram for the upward variation.
        ff_hist_down (ROOT.TH1): Histogram for the downward variation.
        do_mc_subtr_unc (bool): Whether to consider MC subtraction.
        logger (str): Logger name for logging fit information.

    Returns: Dict[str, Union[List, Callable]]: Dictionary with the histograms containing:
        - "nominal": measured histogram.
        - "up": measured histogram including the upper error.
        - "down": measured histogram including the lower error.
        - "mc_up": measured histogram including the upper error from MC subtraction.
        - "mc_down": measured histogram including the lower error from MC subtraction
    """
    log = logging.getLogger(logger)
    if verbose:
        log.info("Measured histograms directly instead of a fit.")
        log.info("-" * 50)

    callable_results = dict(
        zip(
            ["nominal", "unc_up", "unc_down"],
            hist_func(ff_hist, return_callable=True),
        )
    )
    str_results = dict(
        zip(
            ["nominal", "unc_up", "unc_down"],
            hist_func(ff_hist, return_callable=False),
        )
    )
    if do_mc_subtr_unc:
        callable_results.update(
            {
                "mc_subtraction_unc_up": hist_func(ff_hist_up, return_callable=True)[0],
                "mc_subtraction_unc_down": hist_func(
                    ff_hist_down, return_callable=True
                )[0],
            }
        )
        str_results.update(
            {
                "mc_subtraction_unc_up": hist_func(ff_hist_up, return_callable=False)[
                    0
                ],
                "mc_subtraction_unc_down": hist_func(
                    ff_hist_down, return_callable=False
                )[0],
            }
        )

    return callable_results, str_results, str("binwise")
