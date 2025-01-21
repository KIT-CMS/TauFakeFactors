import itertools as itt
import logging
from io import StringIO
from typing import Any, Callable, Dict, List, Literal, Tuple, Union

import numpy as np
import ROOT
from wurlitzer import STDOUT, pipes


def func_yerr(x, par, n, func):
    grad, eps, func_par = [], 1e-6, [par[i] for i in range(n)]
    cov = [[par[n + i * n + j] for j in range(n)] for i in range(n)]
    for i in range(n):
        func_par1 = [func_par[i] for i in range(n)]
        func_par2 = [func_par[i] for i in range(n)]
        func_par1[i] += eps
        func_par2[i] -= eps
        grad.append((func(x, func_par1) - func(x, func_par2)) / (2 * eps))
    variance = sum(sum(grad[i] * cov[i][j] * grad[j] for j in range(n)) for i in range(n))
    variance = (variance ** (2)) ** 0.5
    return (variance) ** 0.5


def str_func_yerr(x, par, n, func):
    grad, eps, func_par = [], 1e-6, [par[i] for i in range(n)]
    cov = [[par[n + i * n + j] for j in range(n)] for i in range(n)]
    for i in range(n):
        func_par1 = [func_par[i] for i in range(n)]
        func_par2 = [func_par[i] for i in range(n)]
        func_par1[i] += eps
        func_par2[i] -= eps
        grad.append(f"(({func(x, func_par1)} - {func(x, func_par2)}) / (2.0 * {eps}))")
    variance = " + ".join(" + ".join(f"(({grad[i]}) * ({cov[i][j]}) * ({grad[j]}))" for j in range(n)) for i in range(n))
    variance = f"pow(pow({variance}, 2), 0.5)"
    return f"pow({variance}, 0.5)"


def poly_n_func(n: int) -> Tuple[Callable, Callable]:

    def nominal(x, par):
        return sum(par[i] * x[0] ** i for i in range(n + 1))

    nominal.__name__ = f"poly_{n}"
    nominal.__n_par__ = n + 1

    def up(x, par):
        nom = nominal(x, [par[i] for i in range(n + 1)])
        unct = func_yerr(x, par, n + 1, nominal)
        return nom + unct / 2.0

    up.__name__ = f"poly_{n}_up"
    up.__n_par__ = (n + 1) + (n + 1) * (n + 1)

    def down(x, par):
        nom = nominal(x, [par[i] for i in range(n + 1)])
        unct = func_yerr(x, par, n + 1, nominal)
        return nom - unct / 2.0

    down.__name__ = f"poly_{n}_down"
    down.__n_par__ = (n + 1) + (n + 1) * (n + 1)

    return nominal, up, down


def poly_n_str_func(n: int) -> Tuple[Callable, Callable]:

    def nominal(x, par_nominal):
        expr = " + ".join((f"{par_nominal[i]} * pow({x}, {i})" for i in range(n + 1)))
        return f"( {expr} )"

    nominal.__name__ = f"poly_{n}"
    nominal.__n_par__ = n + 1

    def up(x, par_up):
        _nom_up = nominal(x, [par_up[i] for i in range(nominal.__n_par__)])
        _up = str_func_yerr(x, par_up, nominal.__n_par__, nominal)
        return f"( ( {_nom_up} ) + ( {_up} ) / 2.0 )"

    up.__name__ = f"poly_{n}_up"
    up.__n_par__ = (n + 1) + (n + 1) * (n + 1)

    def down(x, par_down):
        _nom_down = nominal(x, [par_down[i] for i in range(nominal.__n_par__)])
        _down = str_func_yerr(x, par_down, nominal.__n_par__, nominal)
        return f"( ( {_nom_down} ) - ( {_down} ) / 2.0 )"

    down.__name__ = f"poly_{n}_down"
    down.__n_par__ = (n + 1) + (n + 1) * (n + 1)

    return nominal, up, down


def extract_par_and_cov(fit: Any) -> Tuple[np.ndarray, np.ndarray]:
    par = np.array(
        [
            fit.Get().Parameter(i)
            for i in range(fit.Get().Parameters().size())
        ],
    )

    N = fit.Get().Parameters().size()
    cov, _cov = np.zeros((N, N)), fit.Get().GetCovarianceMatrix()
    for i, j in itt.product(range(N), range(N)):
        cov[i, j] = _cov(i, j)

    return par, cov


generated_functions, generated_str_functions = {}, {}
generated_functions.update({f"poly_{i}": poly_n_func(i) for i in range(0, 6)})
generated_str_functions.update({f"poly_{i}": poly_n_str_func(i) for i in range(0, 6)})


def check_function_validity_within_bounds(
    funcs: Dict[str, Callable[..., List[float]]],
    bounds: Tuple[float, float],
    parameters: Dict[str, List[float]],
) -> bool:
    x = np.linspace(*bounds, 100)
    keys = [
        "nominal",
        # "up",  # TODO: is this needed?
        # "down",  # TODO: is this needed?
    ]
    return all(all(funcs[k]([it], parameters[k]) > 0 for it in x) for k in keys)


def get_wrapped_functions_from_fits(
    graph: ROOT.TGraphAsymmErrors,
    bounds: Tuple[float, float],
    do_mc_subtr_unc: bool,
    ff_hist_up: ROOT.TH1,
    ff_hist_down: ROOT.TH1,
    function_collection: Union[List[Callable], Tuple[Callable, ...]] = (
        "poly_1",
        "poly_2",
        "poly_3",
        "poly_4",
        "poly_5",
    ),
    verbose: bool = True,
    logger: Union[str, None] = None,
    convert_to_callable: bool = True,
    **kwargs: Dict[str, Any],
) -> Dict[str, Union[ROOT.TF1, str, Callable]]:
    """
    Fits a set of functions to a given graph and returns functions
    (nominal, nominal + error, nominal - error) for the best fit.
    Args:
        graph (ROOT.TGraphAsymmErrors): The graph to fit.
        bounds (Tuple[float, float]): The bounds for the fit.
        do_mc_subtr_unc (bool): Whether to consider MC subtraction.
        ff_hist_up (ROOT.TH1): Histogram for the upward variation of the fit.
        ff_hist_down (ROOT.TH1): Histogram for the downward variation of the fit.
        function_collection (Tuple[str, ...], optional):
            Collection of functions to fit. Defaults to polynomial functions of degree 1 to 5.
            If more than one function is given, the one with the best chi2/ndf is chosen.
            function must be defined `generated_functions` and `generated_str_functions`.
        verbose (bool, optional): Whether to print detailed fit information. Defaults to True.
        root_wrap (bool, optional): Whether to wrap the results in ROOT.TF1 objects. Defaults to True.
        logger (Union[str, None], optional): Logger name for logging fit information. Defaults to None.
    Returns: Dict[str, ROOT.TF1]: Dictionary with the best fit function containing:
        - "nominal": The best fit function.
        - "up": The best fit function including the upper error.
        - "down": The best fit function including the lower error.
        - "mc_up": The best fit function including the upper error from MC subtraction.
        - "mc_down": The best fit function including the lower error from MC subtraction
    """
    a, b = bounds
    TF1s, Fits = dict(), dict()
    best_chi2_over_ndf, name = float("inf"), ""
    TF1s_up, TF1s_down, Fits_up, Fits_down = dict(), dict(), dict(), dict()

    for func_collection_name in function_collection:
        func, _, _ = generated_functions[func_collection_name]
        func_name, n_pars = func.__name__, func.__n_par__
        TF1s[func_name] = ROOT.TF1(func_name, func, a, b, n_pars)
        Fits[func_name] = graph.Fit(TF1s[func_name], "SFN")

        try:
            chi2_over_ndf = Fits[func_name].Chi2() / Fits[func_name].Ndf()
        except ZeroDivisionError:
            chi2_over_ndf = float("inf")

        __par, __cov = extract_par_and_cov(Fits[func_name])
        if not check_function_validity_within_bounds(
            funcs=dict(zip(["nominal", "up", "down"], generated_functions[func_name])),
            bounds=bounds,
            parameters=dict(
                zip(
                    ["nominal", "up", "down"],
                    [__par, [*__par, *__cov.flatten()], [*__par, *__cov.flatten()]],
                )
            ),
        ) and len(function_collection) > 1:
            chi2_over_ndf = float("inf")

        if abs(chi2_over_ndf - 1) < abs(best_chi2_over_ndf - 1):
            best_chi2_over_ndf, name = chi2_over_ndf, func_name

            if do_mc_subtr_unc:
                TF1s_up[func_name] = ROOT.TF1(func_name, func, a, b, n_pars)
                TF1s_down[func_name] = ROOT.TF1(func_name, func, a, b, n_pars)

                Fits_up[func_name] = ff_hist_up.Fit(TF1s_up[func_name], "SFN")
                Fits_down[func_name] = ff_hist_down.Fit(TF1s_down[func_name], "SFN")

    if best_chi2_over_ndf == float("inf"):
        raise ValueError("No valid function found for the given bounds.")

    if verbose:
        out = StringIO()
        with pipes(stdout=out, stderr=STDOUT):
            print(f"\n\nUsing: {name}\n")
            Fits[name].Print()
            Fits[name].Get().GetCovarianceMatrix().Print()
            print("-" * 50)
            if do_mc_subtr_unc:
                print("Up\n")
                Fits_up[name].Print()
                Fits_up[name].Get().GetCovarianceMatrix().Print()
                print("-" * 50)
                print("Down\n")
                Fits_down[name].Print()
                Fits_down[name].Get().GetCovarianceMatrix().Print()
                print("-" * 50)

        print(out.getvalue())
        if logger is not None:
            log = logging.getLogger(logger)
            log.info(out.getvalue())
            log.info("-" * 50)

    results = {}

    par, cov = extract_par_and_cov(Fits[name])
    for _k, _func, _vars in zip(
        ["nominal", "up", "down"],
        (generated_functions if convert_to_callable == convert_to_callable else generated_str_functions)[name],
        [par, [*par, *cov.flatten()], [*par, *cov.flatten()]],
    ):
        if convert_to_callable == convert_to_callable:
            results[_k] = ROOT.TF1(f"{_func.__name__}_{_k}", _func, a, b, _func.__n_par__)
            for i, p in enumerate(_vars):
                results[_k].SetParameter(i, p)
        else:
            results[_k] = _func(" ( x ) ", _vars)

    if do_mc_subtr_unc:
        par_up, _ = extract_par_and_cov(Fits_up[name])
        par_down, _ = extract_par_and_cov(Fits_down[name])
        _func, _, _ = (generated_functions if convert_to_callable == convert_to_callable else generated_str_functions)[name]
        for _k, _vars in zip(
            ("mc_up", "mc_down"),
            (par_up, par_down)
        ):
            if convert_to_callable == convert_to_callable:
                results[_k] = ROOT.TF1(f"{_func.__name__}_{_k}", _func, a, b, _func.__n_par__)
                for i, p in enumerate(_vars):
                    results[_k].SetParameter(i, p)
            else:
                results[_k] = _func(" ( x ) ", _vars)

    return results
