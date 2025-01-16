import itertools as itt
import logging
from abc import ABC, ABCMeta, abstractmethod
from io import StringIO
from typing import Any, Callable, Dict, Tuple, Union, Literal

import numdifftools as nd
import numpy as np
from wurlitzer import STDOUT, pipes

try:
    import ROOT
except ImportError:
    print("ROOT not available, functions including ROOT will not work.")

    class FakeROOT:
        TF1 = None
        TGraphAsymmErrors = None
        TH1 = None
    ROOT = FakeROOT


class _StaticMethodMeta(ABCMeta):
    def __getattr__(cls, name):
        if name.startswith('__') and name.endswith('__'):
            raise AttributeError(f"{cls.__name__} has no attribute {name}")
        func = cls._generate_nth_expression(name)
        setattr(cls, name, staticmethod(func))
        return func


class _Base(ABC, metaclass=_StaticMethodMeta):
    @classmethod
    @abstractmethod
    def _generate_nth_expression(cls, name):
        pass

    @staticmethod
    @abstractmethod
    def _nth_expression(x, par, n):
        pass


class FitFunctions(_Base):
    @staticmethod
    def _nth_expression(x, par, n):
        return sum(par[i] * x[0] ** i for i in range(n))

    @classmethod
    def _generate_nth_expression(cls, name):
        n = int(name.split("_")[1])

        def poly_method(x, par):
            return cls._nth_expression(x, par, n + 1)

        poly_method.__name__ = name
        poly_method.__n_par__ = int(n + 1)

        return poly_method


class FitFunctionsStr(_Base):
    @staticmethod
    def _nth_expression(x, par, n):
        expr = " + ".join(f"( {par[i]} * ( x ) ** ({i}) )" for i in range(n))
        return f"( {expr} )"

    @classmethod
    def _generate_nth_expression(cls, name):
        n = int(name.split("_")[1])

        def poly_method(x, par):
            return cls._nth_expression(x, par, n + 1)

        poly_method.__name__ == name

        return poly_method


class FitFunctionsYError(_Base):
    @staticmethod
    def _nth_expression(x, par, n):
        g, h, f = [], 1e-6, getattr(FitFunctions, f"poly_{int(n - 1)}")
        p = [par[i] for i in range(n)]
        c = [[par[n + i * n + j] for j in range(n)] for i in range(n)]
        for i in range(n):
            p1, p2 = [p[i] for i in range(n)], [p[i] for i in range(n)]
            p1[i] += h
            p2[i] -= h
            g.append((f(x, p1) - f(x, p2)) / (2 * h))
        u = sum(sum(g[i] * c[i][j] * g[j] for j in range(n)) for i in range(n))
        u = (u ** (2)) ** 0.5
        return (u) ** 0.5

    @classmethod
    def _generate_nth_expression(cls, name):
        n = int(name.split("_")[1])

        def poly_method(x, par):
            return cls._nth_expression(x, par, n + 1)

        poly_method.__name__ == name

        return poly_method


class FitFunctionsYErrorStr(_Base):
    @staticmethod
    def _nth_expression(x, par, n):
        g, h, f = [], 1e-6, getattr(FitFunctionsStr, f"poly_{int(n - 1)}")
        p = [par[i] for i in range(n)]
        c = [[par[n + i * n + j] for j in range(n)] for i in range(n)]
        for i in range(n):
            p1, p2 = [p[i] for i in range(n)], [p[i] for i in range(n)]
            p1[i] += h
            p2[i] -= h
            g.append(f"( ( {f(x, p1)} - {f(x, p2)} ) / ( 2.0 * {h} ) )")
        u = " + ".join(" + ".join(f"( ( {g[i]} ) * ( {c[i][j]} ) * ( {g[j]} ) )" for j in range(n)) for i in range(n))
        u = f"( ( ( {u} ) ** (2) ) ) ** (0.5)"
        return f"( ( {u} ) ** (0.5) )"

    @classmethod
    def _generate_nth_expression(cls, name):
        n = int(name.split("_")[1])

        def poly_method(x, par):
            return cls._nth_expression(x, par, n + 1)

        poly_method.__name__ == name

        return poly_method


class FitFunctionsUp(_Base):
    @staticmethod
    def _nth_expression(x, par, n):
        _par = [par[i] for i in range(n + 1)]
        return getattr(FitFunctions, f"poly_{int(n)}")(x, _par) + getattr(FitFunctionsYError, f"poly_{int(n)}")(x, par) / 2.0

    @classmethod
    def _generate_nth_expression(cls, name):
        n = int(name.split("_")[1])

        def poly_method(x, par):
            return cls._nth_expression(x, par, n)

        poly_method.__name__ = name
        poly_method.__n_par__ = int((n + 1) + (n + 1) * (n + 1))

        return poly_method


class FitFunctionsDown(_Base):
    @staticmethod
    def _nth_expression(x, par, n):
        _par = [par[i] for i in range(n + 1)]
        return getattr(FitFunctions, f"poly_{int(n)}")(x, _par) - getattr(FitFunctionsYError, f"poly_{int(n)}")(x, par) / 2.0

    @classmethod
    def _generate_nth_expression(cls, name):
        n = int(name.split("_")[1])

        def poly_method(x, par):
            return cls._nth_expression(x, par, n)

        poly_method.__name__ = name
        poly_method.__n_par__ = int((n + 1) + (n + 1) * (n + 1))

        return poly_method


class FitFunctionsUpStr(_Base):
    @staticmethod
    def _nth_expression(x, par, n):
        _par = [par[i] for i in range(n + 1)]
        f_nominal = getattr(FitFunctionsStr, f"poly_{int(n)}")(x, _par)
        f_error = getattr(FitFunctionsYErrorStr, f"poly_{int(n)}")(x, par)
        return f"( ( {f_nominal} ) + ( {f_error} ) / 2.0 )"

    @classmethod
    def _generate_nth_expression(cls, name):
        n = int(name.split("_")[1])

        def poly_method(x, par):
            return cls._nth_expression(x, par, n)

        poly_method.__name__ = name

        return poly_method


class FitFunctionsDownStr(_Base):
    @staticmethod
    def _nth_expression(x, par, n):
        _par = [par[i] for i in range(n + 1)]
        f_nominal = getattr(FitFunctionsStr, f"poly_{int(n)}")(x, _par)
        f_error = getattr(FitFunctionsYErrorStr, f"poly_{int(n)}")(x, par)
        return f"( ( {f_nominal} ) - ( {f_error} ) / 2.0 )"

    @classmethod
    def _generate_nth_expression(cls, name):
        n = int(name.split("_")[1])

        def poly_method(x, par):
            return cls._nth_expression(x, par, n)

        poly_method.__name__ = name

        return poly_method


def _extract_par(fit: Any) -> np.ndarray:
    return np.array(
        [
            fit.Get().Parameter(i)
            for i in range(fit.Get().Parameters().size())
        ],
    )


def _extract_cov(fit: Any) -> np.ndarray:
    N = fit.Get().Parameters().size()
    cov, _cov = np.zeros((N, N)), fit.Get().GetCovarianceMatrix()
    for i, j in itt.product(range(N), range(N)):
        cov[i, j] = _cov(i, j)
    return cov


def extract_par_and_cov(fit: Any) -> Tuple[np.ndarray, np.ndarray]:
    return _extract_par(fit), _extract_cov(fit)


def get_wrapped_error_funcs(func, par, cov):
    def upper_value(x):
        return func(x, par) + np.sqrt(
            np.apply_along_axis(
                lambda J: J @ cov @ J.T,
                1,
                np.atleast_2d(nd.Gradient(lambda par: func(x, par))(par)),
            ),
        ) / 2

    def lower_value(x):
        return func(x, par) - np.sqrt(
            np.apply_along_axis(
                lambda J: J @ cov @ J.T,
                1,
                np.atleast_2d(nd.Gradient(lambda par: func(x, par))(par)),
            ),
        ) / 2

    return upper_value, lower_value


def get_wrapped_func(func, par):
    def _f(x):
        return func(x, par)
    return _f


def get_wrapped_functions_from_fits(
    graph: ROOT.TGraphAsymmErrors,
    bounds: Tuple[float, float],
    do_mc_subtr_unc: bool,
    ff_hist_up: ROOT.TH1,
    ff_hist_down: ROOT.TH1,
    function_collection: Tuple[Callable, ...] = (
        "poly_1",
        "poly_2",
        "poly_3",
        "poly_4",
        "poly_5",
    ),
    verbose: bool = True,
    logger: Union[str, None] = None,
    convert_to: Literal["ROOT", "str"] = "ROOT",
) -> Dict[str, ROOT.TF1]:
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
            function must be defined in termy of FitFunctions (poly_<N>).
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

    for _func_name in function_collection:
        _func = getattr(FitFunctions, _func_name)
        _name, _n_par = _func.__name__, _func.__n_par__
        TF1s[_name] = ROOT.TF1(_name, _func, a, b, _n_par)
        Fits[_name] = graph.Fit(TF1s[_name], "SFN")

        try:
            chi2_over_ndf = Fits[_name].Chi2() / Fits[_name].Ndf()
        except ZeroDivisionError:
            chi2_over_ndf = float("inf")
        if abs(chi2_over_ndf - 1) < abs(best_chi2_over_ndf - 1):
            best_chi2_over_ndf, name = chi2_over_ndf, _name

            if do_mc_subtr_unc:
                TF1s_up[_name] = ROOT.TF1(_name, _func, a, b, _n_par)
                TF1s_down[_name] = ROOT.TF1(_name, _func, a, b, _n_par)

                Fits_up[_name] = ff_hist_up.Fit(TF1s_up[_name], "SFN")
                Fits_down[_name] = ff_hist_down.Fit(TF1s_down[_name], "SFN")

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
    class_collection = (
        FitFunctions,
        FitFunctionsUp,
        FitFunctionsDown,
    ) if convert_to == "ROOT" else (
        FitFunctionsStr,
        FitFunctionsUpStr,
        FitFunctionsDownStr,
    )
    for _k, _cls, _var in zip(
        ["nominal", "up", "down"],
        class_collection,
        [par, [*par, *cov.flatten()], [*par, *cov.flatten()]],
    ):
        _func = getattr(_cls, name)
        if convert_to == "ROOT":
            results[_k] = ROOT.TF1(f"{_func.__name__}_{_k}", _func, a, b, _func.__n_par__)
            for i, p in enumerate(_var):
                results[_k].SetParameter(i, p)
        else:
            results[_k] = lambda x: _func("x", _var)

    if do_mc_subtr_unc:
        par_up, cov_up = extract_par_and_cov(Fits_up[name])
        par_down, cov_down = extract_par_and_cov(Fits_down[name])

        for _k1, _params in zip(
            ("mc_up", "mc_down"),
            (
                (par_up, [*par_up, *cov_up.flatten()], [*par_up, *cov_up.flatten()]),
                (par_down, [*par_down, *cov_down.flatten()], [*par_down, *cov_down.flatten()]),
            )
        ):
            for _k2, _cls, _var in zip(
                ("nominal", "up", "down"),
                class_collection,
                _params,
            ):
                if _k is None:
                    continue

                _func = getattr(_cls, name)

                if convert_to == "ROOT":
                    results[f"{_k1}_{_k2}"] = ROOT.TF1(f"{_func.__name__}_{_k}", _func, a, b, _func.__n_par__)
                    for i, p in enumerate(_var):
                        results[f"{_k1}_{_k2}"].SetParameter(i, p)
                else:
                    results[f"{_k1}_{_k2}"] = lambda x: _func("x", _var)

    return results
