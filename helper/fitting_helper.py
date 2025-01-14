import numpy as np
import itertools as itt
import numdifftools as nd
from typing import Any, Callable, Tuple
try:
    import ROOT
except ImportError:
    class FakeROOT:
        TF1 = None
        TGraphAsymmErrors = None
        TH1 = None
    ROOT = FakeROOT


class FitFunctions:
    @staticmethod
    def poly_1(x, par):
        return par[0] + par[1] * x[0]

    @staticmethod
    def poly_2(x, par):
        return par[0] + par[1] * x[0] + par[2] * x[0] ** 2

    @staticmethod
    def poly_3(x, par):
        return par[0] + par[1] * x[0] + par[2] * x[0] ** 2 + par[3] * x[0] ** 3

    @staticmethod
    def poly_4(x, par):
        return par[0] + par[1] * x[0] + par[2] * x[0] ** 2 + par[3] * x[0] ** 3 + par[4] * x[0] ** 4

    @staticmethod
    def poly_5(x, par):
        return par[0] + par[1] * x[0] + par[2] * x[0] ** 2 + par[3] * x[0] ** 3 + par[4] * x[0] ** 4 + par[5] * x[0] ** 5


FitFunctions.poly_1.__n_par__ = 2
FitFunctions.poly_2.__n_par__ = 3
FitFunctions.poly_3.__n_par__ = 4
FitFunctions.poly_4.__n_par__ = 5
FitFunctions.poly_5.__n_par__ = 6


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
        FitFunctions.poly_1,
        FitFunctions.poly_2,
        FitFunctions.poly_3,
        FitFunctions.poly_4,
        FitFunctions.poly_5,
    ),
    verbose: bool = True,
    root_wrap: bool = True,
) -> Tuple[Tuple[Tuple[Callable, Callable, Callable], ...]]:

    a, b = bounds
    TF1s, Fits = dict(), dict()
    best_chi2_over_ndf, name = float("inf"), ""
    TF1s_up, TF1s_down, Fits_up, Fits_down = dict(), dict(), dict(), dict()

    for _f in function_collection:
        TF1s[_f.__name__] = ROOT.TF1(_f.__name__, _f, a, b, _f.__n_par__)
        Fits[_f.__name__] = graph.Fit(TF1s[_f.__name__], "SFN")

        chi2_over_ndf = Fits[_f.__name__].Chi2() / Fits[_f.__name__].Ndf()
        if abs(chi2_over_ndf - 1) < abs(best_chi2_over_ndf - 1):
            best_chi2_over_ndf = chi2_over_ndf
            name = _f.__name__

            if do_mc_subtr_unc:
                TF1s_up[_f.__name__] = ROOT.TF1(_f.__name__, _f, a, b, _f.__n_par__)
                TF1s_down[_f.__name__] = ROOT.TF1(_f.__name__, _f, a, b, _f.__n_par__)

                Fits_up[_f.__name__] = ff_hist_up.Fit(TF1s_up[_f.__name__], "SFN")
                Fits_down[_f.__name__] = ff_hist_down.Fit(TF1s_down[_f.__name__], "SFN")

    if verbose:
        # TODO add logging
        pass

    _func, (_par, _cov) = getattr(FitFunctions, name), extract_par_and_cov(Fits[name])
    results = [
        (
            get_wrapped_func(_func, _par),
            *get_wrapped_error_funcs(_func, _par, _cov),
        ),
    ]

    if do_mc_subtr_unc:
        _par_up, _cov_up = extract_par_and_cov(Fits_up[name])
        _par_down, _cov_down = extract_par_and_cov(Fits_down[name])

        results.extend(
            [
                (
                    get_wrapped_func(_func, _par_up),
                    *get_wrapped_error_funcs(_func, _par_up, _cov_up),
                ),
                (
                    get_wrapped_func(_func, _par_down),
                    *get_wrapped_error_funcs(_func, _par_down, _cov_down),
                ),
            ]
        )
    if root_wrap:
        results = tuple(
            (
                ROOT.TF1(f"{f.__name__}_{n}", f, a, b, 0),
                ROOT.TF1(f"{f.__name__}_{n}_up", f, a, b, 0),
                ROOT.TF1(f"{f.__name__}_{n}_down", f, a, b, 0),
            )
            for f, n in zip(results, [""] if not do_mc_subtr_unc else ["", "_up", "_down"])
        )
    return tuple(results)
