import ROOT
import numpy as np
import warnings
from typing import List, Tuple, Union, Any

_EXTRA_PARAM_MEANS = "_extra_weighted_means"
_EXTRA_PARAM_COUNTS = "_extra_weighted_counts"
_EXTRA_PARAM_FLAG = "_has_extra_params"


class PatchedRDataFrame:
    def __init__(self, base_rdf: ROOT.RDataFrame):
        self.base_rdf = base_rdf

    def __getattr__(self, attr: str) -> Any:
        return getattr(self.base_rdf, attr)  # any attribute not overridden shall be used

    def Histo1D(self, histo_params, x_expr, weight_expr):
        hist = self.base_rdf.Histo1D(histo_params, x_expr, weight_expr)
        _GetValue = hist.GetValue

        def new_GetValue():
            main_hist = _GetValue()
            if isinstance(histo_params, (list, tuple)):
                if len(histo_params) == 5:
                    _, _, nbins, low, high = histo_params  # dont care about name and title
                    bins = np.linspace(low, high, int(nbins) + 1)
                elif len(histo_params) >= 4:
                    _, _, nbins, bins = histo_params  # dont care about name and title
                else:
                    raise ValueError("Wrapper only supports 4 or 5 parameters")

                arr = self.base_rdf.AsNumpy([x_expr, weight_expr])
                xvals, weights = arr[x_expr], arr[weight_expr]

                counts, _ = np.histogram(xvals, bins=bins, weights=weights)
                sums, _ = np.histogram(xvals, bins=bins, weights=xvals * weights)
                weighted_means = [s / c if c != 0 else 0 for s, c in zip(sums, counts)]
                flag = True
            else:
                counts, sums, flag = None, None, False

            setattr(main_hist, _EXTRA_PARAM_COUNTS, counts)
            setattr(main_hist, _EXTRA_PARAM_MEANS, weighted_means)
            setattr(main_hist, _EXTRA_PARAM_FLAG, flag)

            # Overwrite GetValue so further calls return the histogram directly.
            main_hist.GetValue = lambda: main_hist
            return main_hist

        hist.GetValue = new_GetValue
        return hist


def calc_center_of_mass(
    means: List[List[float]],
    weights: List[List[float]],
    factor: float = 1.0,
) -> Tuple[List[float], List[float]]:
    """
    Calculate the center of mass for the given means and weights.

    Args:
        means (List[List[float]]): List of means for each bin.
        weights (List[List[float]]): List of weights for each bin.
        factor (float): Scaling factor for the weights.

    Returns:
        Tuple[List[float], List[float]]: Tuple containing the updated counts and means.
    """
    assert len(means) == len(weights) == 2, "Means and weights must be lists of two elements each."

    _counts, _means = [], []
    for m1, m2, w1, w2 in zip(*means, *weights):
        tot = w1 + factor * w2
        _counts.append(tot)
        if tot != 0:
            _means.append((w1 * m1 + factor * w2 * m2) / tot)
        else:
            _means.append(0)
    return _counts, _means


# unpatched
_original_Add = object.__getattribute__(ROOT.TH1, "Add")
_original_Multiply = object.__getattribute__(ROOT.TH1, "Multiply")
_original_Divide = object.__getattribute__(ROOT.TH1, "Divide")
_original_Clone = object.__getattribute__(ROOT.TH1, "Clone")


def patched_Add(
    self: ROOT.TH1,
    other: ROOT.TH1,
    factor: float = 1.0,
) -> ROOT.TH1:
    """
    Patched version of TH1.Add that calls the original Add (using the stored _original_Add)
    and updates extra attributes _extra_weighted_counts and _extra_weighted_means from
    patched Histo1D of RDataFrame.

    Bin‐by‐bin values are computed as:
       new_count = count_self + factor * count_other
       new_mean = ( count_self * mean_self + factor * count_other * mean_other ) / new_count   (if new_count != 0)

    Args:
        self (ROOT.TH1): The histogram to which the other histogram is added.
        other (ROOT.TH1): The histogram to be added.
        factor (float): The scaling factor for the other histogram.
    """
    _original_Add(self, other, factor)

    if (hasattr(self, _EXTRA_PARAM_FLAG)):
        _counts, _means = calc_center_of_mass(
            (getattr(self, _EXTRA_PARAM_MEANS), getattr(other, _EXTRA_PARAM_MEANS)),
            (getattr(self, _EXTRA_PARAM_COUNTS), getattr(other, _EXTRA_PARAM_COUNTS)),
            factor,
        )
        setattr(self, _EXTRA_PARAM_COUNTS, _counts)
        setattr(self, _EXTRA_PARAM_MEANS, _means)
        setattr(self, _EXTRA_PARAM_FLAG, True)

    return self


def patched_Multiply(
    self: ROOT.TH1,
    other: ROOT.TH1,
    factor: float = 1.0,
) -> ROOT.TH1:
    """
    Patched version of TH1.Multiply.
    Calls the original Multiply (which accepts only the other histogram),
    and then updates the extra attributes (_extra_weighted_counts and _extra_weighted_means)
    using calc_center_of_mass and our supplied factor.
    """
    # Call original Multiply without factor (because the overload expects a single argument).
    _original_Multiply(self, other)

    if (hasattr(self, _EXTRA_PARAM_FLAG)):
        _counts, _means = calc_center_of_mass(
            (getattr(self, _EXTRA_PARAM_MEANS), getattr(other, _EXTRA_PARAM_MEANS)),
            (getattr(self, _EXTRA_PARAM_COUNTS), getattr(other, _EXTRA_PARAM_COUNTS)),
            factor,
        )
        setattr(self, _EXTRA_PARAM_COUNTS, _counts)
        setattr(self, _EXTRA_PARAM_MEANS, _means)
        setattr(self, _EXTRA_PARAM_FLAG, True)

    return self


def patched_Divide(
    self: ROOT.TH1,
    other: ROOT.TH1,
    factor: float = 1.0,
) -> ROOT.TH1:
    """
    Patched version of TH1.Divide.
    Calls the original Divide (which accepts only the other histogram),
    and then updates extra attributes using calc_center_of_mass with a scaling factor.
    """
    # Call original Divide without factor.
    _original_Divide(self, other)

    if (hasattr(self, _EXTRA_PARAM_FLAG)):
        _counts, _means = calc_center_of_mass(
            (getattr(self, _EXTRA_PARAM_MEANS), getattr(other, _EXTRA_PARAM_MEANS)),
            (getattr(self, _EXTRA_PARAM_COUNTS), getattr(other, _EXTRA_PARAM_COUNTS)),
            factor,
        )
        setattr(self, _EXTRA_PARAM_COUNTS, _counts)
        setattr(self, _EXTRA_PARAM_MEANS, _means)
        setattr(self, _EXTRA_PARAM_FLAG, True)

    return self


def patched_Clone(
    self: ROOT.TH1,
    name: Union[str, None] = None,
) -> ROOT.TH1:
    """
    Patched Clone: copy the extra attributes into the clone, if they exist.

    Args:
        self (ROOT.TH1): The histogram to clone.
        name (str): The name of the new histogram.
    """
    if name is not None:
        clone = _original_Clone(self, name)
    else:
        clone = _original_Clone(self)

    if hasattr(self, _EXTRA_PARAM_FLAG):
        setattr(clone, _EXTRA_PARAM_FLAG, getattr(self, _EXTRA_PARAM_FLAG))
    if hasattr(self, _EXTRA_PARAM_COUNTS):
        setattr(clone, _EXTRA_PARAM_COUNTS, getattr(self, _EXTRA_PARAM_COUNTS))
    if hasattr(self, _EXTRA_PARAM_MEANS):
        setattr(clone, _EXTRA_PARAM_MEANS, getattr(self, _EXTRA_PARAM_MEANS))

    return clone


if not hasattr(ROOT.TH1, "_patched"):
    ROOT.TH1.Add = patched_Add
    ROOT.TH1.Multiply = patched_Multiply
    ROOT.TH1.Divide = patched_Divide
    ROOT.TH1.Clone = patched_Clone
    # Mark the class as patched to avoid multiple patching
    ROOT.TH1._patched = True
