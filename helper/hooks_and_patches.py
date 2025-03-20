from itertools import islice, tee
from typing import Any, Callable, List, Tuple, Union

import numpy as np
import ROOT

_EXTRA_PARAM_MEANS = "_extra_weighted_means"
_EXTRA_PARAM_COUNTS = "_extra_weighted_counts"
_EXTRA_PARAM_FLAG = "_has_extra_params"


class PassThroughWrapper:
    """
    A transparent pass-through wrapper that delegates attribute access,
    assignment, and deletion to the wrapped object.

    This wrapper does not alter the underlying object's methods or attributes.
    If inherited from, any overridden methods or attributes in the subclass
    will be triggered first.
    """

    def __init__(self, obj: Any) -> None:
        """
        Initialize with given object.

        Args:
            obj: Object to wrap.
        """
        object.__setattr__(self, "_obj", obj)

    def __getattribute__(self, name: str) -> Any:
        """
        Retrieve attribute 'name'. If attribute is defined on this wrapper (or its subclass),
        that attribute is returned then. Otherwise, it is delegated to the wrapped object.

        Args:
            name: Name of the attribute to retrieve.

        Returns:
            Attribute value.
        """
        try: # Try to get the attribute from this instance or subclass.
            return object.__getattribute__(self, name)
        except AttributeError:  # Fallback: delegate to the wrapped object.
            _obj = object.__getattribute__(self, "_obj")
            return getattr(_obj, name)

    def __setattr__(self, name: str, value: Any) -> None:
        """
        Delegate setting of attribute 'name' to the wrapped object.

        Args:
            name: Name of the attribute to set.
            value: The value to set.
        """
        _obj = object.__getattribute__(self, "_obj")
        setattr(_obj, name, value)

    def __delattr__(self, name: str) -> None:
        """
        Delegate deletion of attribute 'name' to the wrapped object.

        Args:
            name: Name of the attribute to delete.
        """
        _obj = object.__getattribute__(self, "_obj")
        delattr(_obj, name)


class Histo1DPatchedRDataFrame(PassThroughWrapper):
    """
    A wrapper around ROOT.RDataFrame patching Histo1D method adding extra attributes 
    (weighted counts and weighted means) to the produced histograms.
    """

    def __init__(self, base_rdf: ROOT.RDataFrame):
        """
        Initialization.

        Args:
            base_rdf (ROOT.RDataFrame): The original RDataFrame instance.
        """
        super().__init__(base_rdf)
        self.base_rdf: ROOT.RDataFrame = base_rdf

    def Histo1D(
        self,
        histo_params: Union[List[Any], Tuple[Any, ...]],
        x_expr: str,
        weight_expr: str,
    ) -> Any:
        """
        Create a 1D histogram with extra attributes injecting additional processing
        information into the histogram, replacing the original Histo1D method.

        Args:
            histo_params (Union[List, Tuple]): Parameters for the histogram (expects 4 or 5 elements).
            x_expr (str): Expression to select the x values.
            weight_expr (str): Expression to select the weights.

        Returns:
            ROOT.TH1: Patched histogram with extra attributes.
        """
        hist = self.base_rdf.Histo1D(histo_params, x_expr, weight_expr)
        _GetValue = hist.GetValue

        def new_GetValue():
            main_hist = _GetValue()
            if isinstance(histo_params, (list, tuple)):
                if len(histo_params) == 5:
                    _, _, nbins, low, high = histo_params  # dont care about name and title
                    bins = np.linspace(low, high, int(nbins) + 1)
                elif len(histo_params) == 4:
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
    bin_edges: List[float],
    factor: float = 1.0,
    weights_combination_operation: Callable = lambda a, b: a + b,
) -> Tuple[List[float], List[float]]:
    """
    Calculate the center of mass for the given means and weights.

    Args:
        means (List[List[float]]): List of means for each bin.
        weights (List[List[float]]): List of weights for each bin.
        bin_edges (List[float]): List of bin edges.
        factor (float): Scaling factor for the weights.
        weights_combination_operation (Callable): Function to combine weights, default is addition.

    Returns:
        Tuple[List[float], List[float]]: Tuple containing the updated counts and means.
    """
    assert len(means) == len(weights) == 2, "Means and weights must be lists of two elements each."

    def nwise(g, *, n=2):
        return zip(*(islice(g, i, None) for i, g in enumerate(tee(g, n))))

    _counts, _means = [], []
    for m1, m2, w1, w2, (bin_low, bin_high) in zip(*means, *weights, nwise(bin_edges, n=2)):
        _counts.append(weights_combination_operation(w1, w2))
        if (tot := w1 + factor * w2) != 0:
            _means.append((w1 * m1 + factor * w2 * m2) / tot)
        else:
            _means.append((bin_low + bin_high) / 2)  # Default to bin center if no counts
    return _counts, _means


# unpatched
_original_Add = object.__getattribute__(ROOT.TH1, "Add")
_original_Multiply = object.__getattribute__(ROOT.TH1, "Multiply")
_original_Divide = object.__getattribute__(ROOT.TH1, "Divide")
_original_Clone = object.__getattribute__(ROOT.TH1, "Clone")
_original_Scale = object.__getattribute__(ROOT.TH1, "Scale")


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

    Returns:
        ROOT.TH1: The histogram after addition.
    """
    _original_Add(self, other, factor)

    if (hasattr(self, _EXTRA_PARAM_FLAG)):
        _counts, _means = calc_center_of_mass(
            means=(getattr(self, _EXTRA_PARAM_MEANS), getattr(other, _EXTRA_PARAM_MEANS)),
            weights=(getattr(self, _EXTRA_PARAM_COUNTS), getattr(other, _EXTRA_PARAM_COUNTS)),
            bin_edges=[self.GetBinLowEdge(i) for i in range(1, self.GetNbinsX() + 2)],
            factor=factor,
            weights_combination_operation=lambda a, b: a + factor * b,
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

    Args:
        self (ROOT.TH1): The histogram to be multiplied.
        other (ROOT.TH1): The histogram by which to multiply.
        factor (float): The scaling factor for the other histogram.

    Returns:
        ROOT.TH1: The histogram after multiplication.
    """
    # Call original Multiply without factor (because the overload expects a single argument).
    _original_Multiply(self, other)

    if (hasattr(self, _EXTRA_PARAM_FLAG)):
        _counts, _means = calc_center_of_mass(
            means=(getattr(self, _EXTRA_PARAM_MEANS), getattr(other, _EXTRA_PARAM_MEANS)),
            weights=(getattr(self, _EXTRA_PARAM_COUNTS), getattr(other, _EXTRA_PARAM_COUNTS)),
            bin_edges=[self.GetBinLowEdge(i) for i in range(1, self.GetNbinsX() + 2)],
            factor=factor,
            weights_combination_operation=lambda a, b: a * factor * b,
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

    Args:
        self (ROOT.TH1): The histogram to be divided.
        other (ROOT.TH1): The histogram by which to divide.
        factor (float): The scaling factor for the other histogram.

    Returns:
        ROOT.TH1: The histogram after division.
    """
    # Call original Divide without factor.
    _original_Divide(self, other)

    if (hasattr(self, _EXTRA_PARAM_FLAG)):
        _counts, _means = calc_center_of_mass(
            means=(getattr(self, _EXTRA_PARAM_MEANS), getattr(other, _EXTRA_PARAM_MEANS)),
            weights=(getattr(self, _EXTRA_PARAM_COUNTS), getattr(other, _EXTRA_PARAM_COUNTS)),
            bin_edges=[self.GetBinLowEdge(i) for i in range(1, self.GetNbinsX() + 2)],
            factor=-factor,
            weights_combination_operation=lambda a, b: a / (factor * b) if b != 0 else 0,
        )
        setattr(self, _EXTRA_PARAM_COUNTS, _counts)
        setattr(self, _EXTRA_PARAM_MEANS, _means)
        setattr(self, _EXTRA_PARAM_FLAG, True)

    return self


def patched_Scale(
    self: ROOT.TH1,
    factor: float,
) -> ROOT.TH1:
    """
    Patched version of TH1.Scale, accepting a scaling factor only.
    Calls the original Scale and then updates extra attributes using calc_center_of_mass.

    Args:
        self (ROOT.TH1): The histogram to scale.
        factor (float): The scaling factor.

    Returns:
        ROOT.TH1: The scaled histogram.
    """
    _original_Scale(self, factor)

    if (hasattr(self, _EXTRA_PARAM_FLAG)):
        _counts, _means = calc_center_of_mass(
            means=(getattr(self, _EXTRA_PARAM_MEANS), getattr(self, _EXTRA_PARAM_COUNTS)),
            weights=(getattr(self, _EXTRA_PARAM_COUNTS), getattr(self, _EXTRA_PARAM_COUNTS)),
            bin_edges=[self.GetBinLowEdge(i) for i in range(1, self.GetNbinsX() + 2)],
            factor=factor,
            weights_combination_operation=lambda a, b: a * factor * b,
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

    Returns:
        ROOT.TH1: The cloned histogram.
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
