"""Weighted ``sumw`` / ``sumw2`` accumulators for the b-tag efficiency.

The b-tag efficiency is measured per ``(process, working point, jet flavour)``
combination, binned in jet ``pT`` and jet ``|eta|``. For each such combination
this module stores six per-bin arrays inside a :class:`WeightedEfficiency`
dataclass:

- ``sumw_total`` / ``sumw2_total`` -- sum of event weights (and of squared
  weights) over *all* jets of the flavour;
- ``sumw_pass`` / ``sumw2_pass`` -- the same, restricted to jets passing the
  working-point discriminator cut;
- ``raw_total`` / ``raw_pass`` -- the plain, *unweighted* jet counts (kept for
  the unweighted diagnostic and for the count panels of the plots).

A seventh, scalar field -- ``raw_pt_overflow`` -- carries the real, unweighted
count of jets falling into ROOT's pt-*overflow* bin (i.e. above the top edge
of the finest pt binning): ``Histo2D`` always tracks the overflow bin even
though it is not displayed or included in ``raw_total``, so this is the only
way to see jets above the top pt edge. It is identical for every working
point of a given ``(process, flavour)`` (the total histogram it is read from
does not depend on the working-point cut), and -- unlike the six per-bin
arrays -- it is not touched by bin merging; it is only ever read from the
pre-merge, finest-binning accumulator.

The weighted efficiency (``sumw_pass / sumw_total``) is the production value;
the unweighted, binomial efficiency computed from the raw counts is kept as a
diagnostic. All arrays are plain ``numpy`` (shape ``(n_pt, n_eta)``) so the
dataclass is trivially picklable across the ``ProcessPoolExecutor`` spawn used
by the calculator.

Note on the variance formula
----------------------------
For a per-bin efficiency ``eff = Sp / St`` (``Sp = sumw_pass``,
``St = sumw_total``) the delta method on the two *independent* sums ``Sp`` and
``Sf = St - Sp`` (pass and fail are disjoint fills) gives

    var(eff) = ((1 - eff)^2 * sumw2_pass + eff^2 * sumw2_fail) / sumw_total^2

with ``sumw2_fail = sumw2_total - sumw2_pass``. Because every jet is filled into
``total`` and each *passing* jet is additionally filled into ``pass`` with the
identical weight, ``pass`` is an exact subset of ``total`` and this identity
holds exactly (no rounding slack), so ``sumw2_fail`` is the sum of squared
weights of the failing jets.

All divisions are guarded: empty bins (zero denominator) yield ``nan`` and never
raise here. Range gating (e.g. efficiencies falling outside ``[0, 1]`` because of
negative MC weights, or thresholds on the effective population) is deliberately
*not* applied in this module -- that is handled downstream (Task 18).
"""

from dataclasses import dataclass
from typing import Dict, List, Tuple

import numpy as np


@dataclass
class WeightedEfficiency:
    """Per-(process, WP, flavour) weighted efficiency accumulators.

    Every field is a 2-D ``numpy`` array of shape ``(n_pt, n_eta)`` indexed as
    ``array[i_pt][i_eta]`` (matching the efficiency-grid layout used throughout
    ``btag_efficiency.py``).
    """

    sumw_total: np.ndarray
    sumw2_total: np.ndarray
    sumw_pass: np.ndarray
    sumw2_pass: np.ndarray
    raw_total: np.ndarray
    raw_pass: np.ndarray
    # Scalar: real unweighted count of jets in ROOT's pt-overflow bin (above
    # the top edge of the finest pt binning). Default 0.0 keeps every existing
    # call site (tests included) that does not pass it backward compatible.
    raw_pt_overflow: float = 0.0

    def __post_init__(self) -> None:
        # Coerce anything array-like (lists, tuples, numpy) into float64 arrays
        # so the module functions and pickling behave uniformly.
        self.sumw_total = np.asarray(self.sumw_total, dtype=float)
        self.sumw2_total = np.asarray(self.sumw2_total, dtype=float)
        self.sumw_pass = np.asarray(self.sumw_pass, dtype=float)
        self.sumw2_pass = np.asarray(self.sumw2_pass, dtype=float)
        self.raw_total = np.asarray(self.raw_total, dtype=float)
        self.raw_pass = np.asarray(self.raw_pass, dtype=float)
        self.raw_pt_overflow = float(self.raw_pt_overflow)

        shapes = {
            self.sumw_total.shape,
            self.sumw2_total.shape,
            self.sumw_pass.shape,
            self.sumw2_pass.shape,
            self.raw_total.shape,
            self.raw_pass.shape,
        }
        if len(shapes) != 1:
            raise ValueError(
                f"All WeightedEfficiency arrays must share a shape, got {shapes}."
            )

    @property
    def shape(self) -> Tuple[int, ...]:
        """Return the common ``(n_pt, n_eta)`` shape of the accumulator arrays."""
        return self.sumw_total.shape

    def to_dict(self) -> Dict[str, List]:
        """Return the six raw accumulator arrays plus ``raw_pt_overflow`` as
        nested plain-Python lists (``raw_pt_overflow`` is a plain float)."""
        return {
            "sumw_total": self.sumw_total.tolist(),
            "sumw2_total": self.sumw2_total.tolist(),
            "sumw_pass": self.sumw_pass.tolist(),
            "sumw2_pass": self.sumw2_pass.tolist(),
            "raw_total": self.raw_total.tolist(),
            "raw_pass": self.raw_pass.tolist(),
            "raw_pt_overflow": self.raw_pt_overflow,
        }


def weighted_efficiency(acc: WeightedEfficiency) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(eff, var)`` for the weighted b-tag efficiency.

    ``eff = sumw_pass / sumw_total`` and

        ``var = ((1 - eff)^2 * sumw2_pass + eff^2 * sumw2_fail) / sumw_total^2``

    with ``sumw2_fail = sumw2_total - sumw2_pass``. Bins with a zero
    ``sumw_total`` yield ``nan`` for both ``eff`` and ``var``; nothing is clamped
    here (efficiencies may legitimately fall outside ``[0, 1]`` when negative MC
    weights are involved -- gating is handled downstream).
    """
    sumw_total = acc.sumw_total
    sumw2_fail = acc.sumw2_total - acc.sumw2_pass

    with np.errstate(divide="ignore", invalid="ignore"):
        eff = np.where(sumw_total != 0.0, acc.sumw_pass / sumw_total, np.nan)
        var = np.where(
            sumw_total != 0.0,
            ((1.0 - eff) ** 2 * acc.sumw2_pass + eff**2 * sumw2_fail) / (sumw_total**2),
            np.nan,
        )

    return eff, var


def effective_population(acc: WeightedEfficiency) -> np.ndarray:
    """Return the effective (unweighted-equivalent) population per bin.

    ``n_eff = sumw_total^2 / sumw2_total``; bins with a zero ``sumw2_total``
    (i.e. no fills) yield ``nan``.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(
            acc.sumw2_total != 0.0,
            acc.sumw_total**2 / acc.sumw2_total,
            np.nan,
        )


def unweighted_efficiency(acc: WeightedEfficiency) -> Tuple[np.ndarray, np.ndarray]:
    """Return ``(eff, unc)`` for the unweighted, binomial diagnostic.

    This reproduces exactly the previous, pre-weighting behaviour of
    ``btag_efficiency.py``: from the raw counts,

        ``eff = min(raw_pass / raw_total, 1.0)``
        ``unc = sqrt(eff * (1 - eff) / raw_total)``

    and both are ``0.0`` for empty (``raw_total == 0``) bins. The ``min(.., 1.0)``
    clamp is kept for this *diagnostic* only; the weighted production value is
    never clamped.
    """
    raw_total = acc.raw_total
    raw_pass = acc.raw_pass

    with np.errstate(divide="ignore", invalid="ignore"):
        eff = np.where(raw_total > 0.0, np.minimum(raw_pass / raw_total, 1.0), 0.0)
        unc = np.where(raw_total > 0.0, np.sqrt(eff * (1.0 - eff) / raw_total), 0.0)

    return eff, unc
