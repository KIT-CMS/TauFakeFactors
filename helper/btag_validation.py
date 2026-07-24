"""Strict validation gates, deterministic bin merging, machine-readable report
and atomic payload installation for the b-tag efficiency calculator (Task 18).

The b-tag efficiency payload is a correctionlib ``CorrectionSet`` nested as

    Category(sample_type)
      -> Category(working_point)
           -> Category(jet_flavor)
                -> Binning(jet_eta, flow="error")
                     -> Binning(jet_pt,  flow="clamp")
                          -> float  (weighted efficiency)

This module never produces a *partial* payload: every gate below is fail-fast.
Two very different classes of problem are handled distinctly:

* **Raw pass-set nesting** (``check_wp_nesting``). A jet passing a *tighter*
  working point necessarily passes every *looser* one, so the plain (unweighted)
  pass counts must satisfy ``L >= M >= T >= XT >= XXT`` in every bin. A violation
  can only be a selection / discriminator bug, so it is *immediately fatal*.

* **Weighted quality gates** (``evaluate_gates`` / ``merge_bins``). Because MC
  weights can be negative, the *weighted* efficiency may legitimately fall
  outside ``[0, 1]`` or become non-monotonic across working points, and a bin
  may simply be too sparsely populated (effective denominator too small, or the
  statistical uncertainty too large). None of these are bugs; they are resolved
  deterministically by *merging bins* until every gate passes -- and only a
  violation that survives the coarsest ``1x1`` binning is fatal.

Merging operates on the accumulator *arrays* (summing ``sumw`` / ``sumw2`` / raw
counts), never on efficiencies, and merges all working points of a
``(process, flavor)`` together so the working-point comparability (and the
nesting check) stays meaningful: a bin that fails for *any* working point is
merged for *all* of them.

Thresholds are explicit configuration values (``validation:`` block of the
calculator YAML); the values actually used are recorded verbatim in the
validation report artifact.
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from helper.btag_accumulators import (
    WeightedEfficiency,
    effective_population,
    weighted_efficiency,
    unweighted_efficiency,
)

# The six accumulator arrays that are summed when two bins are merged.
_ACC_FIELDS: Tuple[str, ...] = (
    "sumw_total",
    "sumw2_total",
    "sumw_pass",
    "sumw2_pass",
    "raw_total",
    "raw_pass",
)

# Defaults for the ``validation:`` config block (Task 19 sets the final 2018
# values in config; these are the framework-level fallbacks).
DEFAULT_MIN_EFFECTIVE_DENOMINATOR = 100.0
DEFAULT_MAX_STAT_UNCERTAINTY = 0.25

# Tolerance for the (floating-point) weighted-nesting comparison.
_NEST_TOL = 1e-9

# Merge algorithm version tag stamped into the merge history / report.
MERGE_VERSION = "merge_v1"


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
def get_validation_thresholds(config: Dict) -> Dict[str, float]:
    """Return the ``{min_effective_denominator, max_stat_uncertainty}`` gate.

    Read from the ``validation:`` config block; the framework defaults are used
    for any key that is absent. The returned dict is the value recorded verbatim
    in the validation report (thresholds are never hidden constants).
    """
    block = config.get("validation") or {}
    return {
        "min_effective_denominator": float(
            block.get("min_effective_denominator", DEFAULT_MIN_EFFECTIVE_DENOMINATOR)
        ),
        "max_stat_uncertainty": float(
            block.get("max_stat_uncertainty", DEFAULT_MAX_STAT_UNCERTAINTY)
        ),
    }


# ---------------------------------------------------------------------------
# Input-branch validation (before graph building)
# ---------------------------------------------------------------------------
def validate_input_branches(rdf_columns, required_columns) -> None:
    """Assert that every required column is present before building any graph.

    ``rdf_columns`` may be either

    * a flat iterable of column names (one tree / chain), or
    * a mapping ``{file_path: iterable_of_columns}`` -- the per-file form used by
      the calculator to validate *every* input tree and report which columns are
      missing from which file.

    ``required_columns`` is the list of configured columns (the four probe jet
    columns plus the weight column). Raises :class:`RuntimeError` listing the
    missing names (per file, for the mapping form) if anything is absent.
    """
    required = list(dict.fromkeys(str(c) for c in required_columns))

    if isinstance(rdf_columns, dict):
        missing_by_file: Dict[str, List[str]] = {}
        for file_path, columns in rdf_columns.items():
            present = {str(c) for c in columns}
            missing = [c for c in required if c not in present]
            if missing:
                missing_by_file[str(file_path)] = missing
        if missing_by_file:
            lines = "\n".join(
                f"  {path}: {names}" for path, names in sorted(missing_by_file.items())
            )
            raise RuntimeError(
                f"Missing required input branches in {len(missing_by_file)} "
                f"input file(s):\n{lines}"
            )
        return None

    present = {str(c) for c in rdf_columns}
    missing = [c for c in required if c not in present]
    if missing:
        raise RuntimeError(f"Missing required input branches: {missing}")
    return None


# The throwing per-event check compiled into the RDataFrame graph. Declared once
# per process; the JIT instantiates it for the concrete column types.
_PROBE_DECL_DONE = False
_PROBE_DECL = r"""
namespace btag_validation {
template <typename P, typename E, typename F, typename B, typename W>
bool check_probe_vectors(const ROOT::VecOps::RVec<P> &pt,
                         const ROOT::VecOps::RVec<E> &eta,
                         const ROOT::VecOps::RVec<F> &flav,
                         const ROOT::VecOps::RVec<B> &btag, W weight) {
    const std::size_t n = pt.size();
    if (eta.size() != n || flav.size() != n || btag.size() != n) {
        throw std::runtime_error(
            "btag_validation: probe jet vectors have unequal length in an event");
    }
    if (!std::isfinite(static_cast<double>(weight))) {
        throw std::runtime_error("btag_validation: non-finite event weight");
    }
    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(static_cast<double>(pt[i])) ||
            !std::isfinite(static_cast<double>(eta[i])) ||
            !std::isfinite(static_cast<double>(btag[i]))) {
            throw std::runtime_error("btag_validation: non-finite jet quantity");
        }
    }
    return true;
}
}  // namespace btag_validation
"""


def validate_probe_vectors(rdf, cols: Dict[str, str]):
    """Compile an equal-length + finiteness assertion into the RDF graph.

    ``cols`` maps the roles ``pt``/``eta``/``flavor``/``btag``/``weight`` to the
    configured column names. A ``Filter`` calling a throwing C++ helper is added
    so that the first event with unequal-length probe vectors, a non-finite
    weight or a non-finite jet quantity aborts the event loop (rather than
    silently corrupting the accumulators). Returns the filtered RDataFrame.
    """
    global _PROBE_DECL_DONE
    import ROOT

    if not _PROBE_DECL_DONE:
        ROOT.gInterpreter.Declare(_PROBE_DECL)
        _PROBE_DECL_DONE = True

    expr = (
        f"btag_validation::check_probe_vectors("
        f"{cols['pt']}, {cols['eta']}, {cols['flavor']}, {cols['btag']}, "
        f"{cols['weight']})"
    )
    return rdf.Filter(expr, "btag_probe_vector_validation")


# ---------------------------------------------------------------------------
# Binned accumulator container (one (process, flavor), all working points)
# ---------------------------------------------------------------------------
@dataclass
class BinnedAccumulators:
    """All working-point accumulators for one ``(process, flavor)`` category.

    ``accumulators`` maps working-point name -> :class:`WeightedEfficiency`, in
    loosest-to-tightest order. ``pt_bins`` / ``eta_bins`` are the *shared* bin
    edges (a single binning per category, frozen in the report provenance). All
    working points share the same binning so nesting / comparability checks are
    meaningful.
    """

    accumulators: Dict[str, WeightedEfficiency]
    pt_bins: List[float]
    eta_bins: List[float]

    def copy(self) -> "BinnedAccumulators":
        return BinnedAccumulators(
            accumulators={
                wp: WeightedEfficiency(
                    **{f: np.array(getattr(acc, f), copy=True) for f in _ACC_FIELDS},
                    raw_pt_overflow=acc.raw_pt_overflow,
                )
                for wp, acc in self.accumulators.items()
            },
            pt_bins=list(self.pt_bins),
            eta_bins=list(self.eta_bins),
        )

    @property
    def shape(self) -> Tuple[int, int]:
        return next(iter(self.accumulators.values())).shape


# ---------------------------------------------------------------------------
# Working-point nesting
# ---------------------------------------------------------------------------
@dataclass
class NestingResult:
    """Outcome of :func:`check_wp_nesting`.

    ``raw_ok`` is the *fatal* gate: raw pass counts must nest exactly. The
    weighted fields are informational (they merely trigger merging).
    """

    raw_ok: bool
    raw_violations: List[Dict]
    weighted_ok: bool
    weighted_violations: List[Dict]


def check_wp_nesting(acc_by_wp: Dict[str, WeightedEfficiency]) -> NestingResult:
    """Check working-point pass-set nesting for one ``(process, flavor)``.

    ``acc_by_wp`` must be ordered loosest-to-tightest. For every adjacent
    ``(loose, tight)`` pair and every bin:

    * **raw** ``raw_pass[loose] >= raw_pass[tight]`` (exact subset nesting) --
      a violation is fatal and is reported in ``raw_violations``;
    * **weighted** the weighted efficiency must lie in ``[0, 1]`` and be
      non-increasing from loose to tight -- violations are reported in
      ``weighted_violations`` and merely trigger merging (never fatal).
    """
    wp_names = list(acc_by_wp)
    raw_violations: List[Dict] = []
    weighted_violations: List[Dict] = []

    effs = {wp: weighted_efficiency(acc_by_wp[wp])[0] for wp in wp_names}

    # Range check (per working point).
    for wp in wp_names:
        eff = effs[wp]
        bad = np.argwhere((eff < 0.0) | (eff > 1.0))
        for i_pt, i_eta in bad:
            weighted_violations.append(
                {
                    "kind": "range",
                    "working_point": wp,
                    "i_pt": int(i_pt),
                    "i_eta": int(i_eta),
                    "eff": float(eff[i_pt][i_eta]),
                }
            )

    # Nesting check across adjacent working-point pairs.
    for loose, tight in zip(wp_names, wp_names[1:]):
        raw_loose = acc_by_wp[loose].raw_pass
        raw_tight = acc_by_wp[tight].raw_pass
        bad = np.argwhere(raw_tight > raw_loose)
        for i_pt, i_eta in bad:
            raw_violations.append(
                {
                    "loose": loose,
                    "tight": tight,
                    "i_pt": int(i_pt),
                    "i_eta": int(i_eta),
                    "count_loose": float(raw_loose[i_pt][i_eta]),
                    "count_tight": float(raw_tight[i_pt][i_eta]),
                }
            )

        eff_loose = effs[loose]
        eff_tight = effs[tight]
        bad_w = np.argwhere(eff_tight > eff_loose + _NEST_TOL)
        for i_pt, i_eta in bad_w:
            weighted_violations.append(
                {
                    "kind": "nesting",
                    "loose": loose,
                    "tight": tight,
                    "i_pt": int(i_pt),
                    "i_eta": int(i_eta),
                    "eff_loose": float(eff_loose[i_pt][i_eta]),
                    "eff_tight": float(eff_tight[i_pt][i_eta]),
                }
            )

    return NestingResult(
        raw_ok=not raw_violations,
        raw_violations=raw_violations,
        weighted_ok=not weighted_violations,
        weighted_violations=weighted_violations,
    )


# ---------------------------------------------------------------------------
# Weighted quality gates
# ---------------------------------------------------------------------------
@dataclass
class GateStatus:
    """Per-bin outcome of :func:`evaluate_gates` for one ``BinnedAccumulators``."""

    fail_mask: np.ndarray  # bool, shape (n_pt, n_eta)
    reasons: Dict[Tuple[int, int], List[str]] = field(default_factory=dict)
    failing_working_points: Dict[Tuple[int, int], List[str]] = field(
        default_factory=dict
    )

    @property
    def any_fail(self) -> bool:
        return bool(self.fail_mask.any())

    def failing_pt_rows(self) -> List[int]:
        rows = np.argwhere(self.fail_mask.any(axis=1)).ravel()
        return sorted(int(r) for r in rows)


def evaluate_gates(
    binned: BinnedAccumulators, thresholds: Dict[str, float]
) -> GateStatus:
    """Evaluate every weighted quality gate for all working points of a category.

    A bin ``(i_pt, i_eta)`` *fails* if, for **any** working point, at least one of

    * the effective denominator ``sumw_total^2 / sumw2_total`` is below
      ``min_effective_denominator`` (an empty bin, i.e. ``nan``, always fails);
    * the statistical uncertainty ``sqrt(var)`` exceeds ``max_stat_uncertainty``;
    * the weighted efficiency lies outside ``[0, 1]``;
    * the weighted efficiency breaks working-point nesting (``eff_tight >
      eff_loose``).

    holds. The failing working points are recorded per bin so that
    :func:`merge_bins` (and the report) can attribute the merge.
    """
    min_denom = float(thresholds["min_effective_denominator"])
    max_unc = float(thresholds["max_stat_uncertainty"])

    wp_names = list(binned.accumulators)
    n_pt, n_eta = binned.shape
    fail = np.zeros((n_pt, n_eta), dtype=bool)
    reasons: Dict[Tuple[int, int], List[str]] = {}
    failing_wps: Dict[Tuple[int, int], List[str]] = {}

    def _flag(i_pt: int, i_eta: int, reason: str, wp: str) -> None:
        fail[i_pt][i_eta] = True
        reasons.setdefault((i_pt, i_eta), [])
        if reason not in reasons[(i_pt, i_eta)]:
            reasons[(i_pt, i_eta)].append(reason)
        failing_wps.setdefault((i_pt, i_eta), [])
        if wp not in failing_wps[(i_pt, i_eta)]:
            failing_wps[(i_pt, i_eta)].append(wp)

    effs: Dict[str, np.ndarray] = {}
    for wp in wp_names:
        acc = binned.accumulators[wp]
        eff, var = weighted_efficiency(acc)
        effs[wp] = eff
        unc = np.sqrt(var)
        n_eff = effective_population(acc)

        denom_fail = ~(n_eff >= min_denom)  # nan -> fail
        unc_fail = unc > max_unc  # nan comparisons are False
        range_fail = (eff < 0.0) | (eff > 1.0)

        for i_pt, i_eta in np.argwhere(denom_fail):
            _flag(int(i_pt), int(i_eta), "effective_denominator", wp)
        for i_pt, i_eta in np.argwhere(unc_fail):
            _flag(int(i_pt), int(i_eta), "stat_uncertainty", wp)
        for i_pt, i_eta in np.argwhere(range_fail):
            _flag(int(i_pt), int(i_eta), "eff_out_of_range", wp)

    for loose, tight in zip(wp_names, wp_names[1:]):
        nest_fail = effs[tight] > effs[loose] + _NEST_TOL
        for i_pt, i_eta in np.argwhere(nest_fail):
            _flag(int(i_pt), int(i_eta), "weighted_nesting", tight)

    return GateStatus(
        fail_mask=fail, reasons=reasons, failing_working_points=failing_wps
    )


# ---------------------------------------------------------------------------
# Array merging primitives
# ---------------------------------------------------------------------------
def _merge_pt_rows(binned: BinnedAccumulators, lo: int) -> None:
    """Merge pt rows ``lo`` and ``lo+1`` (in place), summing every array."""
    hi = lo + 1
    for wp, acc in binned.accumulators.items():
        for fname in _ACC_FIELDS:
            arr = getattr(acc, fname)
            merged_row = arr[lo] + arr[hi]
            new_arr = np.delete(arr, hi, axis=0)
            new_arr[lo] = merged_row
            setattr(acc, fname, new_arr)
    binned.pt_bins = binned.pt_bins[:hi] + binned.pt_bins[hi + 1 :]


def _merge_eta_all(binned: BinnedAccumulators) -> None:
    """Collapse all eta columns into one (in place), summing every array."""
    for wp, acc in binned.accumulators.items():
        for fname in _ACC_FIELDS:
            arr = getattr(acc, fname)
            setattr(acc, fname, arr.sum(axis=1, keepdims=True))
    binned.eta_bins = [binned.eta_bins[0], binned.eta_bins[-1]]


def _row_metrics(
    binned: BinnedAccumulators, i_pt: int, thresholds: Dict[str, float]
) -> Dict:
    """Compact per-pt-row gate metrics used for the merge-history pre/post record."""
    min_eff_denom = np.inf
    max_unc = 0.0
    eff_min = np.inf
    eff_max = -np.inf
    for wp, acc in binned.accumulators.items():
        eff, var = weighted_efficiency(acc)
        unc = np.sqrt(var)
        n_eff = effective_population(acc)
        row_denom = n_eff[i_pt]
        row_unc = unc[i_pt]
        row_eff = eff[i_pt]
        finite_denom = row_denom[np.isfinite(row_denom)]
        if finite_denom.size:
            min_eff_denom = min(min_eff_denom, float(np.min(finite_denom)))
        finite_unc = row_unc[np.isfinite(row_unc)]
        if finite_unc.size:
            max_unc = max(max_unc, float(np.max(finite_unc)))
        finite_eff = row_eff[np.isfinite(row_eff)]
        if finite_eff.size:
            eff_min = min(eff_min, float(np.min(finite_eff)))
            eff_max = max(eff_max, float(np.max(finite_eff)))
    return {
        "min_effective_denominator": None if np.isinf(min_eff_denom) else min_eff_denom,
        "max_stat_uncertainty": max_unc,
        "min_weighted_eff": None if np.isinf(eff_min) else eff_min,
        "max_weighted_eff": None if np.isinf(eff_max) else eff_max,
    }


def merge_bins(
    acc: BinnedAccumulators,
    thresholds: Dict[str, float],
    version: str = MERGE_VERSION,
) -> Tuple[BinnedAccumulators, List[Dict]]:
    """Deterministically merge bins until every weighted gate passes.

    While any bin fails (:func:`evaluate_gates`), the lowest-index offending pt
    bin is merged into its lower-pt neighbour (the first pt bin merges upward).
    Only once pt merging is exhausted (a single pt bin remains) are the eta bins
    collapsed. A violation that survives the coarsest ``1x1`` binning is fatal
    (:class:`RuntimeError`). Merging sums the accumulator arrays (never
    efficiencies) and merges all working points of the category together.

    Returns ``(merged, history)`` where ``history`` is the full list of merge
    operations with pre/post gate metrics.
    """
    binned = acc.copy()
    history: List[Dict] = []

    while True:
        status = evaluate_gates(binned, thresholds)
        if not status.any_fail:
            break

        n_pt = len(binned.pt_bins) - 1
        n_eta = len(binned.eta_bins) - 1

        if n_pt > 1:
            p = status.failing_pt_rows()[0]
            lo = 0 if p == 0 else p - 1
            pre = _row_metrics(binned, p, thresholds)
            reasons = sorted(
                {r for k, rs in status.reasons.items() if k[0] == p for r in rs}
            )
            pt_before = list(binned.pt_bins)
            _merge_pt_rows(binned, lo)
            history.append(
                {
                    "op": "merge_pt",
                    "version": version,
                    "offending_pt_index": p,
                    "merged_pt_indices": [lo, lo + 1],
                    "reasons": reasons,
                    "pt_bins_before": pt_before,
                    "pt_bins_after": list(binned.pt_bins),
                    "eta_bins": list(binned.eta_bins),
                    "pre": pre,
                    "post": _row_metrics(binned, lo, thresholds),
                }
            )
        elif n_eta > 1:
            pre = _row_metrics(binned, 0, thresholds)
            reasons = sorted({r for rs in status.reasons.values() for r in rs})
            eta_before = list(binned.eta_bins)
            _merge_eta_all(binned)
            history.append(
                {
                    "op": "merge_eta",
                    "version": version,
                    "reasons": reasons,
                    "eta_bins_before": eta_before,
                    "eta_bins_after": list(binned.eta_bins),
                    "pt_bins": list(binned.pt_bins),
                    "pre": pre,
                    "post": _row_metrics(binned, 0, thresholds),
                }
            )
        else:
            failing = sorted(status.reasons.items())
            raise RuntimeError(
                "b-tag validation gate cannot be satisfied at the coarsest "
                "1x1 (1 pt x 1 eta) binning; unrecoverable failure. "
                f"Failing reasons: {failing}"
            )

    return binned, history


# ---------------------------------------------------------------------------
# SF-flow / lookup-overflow diagnostics
# ---------------------------------------------------------------------------
def raw_flow_counts(binned: BinnedAccumulators) -> Dict[str, float]:
    """Raw jet populations relevant to b-tag SF validity / lookup flow.

    From the (finest, pre-merge) raw ``raw_total`` array and its pt binning:

    * ``pt_lt_30``   -- jets in bins fully below 30 GeV (below SF validity);
    * ``pt_ge_300``  -- jets in bins at or above 300 GeV (SF extrapolation);
    * ``pt_gt_1000`` -- jets at or above 1000 GeV (lookup clamp/overflow): the
      sum of any in-range bins whose lower edge reaches 1000 GeV *plus* the
      real ``raw_pt_overflow`` scalar (ROOT's pt-overflow bin, jets above the
      top edge of the finest binning -- the common case, since every
      configured binning tops out at 1000 GeV).

    All three counts are read from the first working point's accumulator:
    ``raw_total`` and ``raw_pt_overflow`` are identical across working points
    for a given ``(process, flavour)`` (they come from histograms filled from
    the flavour mask only, independent of the working-point cut).
    """
    acc = next(iter(binned.accumulators.values()))
    row_total = acc.raw_total.sum(axis=1)
    pt = binned.pt_bins
    counts = {"pt_lt_30": 0.0, "pt_ge_300": 0.0, "pt_gt_1000": 0.0}
    for i in range(len(pt) - 1):
        lo, hi = pt[i], pt[i + 1]
        if hi <= 30.0:
            counts["pt_lt_30"] += float(row_total[i])
        if lo >= 300.0:
            counts["pt_ge_300"] += float(row_total[i])
        if lo >= 1000.0:
            counts["pt_gt_1000"] += float(row_total[i])
    counts["pt_gt_1000"] += float(acc.raw_pt_overflow)
    return counts


# ---------------------------------------------------------------------------
# JSON helpers
# ---------------------------------------------------------------------------
def _finite_grid(array: np.ndarray) -> List:
    """``array.tolist()`` with non-finite entries -> ``None`` (JSON has no NaN)."""
    return np.where(np.isfinite(array), array, None).tolist()


def _failing_bins(status: GateStatus) -> List[Dict]:
    bins = []
    for (i_pt, i_eta), rs in sorted(status.reasons.items()):
        bins.append(
            {
                "i_pt": i_pt,
                "i_eta": i_eta,
                "reasons": rs,
                "working_points": status.failing_working_points.get((i_pt, i_eta), []),
            }
        )
    return bins


def summarize_category(
    pre: BinnedAccumulators,
    post: BinnedAccumulators,
    history: List[Dict],
    thresholds: Dict[str, float],
    nesting: NestingResult,
) -> Dict:
    """Assemble the machine-readable per-category report entry.

    Records pre/post binning, pre/post gate results, the full merge history,
    per-working-point effective populations and weighted-vs-unweighted deltas
    (post-merge), the SF-flow / overflow raw populations (pre-merge, finest
    binning) and the nesting outcome.
    """
    pre_status = evaluate_gates(pre, thresholds)
    post_status = evaluate_gates(post, thresholds)

    eff_pop: Dict[str, List] = {}
    delta: Dict[str, List] = {}
    for wp, acc in post.accumulators.items():
        eff_pop[wp] = _finite_grid(effective_population(acc))
        w_eff, _ = weighted_efficiency(acc)
        u_eff, _ = unweighted_efficiency(acc)
        delta[wp] = _finite_grid(w_eff - u_eff)

    return {
        "pt_bins_pre": list(pre.pt_bins),
        "eta_bins_pre": list(pre.eta_bins),
        "pt_bins_post": list(post.pt_bins),
        "eta_bins_post": list(post.eta_bins),
        "gates_pre": {
            "n_failing_bins": int(pre_status.fail_mask.sum()),
            "failing_bins": _failing_bins(pre_status),
        },
        "gates_post": {
            "n_failing_bins": int(post_status.fail_mask.sum()),
            "failing_bins": _failing_bins(post_status),
        },
        "merge_history": history,
        "effective_population_post": eff_pop,
        "weighted_minus_unweighted_eff_post": delta,
        "raw_flow_counts": raw_flow_counts(pre),
        "nesting": {
            "raw_ok": nesting.raw_ok,
            "weighted_ok": nesting.weighted_ok,
            "n_raw_violations": len(nesting.raw_violations),
            "n_weighted_violations": len(nesting.weighted_violations),
        },
    }


# ---------------------------------------------------------------------------
# correctionlib round-trip gate
# ---------------------------------------------------------------------------
def roundtrip_check(
    json_path: str,
    all_efficiencies: Dict[str, Dict[str, Dict[str, List[List[float]]]]],
    all_flavor_bins: Dict[str, Dict[str, Tuple[List[float], List[float]]]],
    wps: Dict[str, float],
    flavors: Dict[str, int],
) -> Tuple[bool, int, List[Dict]]:
    """Evaluate the produced JSON via correctionlib and compare to stored floats.

    Every ``(sample_type, working_point, flavor)`` leaf is probed at each bin
    centre and at each (in-range) lower bin edge in both pt and eta, and the
    value returned by ``correctionlib`` is required to match the efficiency
    stored in ``all_efficiencies`` to within floating-point serialization
    tolerance (``math.isclose`` with ``rel_tol=1e-9``, ``abs_tol=1e-12``).

    An exact ``==`` comparison is wrong here: correctionlib serializes the
    efficiency to JSON text and re-parses it on lookup, so a value can come back
    differing in its last unit-in-the-last-place (~1e-16 for a float64 in
    ``[0, 1]``). The tolerance is far tighter than any physically meaningful
    efficiency difference, so genuine build errors (wrong bin, transposed axis,
    stale value) still fail the gate, while pure round-off does not. Returns
    ``(ok, n_checked, mismatches)``.
    """
    import correctionlib

    cset = correctionlib.CorrectionSet.from_file(json_path)
    corr = cset["btag_efficiency"]

    mismatches: List[Dict] = []
    n_checked = 0

    for sample_type, proc_eff in all_efficiencies.items():
        for wp in wps:
            for flavor_name, flavor_id in flavors.items():
                pt_bins, eta_bins = all_flavor_bins[sample_type][flavor_name]
                eff_grid = proc_eff[wp][flavor_name]
                n_pt = len(pt_bins) - 1
                n_eta = len(eta_bins) - 1
                for i_eta in range(n_eta):
                    eta_probes = [
                        (eta_bins[i_eta] + eta_bins[i_eta + 1]) / 2.0,  # centre
                        eta_bins[i_eta],  # left (in-range) edge
                    ]
                    for i_pt in range(n_pt):
                        pt_probes = [
                            (pt_bins[i_pt] + pt_bins[i_pt + 1]) / 2.0,
                            pt_bins[i_pt],
                        ]
                        expected = eff_grid[i_pt][i_eta]
                        for eta_probe in eta_probes:
                            for pt_probe in pt_probes:
                                got = corr.evaluate(
                                    sample_type,
                                    wp,
                                    int(flavor_id),
                                    float(eta_probe),
                                    float(pt_probe),
                                )
                                n_checked += 1
                                if expected is None or got is None:
                                    is_match = (expected is None and got is None)
                                else:
                                    is_match = math.isclose(
                                        got, expected, rel_tol=1e-9, abs_tol=1e-12
                                    )
                                if not is_match:
                                    mismatches.append(
                                        {
                                            "sample_type": sample_type,
                                            "working_point": wp,
                                            "flavor": flavor_name,
                                            "i_pt": i_pt,
                                            "i_eta": i_eta,
                                            "eta": float(eta_probe),
                                            "pt": float(pt_probe),
                                            "expected": expected,
                                            "got": got,
                                        }
                                    )
    return (not mismatches), n_checked, mismatches


# ---------------------------------------------------------------------------
# Validation report
# ---------------------------------------------------------------------------
def write_validation_report(
    path: str,
    *,
    status: str,
    thresholds: Dict[str, float],
    merge_version: str,
    categories: Dict,
    roundtrip: Dict,
    era: str = "",
    channel: str = "",
    channels_present: Optional[List[str]] = None,
    required_channels: Optional[List[str]] = None,
    extra: Optional[Dict] = None,
) -> Dict:
    """Assemble and write the machine-readable validation report JSON.

    Returns the assembled report dict. ``allow_nan=False`` guarantees the file is
    strictly valid JSON (all non-finite values must already be ``None``).
    """
    report: Dict = {
        "status": status,
        "thresholds": dict(thresholds),
        "merge_version": merge_version,
        "era": era,
        "channel": channel,
        "channels_present": list(channels_present or []),
        "required_channels": list(required_channels or []),
        "categories": categories,
        "roundtrip": roundtrip,
    }
    if extra:
        report.update(extra)

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fout:
        json.dump(report, fout, indent=4, allow_nan=False)
    return report


# ---------------------------------------------------------------------------
# Atomic payload installation
# ---------------------------------------------------------------------------
def _sha256(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fin:
        for chunk in iter(lambda: fin.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _atomic_replace_dir(tmp_dir: str, target_dir: str) -> None:
    """Atomically move ``tmp_dir`` onto ``target_dir`` (restoring on failure)."""
    parent = os.path.dirname(os.path.abspath(target_dir))
    os.makedirs(parent, exist_ok=True)
    if os.path.exists(target_dir):
        backup = f"{target_dir}.bak-{os.getpid()}"
        shutil.rmtree(backup, ignore_errors=True)
        os.replace(target_dir, backup)
        try:
            os.replace(tmp_dir, target_dir)
        except Exception:
            os.replace(backup, target_dir)  # restore the previous payload
            raise
        shutil.rmtree(backup, ignore_errors=True)
    else:
        os.replace(tmp_dir, target_dir)


def stage_payload(payload_files: List[str], provenance: Dict, target_dir: str) -> str:
    """Stage and checksum-verify a payload into a fresh ``target_dir.tmp-<pid>``.

    This is the first phase of a two-phase install: it performs the same fast
    refusals and checksum verification as :func:`atomic_install` but stops
    short of touching ``target_dir`` itself. On success, returns the staged
    tmp directory's path; the caller must eventually either commit it
    (:func:`commit_staged`) or discard it (``shutil.rmtree``). On any failure
    the tmp directory is cleaned up before the exception propagates, so a
    failed call never leaves a stray ``.tmp-<pid>`` directory behind.

    Refusals (raised as :class:`RuntimeError`, before anything is staged):

    * ``provenance["status"] != "passed"``;
    * some entry of ``provenance["required_channels"]`` is missing from
      ``provenance["channels_present"]``.

    This split enables a multi-channel, all-or-nothing publish: every
    channel's payload can be staged and verified before *any* channel's
    ``target_dir`` is touched (see ``btag_efficiency.py:publish_all_channels``).
    """
    status = provenance.get("status")
    present = set(provenance.get("channels_present", []))
    required = set(provenance.get("required_channels", []))
    checksums = provenance.get("checksums", {})

    # Fast refusals (before staging anything).
    if status != "passed":
        raise RuntimeError(
            f"stage_payload refuses to stage: report status is '{status}', "
            f"not 'passed'."
        )
    missing_channels = sorted(required - present)
    if missing_channels:
        raise RuntimeError(
            f"stage_payload refuses to stage: required channels "
            f"{missing_channels} not present (have {sorted(present)})."
        )

    tmp_dir = f"{target_dir}.tmp-{os.getpid()}"
    shutil.rmtree(tmp_dir, ignore_errors=True)
    os.makedirs(tmp_dir)
    try:
        for src in payload_files:
            base = os.path.basename(src)
            shutil.copy2(src, os.path.join(tmp_dir, base))
        # Verify checksums of the staged copies (detects truncation/corruption).
        for src in payload_files:
            base = os.path.basename(src)
            expected = checksums.get(base)
            if expected is None:
                raise RuntimeError(
                    f"stage_payload: no checksum recorded for payload file '{base}'."
                )
            actual = _sha256(os.path.join(tmp_dir, base))
            if actual != expected:
                raise RuntimeError(
                    f"stage_payload: checksum mismatch for '{base}' "
                    f"(expected {expected}, got {actual}); aborting install."
                )
    except Exception:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise

    return tmp_dir


def commit_staged(tmp_dir: str, target_dir: str) -> str:
    """Commit an already-staged, checksum-verified tmp dir onto ``target_dir``.

    Second phase of the two-phase install (see :func:`stage_payload`): performs
    the ``os.replace`` swap via :func:`_atomic_replace_dir`. On failure the tmp
    directory is cleaned up before the exception propagates. Returns
    ``target_dir`` on success.
    """
    try:
        _atomic_replace_dir(tmp_dir, target_dir)
    except Exception:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise
    return target_dir


def atomic_install(payload_files: List[str], provenance: Dict, target_dir: str) -> str:
    """Install the payload into ``target_dir`` all-or-nothing (single channel).

    Composes :func:`stage_payload` (staging + checksum verification) followed
    by :func:`commit_staged` (the ``os.replace`` swap). The staging directory
    is ``os.replace``-d into place **only if**

    * ``provenance["status"] == "passed"``, and
    * every ``provenance["required_channels"]`` is in
      ``provenance["channels_present"]``.

    Any refusal or checksum error cleans up the staging directory and leaves an
    existing ``target_dir`` untouched. Returns ``target_dir`` on success.
    """
    tmp_dir = stage_payload(payload_files, provenance, target_dir)
    return commit_staged(tmp_dir, target_dir)
