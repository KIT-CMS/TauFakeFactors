#!/usr/bin/env python3
"""Check a TauFakeFactors b-tag efficiency payload against the strict consumer.

The consumer is bbtautau CROWN's
cpp_addons/src/btag_sf_strict.cxx::multi_wp_event_weight, used for all four
Run-2 eras (2016preVFP, 2016postVFP, 2017, 2018) with the working points the
analysis profile configures ([T, M] since v7). It throws unless every
efficiency it evaluates is finite and in (0, 1], the efficiencies decrease
strictly from the loosest to the tightest working point, and 1 - eff(loosest)
is positive. This script walks the correctionlib JSON and reports every
violation.

The working points are read from the payload (its working_point keys); with
--era or --sf-check they are ordered by the UParTAK4_wp_values thresholds of
the BTV payload, otherwise by their names (L < M < T < XT < XXT). CROWN only
evaluates the working points of the analysis profile
(analysis_profiles.btag_sf_working_points = T, M), so with --era (or --wps)
the payload's working-point set must equal that list (a payload with a
leftover L would otherwise be emulated with a "fails L" region the consumer
never builds), and --sf-check emulates exactly that list.

For each file it prints, next to the verdict, an "inspected:" line naming what
was actually walked (sample types, working points, flavours, the eta x pt grid,
the resulting cells and the number of adjacent-working-point comparisons made),
so that a payload short a sample type or a working point is self-evident. When
the payload's grid is ragged, that line says so instead of quoting one shape.
With --era the era's SM sample types are required (a missing one is a
violation, the consumer would fail on the first job of that process); without
it, structural completeness is only reported. Violations are tallied by class
before the truncated list, and bins are named in GeV and |eta| rather than by
index, because the remedy is coarser per-process jet_pt_bins / jet_eta_bins.

--sf-check additionally recomputes the consumer's jet weight with the BTV scale
factors of the era (pinned 2026-06-18 payload, or --sf-payload), including the
consumer's treatment of b/c jets below the calibrated pT range (SF of the first
bin with twice the deviation) and the BTV rule for a negative numerator of an
intermediate factor (set to 1, counted, not a violation). It also prints, per
group of sample types with identical values (one pool), the smallest jet weight
of a jet failing the loosest working point next to its efficiency and the SF
ceiling 1/SF there, and warns (no violation) below --warn-below.

Usage (from the TauFakeFactors directory; run_eff_<era>.sh calls it after the measurement):
    python3 scripts/check_btag_payload.py --era 2017 --sf-check \\
        workdir/btag_eff_2017/2017/*/btag_efficiency.json.gz
"""
import argparse
import gzip
import json
import math
import sys

# Canonical name order, loosest to tightest; only used when no BTV payload is
# at hand to order the working points by their score thresholds.
WP_NAME_ORDER = ["L", "M", "T", "XT", "XXT"]

# Violation classes, in the order they are tallied.
KINDS = ["NOT FINITE", "eff <= 0", "eff > 1", "NOT MONOTONIC", "DEGENERATE",
         "JET WEIGHT <= 0", "SF EVAL FAILED", "MISSING SAMPLE TYPE",
         "UNKNOWN WP", "WP SET"]

RUN2_ERAS = ["2016preVFP", "2016postVFP", "2017", "2018"]

# The working points the consumer evaluates in all four Run-2 eras: CROWN
# bbtautau analysis_profiles.btag_sf_working_points (v7 E14, C14).
CONSUMER_WPS = ("T", "M")

# Pinned BTV UParTAK4 payloads, as CROWN bbtautau btag_payloads.py:22-28.
BTV_PIN = (
    "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/"
    "Run2-{era}-UL-NanoAODv15/2026-06-18/btagging.json.gz"
)
PINNED_BTV = BTV_PIN.format(era="2018")  # default of --sf-payload without --era

# SM sample types CROWN looks up in the payload of each era (sm_config /
# sm_btag_efficiency_config AVAILABLE_SAMPLES; v7 adds rem_htautau). 2016postVFP
# has no rem_htautau sample, but its payload carries the key as well (written
# from the diboson pool), so that it covers every type the config accepts.
SM_SAMPLE_TYPES = [
    "hh2b2tau", "dyjets", "wjets", "ttbar", "singletop", "diboson",
    "electroweak_boson", "ggh_htautau", "vbf_htautau", "vbf_hbb",
    "rem_hbb", "rem_higgs", "rem_htautau", "rem_ttbar",
]
REQUIRED_SAMPLE_TYPES = {era: SM_SAMPLE_TYPES for era in RUN2_ERAS}
# --require-sm-2018 (deprecated): the 13 types of the 2018 v1-v6 payloads.
SM_2018_SAMPLE_TYPES = [name for name in SM_SAMPLE_TYPES if name != "rem_htautau"]


def name_sort_key(wp):
    """Order known working points loosest to tightest by name, unknown ones last."""
    known = WP_NAME_ORDER.index(wp) if wp in WP_NAME_ORDER else len(WP_NAME_ORDER)
    return (known, str(wp))


def order_wps(wps, thresholds=None):
    """Working points loosest to tightest: by BTV threshold when known, else by name."""
    if thresholds:
        return sorted(wps, key=lambda wp: (thresholds.get(wp, math.inf), name_sort_key(wp)))
    return sorted(wps, key=name_sort_key)


def edge_str(value):
    """Render a bin edge without a trailing '.0' (20.0 -> '20', 1.4 -> '1.4')."""
    try:
        return f"{float(value):g}"
    except (TypeError, ValueError):
        return str(value)


def bin_label(edges, index, name):
    """Name one bin in physical units, e.g. 'pt=[140,200)'.

    Falls back to the bin index when the node carries no usable edge list, so a
    hand-made or uniform-binning payload is reported rather than crashing.
    """
    if isinstance(edges, dict):  # correctionlib uniform binning
        try:
            width = (float(edges["high"]) - float(edges["low"])) / int(edges["n"])
            low = float(edges["low"]) + index * width
            return f"{name}=[{edge_str(low)},{edge_str(low + width)})"
        except (KeyError, TypeError, ValueError, ZeroDivisionError):
            return f"{name}_bin={index}"
    try:
        return f"{name}=[{edge_str(edges[index])},{edge_str(edges[index + 1])})"
    except (IndexError, KeyError, TypeError):
        return f"{name}_bin={index}"


def collect(node):
    """Walk one sample_type.

    Returns (effs, labels): {(wp, flavor, eta_index, pt_index): efficiency} and
    the matching {key: 'eta=[..,..)/pt=[..,..)'} physical-unit bin names.
    """
    effs = {}
    labels = {}
    for wp_item in node["content"]:
        wp = wp_item["key"]
        for flavor_item in wp_item["value"]["content"]:
            flavor = flavor_item["key"]
            eta_binning = flavor_item["value"]
            eta_edges = eta_binning.get("edges")
            for eta_index, pt_binning in enumerate(eta_binning["content"]):
                pt_edges = pt_binning.get("edges")
                for pt_index, value in enumerate(pt_binning["content"]):
                    key = (wp, flavor, eta_index, pt_index)
                    effs[key] = float(value)
                    labels[key] = (
                        bin_label(eta_edges, eta_index, "eta")
                        + "/"
                        + bin_label(pt_edges, pt_index, "pt")
                    )
    return effs, labels


def load_correctionset(path):
    """Read a correctionlib payload, gzipped or not.

    The payload is shipped as ``.json.gz`` and installed under that name, so a
    checker that only reads plain JSON silently ends up validating a different
    file than the one CROWN loads.
    """
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def payload_wps(cset):
    """The working_point keys of an efficiency payload (union over sample types)."""
    return {
        wp_item["key"]
        for correction in cset["corrections"]
        for sample_item in correction["data"]["content"]
        for wp_item in sample_item["value"]["content"]
    }


def btv_thresholds(sf_path):
    """{WP: score cut} from the UParTAK4_wp_values correction of a BTV payload."""
    for correction in load_correctionset(sf_path)["corrections"]:
        if correction["name"] == WP_VALUES:
            return {item["key"]: float(item["value"]) for item in correction["data"]["content"]}
    raise SystemExit(f"{sf_path} has no {WP_VALUES} correction")


def check_file(path, thresholds=None, consumer_wps=None):
    """Return (problems, stats): [(kind, detail)] and what was inspected.

    With *consumer_wps* every sample type must carry exactly those working
    points for every flavour ("WP SET" otherwise).
    """
    cset = load_correctionset(path)

    problems = []
    stats = {
        "sample_types": [],
        "wp_sets": set(),  # one entry per (sample_type, flavour) group
        "flavour_sets": set(),  # one entry per sample_type
        "all_wps": set(),
        "all_flavours": set(),
        "eta_counts": set(),  # eta-bin count per (sample_type, wp, flavour)
        "pt_counts": set(),  # pt-bin count per (sample_type, wp, flavour, eta)
        "cells": 0,
        "comparisons": 0,
    }
    wp_order = order_wps(payload_wps(cset), thresholds)
    stats["wp_order"] = wp_order
    if thresholds:
        for wp in wp_order:
            if wp not in thresholds:
                problems.append(("UNKNOWN WP", f"'{wp}' is not a working point of the BTV payload"))

    for correction in cset["corrections"]:
        for sample_item in correction["data"]["content"]:
            sample_type = sample_item["key"]
            effs, labels = collect(sample_item["value"])

            # Record the shape of what was actually walked for this sample type.
            stats["sample_types"].append(sample_type)
            stats["cells"] += len(effs)
            stats["all_wps"].update(wp for (wp, _, _, _) in effs)
            flavours = {flavor for (_, flavor, _, _) in effs}
            stats["all_flavours"].update(flavours)
            stats["flavour_sets"].add(tuple(sorted(map(str, flavours))))
            per_flavour_wps = {}
            grid = {}
            for (wp, flavor, i_eta, i_pt) in effs:
                per_flavour_wps.setdefault(flavor, set()).add(wp)
                grid.setdefault((wp, flavor), {}).setdefault(i_eta, set()).add(i_pt)
            for wps in per_flavour_wps.values():
                stats["wp_sets"].add(tuple(w for w in wp_order if w in wps))
            if consumer_wps is not None:
                wrong = {flavor: sorted(map(str, wps))
                         for flavor, wps in sorted(per_flavour_wps.items(), key=lambda item: str(item[0]))
                         if wps != set(consumer_wps)}
                if wrong or not per_flavour_wps:
                    problems.append((
                        "WP SET",
                        f"{sample_type}: working points per flavour {wrong or 'none'} "
                        f"!= the consumer's {','.join(consumer_wps)}",
                    ))
            for eta_map in grid.values():
                stats["eta_counts"].add(len(eta_map))
                for pt_indices in eta_map.values():
                    stats["pt_counts"].add(len(pt_indices))

            # Physical-unit name of each (flavour, eta, pt) cell, from whichever
            # working point carries it first.
            cell_labels = {}
            for (wp, flavor, i_eta, i_pt), label in labels.items():
                cell_labels.setdefault((flavor, i_eta, i_pt), label)

            # Invariant 1: finite and in (0, 1]
            for (wp, flavor, i_eta, i_pt), eff in sorted(effs.items()):
                where = f"{sample_type}/{wp}/flavor={flavor}/{labels[(wp, flavor, i_eta, i_pt)]}"
                if not math.isfinite(eff):
                    problems.append(("NOT FINITE", f"{where}: {eff}"))
                elif eff <= 0.0:
                    problems.append(("eff <= 0", f"{where}: {eff}"))
                elif eff > 1.0:
                    problems.append(("eff > 1", f"{where}: {eff}"))

            # Invariants 2 and 3: strict ordering, and 1 - eff(loosest) > 0
            cells = sorted({(f, e, p) for (_, f, e, p) in effs})
            for flavor, i_eta, i_pt in cells:
                where = f"{sample_type}/flavor={flavor}/{cell_labels[(flavor, i_eta, i_pt)]}"
                chain = [(wp, effs[(wp, flavor, i_eta, i_pt)])
                         for wp in wp_order
                         if (wp, flavor, i_eta, i_pt) in effs]
                stats["comparisons"] += max(len(chain) - 1, 0)
                for (loose_wp, loose), (tight_wp, tight) in zip(chain, chain[1:]):
                    if loose - tight <= 0.0:
                        problems.append((
                            "NOT MONOTONIC",
                            f"{where}: eff({loose_wp})={loose} "
                            f"must exceed eff({tight_wp})={tight}",
                        ))
                if chain and 1.0 - chain[0][1] <= 0.0:
                    problems.append((
                        "DEGENERATE",
                        f"{where}: 1 - eff({chain[0][0]})={1.0 - chain[0][1]} "
                        f"must be positive",
                    ))
    return problems, stats


def counted(number, noun):
    return f"{number} {noun}" if number == 1 else f"{number} {noun}s"


def describe(stats):
    """One line naming what was walked, honest about a ragged payload."""
    parts = [counted(len(stats["sample_types"]), "sample type")]

    wp_sets = stats["wp_sets"]
    if len(wp_sets) == 1:
        wps = next(iter(wp_sets))
        parts.append(f"{counted(len(wps), 'WP')} ({','.join(map(str, wps))})")
    elif not wp_sets:
        parts.append("0 WPs")
    else:
        sizes = sorted({len(wps) for wps in wp_sets})
        union = ",".join(w for w in stats["wp_order"] if w in stats["all_wps"])
        parts.append(
            f"WPs RAGGED {sizes[0]}-{sizes[-1]} per flavour (union of "
            f"{len(stats['all_wps'])}: {union})"
        )

    flavour_sets = stats["flavour_sets"]
    if len(flavour_sets) == 1:
        parts.append(counted(len(next(iter(flavour_sets))), "flavour"))
    elif not flavour_sets:
        parts.append("0 flavours")
    else:
        sizes = sorted({len(flavours) for flavours in flavour_sets})
        union = ",".join(sorted(map(str, stats["all_flavours"])))
        parts.append(
            f"flavours RAGGED {sizes[0]}-{sizes[-1]} per sample type (union of "
            f"{len(stats['all_flavours'])}: {union})"
        )

    eta_counts, pt_counts = stats["eta_counts"], stats["pt_counts"]
    if not eta_counts or not pt_counts:
        grid = "no eta/pt bins"
    elif len(eta_counts) == 1 and len(pt_counts) == 1:
        grid = f"{next(iter(eta_counts))} eta x {next(iter(pt_counts))} pt"
    else:
        etas, pts = sorted(eta_counts), sorted(pt_counts)
        grid = (f"RAGGED grid, eta {etas[0]}-{etas[-1]} x pt {pts[0]}-{pts[-1]}")
    parts.append(f"{grid} = {counted(stats['cells'], 'cell')}")

    parts.append(counted(stats["comparisons"], "WP comparison"))
    return "inspected: " + ", ".join(parts)


# ---------------------------------------------------------------------------
# SF-weighted jet-weight check (the fourth failure class)
# ---------------------------------------------------------------------------
# The three value invariants above are necessary but NOT sufficient. The strict
# consumer's last gate is on the jet weight itself: jet_w = num / denom must be
# finite and > 0. For the disjoint score regions the consumer builds
#   passes tightest      : num = SF(tight)*eff(tight),        denom = eff(tight)
#   passes low, not high : num = SF(l)*eff(l) - SF(h)*eff(h), denom = eff(l)-eff(h)
#   fails loosest        : num = 1 - SF(loose)*eff(loose),    denom = 1 - eff(loose)
# so a payload with perfectly monotonic efficiencies still aborts the job when
# the scale factors push a numerator non-positive (that is what the 156
# "non-positive/non-finite jet contribution" failures of sm2018_binned_v1 were).
# The one exception (v7 K1, BTV fixedWP "set these intermediate factors to
# unity"): a NEGATIVE numerator of a "passes low, not high" factor sets the
# factor to 1 and is counted; it is reported here as k1-unity, not as a
# violation.
#
# Below the calibrated pT range of the comb SFs (bc_min_pt, the first pT edge of
# UParTAK4_comb, 30 GeV) the consumer evaluates b/c jets at bc_min_pt + 1e-3
# with twice the deviation of a non-central variation from central
# (btag_sf_strict.cxx, below_bc_range); evaluating at the jet pT instead would
# return the payload's flow value, i.e. the >300 GeV SF.
#
# Flavour dispatch mirrors btag_sf_strict.cxx exactly: b/c (5, 4) use
# UParTAK4_comb, light (0) uses UParTAK4_light; |eta| support is [0, 2.4).

COMB_SF, LIGHT_SF, WP_VALUES = "UParTAK4_comb", "UParTAK4_light", "UParTAK4_wp_values"
ETA_MAX = 2.4


def _payload_edges(cset, axis):
    """Collect every eta (axis=0) or pt (axis=1) bin edge in the payload."""
    edges = set()
    for correction in cset["corrections"]:
        for sample_item in correction["data"]["content"]:
            for wp_item in sample_item["value"]["content"]:
                for flavor_item in wp_item["value"]["content"]:
                    eta_binning = flavor_item["value"]
                    if axis == 0:
                        edges.update(float(e) for e in eta_binning.get("edges", []))
                    else:
                        for pt_binning in eta_binning["content"]:
                            edges.update(float(e) for e in pt_binning.get("edges", []))
    return sorted(edges)


def _sf_binnings(sf_raw, axis, names=(COMB_SF, LIGHT_SF)):
    """Every *axis* binning node of the named SF corrections."""
    nodes = []

    def walk(node):
        if not isinstance(node, dict):
            return
        if node.get("nodetype") == "binning":
            if node.get("input") == axis:
                nodes.append(node)
            for child in node.get("content", []):
                walk(child)
        else:
            for child in node.get("content", []) or []:
                walk(child)
            walk(node.get("value"))

    for correction in sf_raw.get("corrections", []):
        if correction.get("name") in names:
            walk(correction["data"])
    return nodes


def _sf_edges(sf_raw, axis):
    """Collect every *axis* bin edge used by the two UParTAK4 SF corrections."""
    return sorted({float(e) for node in _sf_binnings(sf_raw, axis) for e in node.get("edges", [])})


def bc_min_pt(sf_raw):
    """Lowest calibrated pT of UParTAK4_comb (as btag_payloads.load_upart_comb_min_pt)."""
    lows = {float(node["edges"][0]) for node in _sf_binnings(sf_raw, "pt", (COMB_SF,))}
    if len(lows) != 1:
        raise SystemExit(f"{COMB_SF} has no common lowest pT edge: {sorted(lows)}")
    return lows.pop()


def probe_points(edges, lo, hi):
    """Midpoints of every interval the union of *edges* cuts [lo, hi) into."""
    inner = sorted({lo, hi} | {e for e in edges if lo < e < hi})
    return [(a + b) / 2.0 for a, b in zip(inner[:-1], inner[1:])]


def jet_weights(e, s, wps):
    """Consumer jet weights per score region for one probe point.

    *e*, *s*: {wp: efficiency}, {wp: SF}; *wps*: loosest to tightest. Returns
    [(region, jet_w, k1)] with k1 = True where the K1 rule set the factor to 1.
    """
    tight_to_loose = list(reversed(wps))
    weights = [("pass " + tight_to_loose[0], s[tight_to_loose[0]], False)]
    for high, low in zip(tight_to_loose[:-1], tight_to_loose[1:]):
        denom = e[low] - e[high]
        num = s[low] * e[low] - s[high] * e[high]
        if denom > 0.0 and num < 0.0:
            weights.append((f"{low} not {high}", 1.0, True))
        else:
            weights.append((f"{low} not {high}", num / denom if denom else math.inf, False))
    loosest = tight_to_loose[-1]
    denom = 1.0 - e[loosest]
    weights.append((f"fails {loosest}",
                    (1.0 - s[loosest] * e[loosest]) / denom if denom else math.inf, False))
    return weights


def check_sf_weights(path, sf_path=PINNED_BTV, variations=None, wps=None):
    """Recompute the consumer's jet weight on a grid and report every failure.

    Returns ``(problems, stats)``; *problems* uses the same ``(kind, detail)``
    shape as :func:`check_file` so both feed one report.
    """
    try:
        # Imported lazily so the value checks still run without correctionlib.
        import correctionlib
    except ImportError as error:
        raise SystemExit(
            "--sf-check needs correctionlib; activate the environment that has "
            "it (e.g. the LCG_108 view or `conda activate fake_factors`) or drop "
            "--sf-check to run the efficiency-value checks alone."
        ) from error

    eff_set = correctionlib.CorrectionSet.from_file(path)
    sf_set = correctionlib.CorrectionSet.from_file(sf_path)
    eff = eff_set["btag_efficiency"]

    raw = load_correctionset(path)
    sf_raw = load_correctionset(sf_path)
    thresholds = btv_thresholds(sf_path)
    present = payload_wps(raw)
    wps = order_wps(present if wps is None else wps, thresholds)
    sample_types = [
        item["key"]
        for correction in raw["corrections"]
        for item in correction["data"]["content"]
    ]
    absent = [wp for wp in wps if wp not in present]
    if absent:
        # every probe would fail to evaluate; one violation says it all
        return [("WP SET", f"cannot emulate {','.join(wps)}: the payload has no "
                           f"{','.join(absent)}")], {
            "sample_types": sample_types, "variations": {}, "eta_probes": [],
            "pt_probes": [], "weights_checked": 0, "regions": len(wps) + 1,
            "k1_unity": 0, "lowest": (math.inf, ""), "fails_loosest": {},
            "bc_min_pt": math.nan, "wps": wps,
        }
    min_pt = bc_min_pt(sf_raw)
    # Probe the union of BOTH payloads' bin edges plus bc_min_pt. The consumer
    # evaluates the scale factor at the jet's own pT, so an efficiency bin that
    # straddles two SF bins hides the SF bin with the worse ratio from a grid
    # built out of efficiency edges alone.
    eta_probes = probe_points(
        _payload_edges(raw, 0) + _sf_edges(sf_raw, "abseta"), 0.0, ETA_MAX
    )
    pt_probes = probe_points(
        _payload_edges(raw, 1) + _sf_edges(sf_raw, "pt") + [min_pt], 20.0, 1000.0
    )

    if variations is None:
        variations = {}
        for name in (COMB_SF, LIGHT_SF):
            keys = {
                item["key"]
                for correction in sf_raw["corrections"]
                if correction["name"] == name
                for item in correction["data"]["content"]
                if "key" in item
            }
            variations[name] = sorted(keys)

    problems, checked, k1 = [], 0, 0
    lowest = (math.inf, "")  # smallest positive jet weight and where
    # per sample type: the smallest jet weight of a jet failing the loosest WP,
    # (jet_w, where, eff(loosest), SF ceiling 1/SF(loosest)) -- the margin of
    # the SF ceiling eff < 1/SF, e.g. eff_M(b, > 300 GeV) in 2016preVFP
    fails_loosest = {}
    loosest = wps[0]
    for sample_type in sample_types:
        for flavor in (5, 4, 0):
            name = LIGHT_SF if flavor == 0 else COMB_SF
            sf = sf_set[name]
            for variation in variations[name]:
                for abs_eta in eta_probes:
                    for pt in pt_probes:
                        below = flavor in (5, 4) and pt < min_pt
                        sf_pt = min_pt + 1e-3 if below else pt
                        try:
                            e = {
                                wp: eff.evaluate(sample_type, wp, flavor, abs_eta, pt)
                                for wp in wps
                            }
                            s = {}
                            for wp in wps:
                                value = sf.evaluate(variation, wp, flavor, abs_eta, sf_pt)
                                if below and variation != "central":
                                    central = sf.evaluate("central", wp, flavor, abs_eta, sf_pt)
                                    value = central + 2.0 * (value - central)
                                s[wp] = value
                        except Exception as error:  # unknown key / out of support
                            problems.append((
                                "SF EVAL FAILED",
                                f"{sample_type}/{variation}/flavor={flavor}/"
                                f"eta={abs_eta:g}/pt={pt:g}: {error}",
                            ))
                            continue

                        for region, jet_w, unity in jet_weights(e, s, wps):
                            checked += 1
                            k1 += unity
                            if (region == f"fails {loosest}" and math.isfinite(jet_w)
                                    and jet_w < fails_loosest.get(sample_type, (math.inf,))[0]):
                                fails_loosest[sample_type] = (
                                    jet_w,
                                    f"{variation}/flavor={flavor}/eta={abs_eta:g}/pt={pt:g}",
                                    e[loosest],
                                    1.0 / s[loosest] if s[loosest] > 0.0 else math.inf,
                                )
                            if math.isfinite(jet_w) and 0.0 < jet_w < lowest[0] and not unity:
                                lowest = (jet_w, f"{sample_type}/{variation}/flavor={flavor}/"
                                                 f"eta={abs_eta:g}/pt={pt:g}/[{region}]")
                            if not math.isfinite(jet_w) or jet_w <= 0.0:
                                problems.append((
                                    "JET WEIGHT <= 0",
                                    f"{sample_type}/{variation}/flavor={flavor}/"
                                    f"eta={abs_eta:g}/pt={pt:g}/[{region}]: "
                                    f"jet_w={jet_w:g}",
                                ))
    stats = {
        "sample_types": sample_types,
        "variations": variations,
        "eta_probes": eta_probes,
        "pt_probes": pt_probes,
        "weights_checked": checked,
        "regions": len(wps) + 1,
        "k1_unity": k1,
        "lowest": lowest,
        "fails_loosest": fails_loosest,
        "bc_min_pt": min_pt,
        "wps": wps,
    }
    return problems, stats


def fails_loosest_report(fails_loosest, loosest, warn_below):
    """Lines of the fails-loosest margin, one per group of identical values (a pool)."""
    groups = {}
    for sample_type, entry in fails_loosest.items():
        groups.setdefault(entry, []).append(sample_type)
    lines = []
    for (jet_w, where, eff_value, ceiling), members in sorted(groups.items()):
        flag = "WARNING " if jet_w < warn_below else ""
        lines.append(
            f"    {flag}{','.join(sorted(members))}: jet_w={jet_w:.4g} at {where} "
            f"(eff({loosest})={eff_value:.4g} vs ceiling 1/SF={ceiling:.4g})"
        )
    return lines


def tally(problems):
    counts = {}
    for kind, _ in problems:
        counts[kind] = counts.get(kind, 0) + 1
    ordered = [k for k in KINDS if k in counts] + [k for k in counts if k not in KINDS]
    return "by class: " + ", ".join(f"{kind}: {counts[kind]}" for kind in ordered)


def check_sample_types(path, required):
    """Report every required sample_type the payload does not carry.

    A payload short a sample type does not fail any value invariant -- it fails
    at run time, in the middle of the production, on the first job that runs
    that process. Structural completeness is therefore enforced here rather
    than merely reported.
    """
    present = {
        item["key"]
        for correction in load_correctionset(path)["corrections"]
        for item in correction["data"]["content"]
    }
    return [
        ("MISSING SAMPLE TYPE", f"'{name}' is not in the payload")
        for name in required if name not in present
    ]


def main(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("paths", nargs="*", help="payload .json or .json.gz files")
    parser.add_argument(
        "--era", choices=RUN2_ERAS, default=None,
        help="era of the payload: selects the pinned BTV payload "
             f"({BTV_PIN.format(era='<era>')}), orders the working points by its "
             "thresholds and requires the era's SM sample types",
    )
    parser.add_argument(
        "--sf-check", action="store_true",
        help="also recompute the consumer's jet weight against the BTV scale "
             "factors and fail on any non-positive/non-finite value "
             "(needs correctionlib and cvmfs)",
    )
    parser.add_argument(
        "--sf-payload", default=None,
        help="BTV scale-factor payload for --sf-check (default: the --era pin, "
             f"else {PINNED_BTV})",
    )
    parser.add_argument(
        "--wps", default=None,
        help="comma-separated working points the consumer evaluates; the "
             "payload must carry exactly these and --sf-check emulates them "
             f"(default with --era: {','.join(CONSUMER_WPS)}, the CROWN profile's "
             "btag_sf_working_points; without --era: the payload's own)",
    )
    parser.add_argument(
        "--warn-below", type=float, default=0.05,
        help="--sf-check: flag pools whose smallest jet weight of a jet failing "
             "the loosest WP stays below this (a warning, not a violation; "
             "default 0.05)",
    )
    parser.add_argument(
        "--require-sm-2018", action="store_true",
        help="deprecated (use --era): fail unless the 13 SM sample types of the "
             "2018 v1-v6 payloads are present",
    )
    args = parser.parse_args(argv)

    if not args.paths:
        parser.print_help()
        return 2

    sf_payload = args.sf_payload or (BTV_PIN.format(era=args.era) if args.era else PINNED_BTV)
    thresholds = btv_thresholds(sf_payload) if (args.era or args.sf_check) else None
    required = []
    if args.era:
        required = REQUIRED_SAMPLE_TYPES[args.era]
    elif args.require_sm_2018:
        required = SM_2018_SAMPLE_TYPES
    if args.wps:
        consumer_wps = tuple(wp.strip() for wp in args.wps.split(",") if wp.strip())
    else:
        consumer_wps = CONSUMER_WPS if args.era else None
    print(f"era: {args.era or 'not given'}; BTV payload: "
          f"{sf_payload if thresholds else 'not used'}; consumer WPs: "
          f"{','.join(consumer_wps) if consumer_wps else 'those of the payload'}")

    total = 0
    for path in args.paths:
        problems, stats = check_file(path, thresholds, consumer_wps)
        if required:
            problems = check_sample_types(path, required) + problems
        sf_stats = None
        if args.sf_check:
            sf_problems, sf_stats = check_sf_weights(
                path, sf_payload, wps=list(consumer_wps) if consumer_wps else stats["wp_order"]
            )
            problems = problems + sf_problems

        total += len(problems)
        print(f"=== {path}: {len(problems)} violation(s) ===")
        print(describe(stats))
        if sf_stats:
            n_var = sum(len(v) for v in sf_stats["variations"].values())
            print(
                f"sf-check: {sf_stats['weights_checked']} jet weights over "
                f"{len(sf_stats['sample_types'])} sample types x {n_var} "
                f"variations x {len(sf_stats['eta_probes'])} eta x "
                f"{len(sf_stats['pt_probes'])} pt probe points x "
                f"{sf_stats['regions']} score regions "
                f"(WPs {','.join(sf_stats['wps'])}, b/c below "
                f"{sf_stats['bc_min_pt']:g} GeV as the consumer); "
                f"k1-unity: {sf_stats['k1_unity']}"
            )
            print(f"sf-check: smallest positive jet weight {sf_stats['lowest'][0]:.4g} "
                  f"at {sf_stats['lowest'][1]}")
            if sf_stats["fails_loosest"]:
                loosest = sf_stats["wps"][0]
                print(f"sf-check: smallest jet weight failing {loosest} per pool "
                      f"(WARNING below {args.warn_below:g}):")
                for line in fails_loosest_report(sf_stats["fails_loosest"], loosest, args.warn_below):
                    print(line)
        if problems:
            print(tally(problems))
        for kind, detail in problems[:50]:
            print(f"    {kind:<19} {detail}")
        if len(problems) > 50:
            print(f"    ... and {len(problems) - 50} more")

    print("\nPAYLOAD IS SAFE TO SHIP" if total == 0
          else f"\nDO NOT SHIP: {total} violation(s) the strict consumer would throw on")
    return 1 if total else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
