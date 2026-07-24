"""
Script for calculating b-tagging efficiency from preselected ROOT RDataFrames.

For each sample type (process) the efficiency is computed as

    eff(wp, flavor, pt, eta) = N(jets passing WP) / N(total jets)

binned in jet pT and jet |eta|, for each b-tagging working point and jet
flavour category.

The result is written as a correctionlib JSON (+ gzipped version) with the
following nested lookup structure:

    Category(sample_type)
      -> Category(working_point)
           -> Category(jet_flavor)
                -> Binning(jet_eta)
                     -> Binning(jet_pt)
                          -> float  (efficiency)

Usage
-----
    python btag_efficiency.py --config-file configs/btag_efficiency/2024/btag_efficiency.yaml

The config may define either a single ``channel`` or a shared ``channels`` list.
"""

import argparse
import concurrent.futures
import copy
import glob
import gzip
import json
import logging
import multiprocessing
import os
import shutil
import subprocess
from typing import Dict, List, Optional, Tuple

import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import correctionlib.schemav2 as cs

import CustomLogging as logging_helper
import helper.functions as func
from helper.btag_accumulators import (
    WeightedEfficiency,
    effective_population,
    unweighted_efficiency,
    weighted_efficiency,
)
from helper import btag_validation as bval
from helper.correctionlib_json import write_json

hep.style.use(hep.style.CMS)
plt.rcParams["axes.linewidth"] = 1.0 # set non bold axes lines

# Track per-bin sum of squared weights on every histogram so that
# GetBinError()**2 == sumw2 for the weighted efficiency accumulators below.
ROOT.TH1.SetDefaultSumw2(True)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
def build_arg_parser() -> argparse.ArgumentParser:
    """Build the command line argument parser for the btag efficiency calculator."""
    parser = argparse.ArgumentParser(
        description="Calculate b-tagging efficiency and write a correctionlib JSON."
    )
    parser.add_argument(
        "--config-file",
        required=True,
        help="Path to the btag efficiency config YAML file.",
    )
    parser.add_argument(
        "--file-path",
        default=None,
        help="Override the 'file_path' config setting (input preselection directory). Also "
        "settable via the TFF_FILE_PATH environment variable. Precedence: this CLI argument "
        "> TFF_FILE_PATH > config file.",
    )
    parser.add_argument(
        "--output-path",
        default=None,
        help="Override the 'output_base' config setting (base directory for calculator "
        "outputs, replacing the literal 'workdir'; default: 'workdir'). Also settable via "
        "the TFF_OUTPUT_PATH environment variable. Precedence: this CLI argument > "
        "TFF_OUTPUT_PATH > config file > default.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (default: INFO).",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of parallel worker processes for inter-sample parallelization (default: 1).",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=4,
        help="ROOT threads per worker for RDataFrame implicit MT (0 = all cores, default: 0).",
    )
    return parser


parser = build_arg_parser()


# ---------------------------------------------------------------------------
# B-tag working-point resolution (pinned BTV payload validation)
# ---------------------------------------------------------------------------
# Tolerance for comparing the payload working points against the config block.
_WP_TOL = 1e-9

# Name of the UParTAK4 working-point correction inside the BTV payload.
WP_VALUES_CORRECTION = "UParTAK4_wp_values"


def load_working_points_from_source(source_path: str) -> Dict[str, float]:
    """Read ``UParTAK4_wp_values`` from a gzipped BTV correctionlib payload.

    Uses only the standard library (``gzip`` + ``json``); ``correctionlib`` is
    intentionally not needed just to read the frozen working-point cut values.
    Walks ``corrections[name == "UParTAK4_wp_values"].data.content`` collecting
    the ``key`` -> ``value`` pairs (e.g. ``{"L": 0.0308, "M": 0.161, ...}``).

    Raises:
        RuntimeError: The payload is missing, unreadable, or contains no
            ``UParTAK4_wp_values`` correction.
    """
    try:
        with gzip.open(source_path, "rt") as handle:
            payload = json.load(handle)
    except FileNotFoundError as error:
        raise RuntimeError(
            f"btag_working_points_source '{source_path}' not found; cannot "
            f"validate the configured working points against the pinned BTV payload."
        ) from error
    except OSError as error:
        raise RuntimeError(
            f"btag_working_points_source '{source_path}' could not be read as a "
            f"gzipped correctionlib payload: {error}"
        ) from error

    for correction in payload.get("corrections", []):
        if correction.get("name") == WP_VALUES_CORRECTION:
            return {
                item["key"]: float(item["value"])
                for item in correction["data"]["content"]
            }
    found = [c.get("name") for c in payload.get("corrections", [])]
    raise RuntimeError(
        f"correction '{WP_VALUES_CORRECTION}' not found in "
        f"btag_working_points_source '{source_path}'; found {found}."
    )


def resolve_working_points(config: Dict) -> Dict[str, float]:
    """Return the validated b-tag working points for the run.

    When ``btag_working_points_source`` is set in the config, the
    ``UParTAK4_wp_values`` read from that pinned BTV payload are REQUIRED to
    agree with the ``btag_working_points`` config block: both must be present
    and every working point must match within :data:`_WP_TOL`, else the run is
    aborted fatally (a silent upstream re-derivation of the discriminator cut
    values must never slip into the efficiency payload). Without a source the
    configured ``btag_working_points`` are returned unchanged (Run-3 behaviour).

    Returns the configured working-point dict (unchanged) on success.

    Raises:
        RuntimeError: A source is set but ``btag_working_points`` is missing, or
            the working-point sets/values disagree with the pinned payload.
    """
    configured = config.get("btag_working_points")
    source = config.get("btag_working_points_source")
    if not source:
        return configured

    if not configured:
        raise RuntimeError(
            "btag_working_points_source is set but no 'btag_working_points' block "
            "is present in the config; both must be present so the configured cut "
            "values can be validated against the pinned BTV payload."
        )

    payload_wps = load_working_points_from_source(source)

    configured_keys = set(configured)
    payload_keys = set(payload_wps)
    if configured_keys != payload_keys:
        raise RuntimeError(
            f"b-tag working-point sets disagree with the pinned BTV payload "
            f"'{source}': config has {sorted(configured_keys)}, payload has "
            f"{sorted(payload_keys)}. Revalidate the working points before use."
        )

    mismatches = {
        wp: (float(configured[wp]), payload_wps[wp])
        for wp in configured
        if abs(float(configured[wp]) - payload_wps[wp]) > _WP_TOL
    }
    if mismatches:
        raise RuntimeError(
            f"b-tag working points disagree with the pinned BTV payload '{source}' "
            f"(config vs payload): {mismatches}. Revalidate the working points and "
            f"update the config block before use."
        )

    return configured


# ---------------------------------------------------------------------------
# File discovery
# ---------------------------------------------------------------------------
def get_process_files(
    config: Dict,
    process: str,
) -> List[str]:
    """Return all ROOT files that belong to *process*.

    The function searches for files matching
        {file_path}/preselection/{era}/{channel}/{process}*.root
    for every configured channel.

    Args:
        config: Loaded configuration dictionary.
        process: Process / sample-type name (e.g. ``"ttbar"``).

    Returns:
        List of matching file paths (may be empty).
    """
    file_path = config["file_path"]
    channels = get_config_channels(config)
    era = config.get("era", None)

    found: List[str] = []
    if era is None:
        raise ValueError(
            "Config must specify an 'era' for file discovery."
        )
    for channel in channels:
        pattern = os.path.join(
            file_path,
            "preselection",
            era,
            channel,
            f"{process}*.root",
        )
        found.extend(glob.glob(pattern))

    return sorted(set(found))


def get_config_channels(config: Dict) -> List[str]:
    """Return the list of channels configured for this run.

    Supports both the legacy single-channel key ``channel`` and the new
    multi-channel key ``channels`` for a shared config file.
    """
    channels = config.get("channels")
    if channels is not None:
        if not isinstance(channels, list) or not channels:
            raise ValueError("Config key 'channels' must be a non-empty list.")
        return [str(channel) for channel in channels]

    channel = config.get("channel")
    if channel is None:
        raise ValueError("Config must specify either 'channel' or 'channels'.")
    return [str(channel)]


def get_channel_label(config: Dict) -> str:
    """Return a compact channel label for logging, plots, and output paths."""
    channels = get_config_channels(config)
    return channels[0] if len(channels) == 1 else "all_channels"


def get_channel_display(config: Dict) -> str:
    """Return a human-readable channel label."""
    return ", ".join(get_config_channels(config))


def get_channel_configs(config: Dict) -> List[Dict]:
    """Return single-channel config views for the configured channels.

    A shared config containing ``channels`` is expanded into one config per
    channel so that each run writes to its own output directory.
    """
    channel_configs: List[Dict] = []
    for channel in get_config_channels(config):
        channel_config = copy.deepcopy(config)
        channel_config["channel"] = channel
        channel_config.pop("channels", None)
        channel_configs.append(channel_config)
    return channel_configs


# ---------------------------------------------------------------------------
# Binning helpers
# ---------------------------------------------------------------------------
def get_process_bins(
    config: Dict,
    process: str,
) -> Tuple[List[float], List[float]]:
    """Return (pt_bins, eta_bins) for *process*.

    Per-process overrides in ``processes[process]`` take
    precedence; otherwise the global ``jet_pt_bins`` / ``jet_eta_bins`` are used.
    """
    proc_conf = config["processes"].get(process, {})
    pt_bins = proc_conf.get("jet_pt_bins", config["jet_pt_bins"])
    eta_bins = proc_conf.get("jet_eta_bins", config["jet_eta_bins"])
    return list(pt_bins), list(eta_bins)


# ---------------------------------------------------------------------------
# Input-branch validation
# ---------------------------------------------------------------------------
def required_input_columns(config: Dict) -> List[str]:
    """Return the configured columns that must be present in every input tree.

    The four probe jet columns plus the per-event weight column (the weighted
    efficiency is the production value, so the weight column is mandatory).
    """
    return [
        config["jet_pt_column"],
        config["jet_eta_column"],
        config["jet_flavor_column"],
        config["jet_btag_column"],
        config.get("weight_column", "weight"),
    ]


def _tree_columns(file_path: str, tree_name: str) -> List[str]:
    """Return the branch names of ``tree_name`` in ``file_path`` (no event loop)."""
    fin = ROOT.TFile.Open(file_path)
    if not fin or fin.IsZombie():
        raise RuntimeError(f"Could not open input file '{file_path}'.")
    try:
        tree = fin.Get(tree_name)
        if not tree:
            raise RuntimeError(
                f"Tree '{tree_name}' not found in input file '{file_path}'."
            )
        return [str(branch.GetName()) for branch in tree.GetListOfBranches()]
    finally:
        fin.Close()


def validate_process_branches(files: List[str], config: Dict) -> None:
    """Validate that every configured column is present in *every* input file.

    Wired into the per-process file loop before any RDataFrame graph is built;
    raises :class:`RuntimeError` listing the missing names per file.
    """
    tree_name = config["tree"]
    columns_by_file = {f: _tree_columns(f, tree_name) for f in files}
    bval.validate_input_branches(columns_by_file, required_input_columns(config))


# ---------------------------------------------------------------------------
# Efficiency calculation
# ---------------------------------------------------------------------------
def calculate_efficiency_histograms(
    chain: ROOT.TChain,
    config: Dict,
    process: str,
    pt_bins: List[float],
    eta_bins: List[float],
) -> Tuple[
    Dict[str, Dict[str, Tuple[ROOT.TH2D, ROOT.TH2D]]],
    Dict[str, Dict[str, WeightedEfficiency]],
]:
    """Fill 2-D (pT, |eta|) histograms for every (WP, flavour) combination.

    For each combination two flavours of histogram are filled from the same
    masked jet columns:

    - *unweighted* ``(h_pass, h_total)`` -- plain jet counts, kept for the count
      panels of the plots and for the unweighted diagnostic;
    - *weighted* ``(h_pass_w, h_total_w)`` -- filled with the per-event weight
      (broadcast onto every jet of the event), tracking ``sumw`` in the bin
      content and ``sumw2`` in ``GetBinError()**2`` (``TH1::SetDefaultSumw2`` is
      enabled at import).

    The weighted sums, together with the raw counts, are packed into a
    :class:`WeightedEfficiency` accumulator per (WP, flavour). The weighted
    efficiency is the production value; the unweighted counts drive the
    diagnostic.

    Args:
        chain:   ROOT TChain with all files for this process already attached.
        config:  Loaded configuration dictionary.
        process: Process name (used only for histogram titles / logging).

    Returns:
        Tuple ``(histograms, accumulators)`` where
        ``histograms[wp][flavor] = (h_pass, h_total)`` are the *unweighted*
        count histograms and ``accumulators[wp][flavor]`` is a
        :class:`WeightedEfficiency`.

    Raises:
        RuntimeError: The RDataFrame is empty, or the configured weight column
            is missing from the input tree (the weighted path never silently
            falls back to unweighted).
    """
    log = logging.getLogger("btag_efficiency")

    jet_pt_col = config["jet_pt_column"]
    jet_eta_col = config["jet_eta_column"]
    jet_flavor_col = config["jet_flavor_column"]
    jet_btag_col = config["jet_btag_column"]
    weight_col = config.get("weight_column", "weight")
    jet_sel = config.get("jet_selection", "").strip()

    wps: Dict[str, float] = config["btag_working_points"]
    flavors: Dict[str, int] = config["jet_flavor_categories"]

    n_pt = len(pt_bins) - 1
    n_eta = len(eta_bins) - 1
    pt_arr = np.array(pt_bins, dtype=float)
    eta_arr = np.array(eta_bins, dtype=float)

    rdf = ROOT.RDataFrame(chain)

    if func.rdf_is_empty(rdf):
        raise RuntimeError(f"Empty RDataFrame for process '{process}'.")

    # The weighted efficiency is the production value, so a missing weight
    # column is a hard error rather than a silent fall-back to unweighted
    # counts. (Task 18 adds full branch validation; until then this is the
    # guard.)
    column_names = {str(name) for name in rdf.GetColumnNames()}
    if weight_col not in column_names:
        raise RuntimeError(
            f"Weight column '{weight_col}' not found in the input tree for "
            f"process '{process}'. The weighted b-tag efficiency requires this "
            f"per-event weight column (snapshotted by preselection.py); refusing "
            f"to fall back to unweighted counts."
        )

    # Compile an equal-length + finiteness assertion into the graph: the first
    # event with unequal-length probe vectors or a non-finite weight/jet value
    # aborts the event loop rather than silently corrupting the accumulators.
    rdf = bval.validate_probe_vectors(
        rdf,
        {
            "pt": jet_pt_col,
            "eta": jet_eta_col,
            "flavor": jet_flavor_col,
            "btag": jet_btag_col,
            "weight": weight_col,
        },
    )

    # Optionally apply a per-jet pre-selection mask so that only jets within
    # the configured pT / eta acceptance are counted.
    #
    # Because the jet columns are RVec<float> / RVec<int> we use element-wise
    # masking rather than a row-level Filter.
    if jet_sel:
        # Build a boolean mask RVec and apply it to every jet column we use.
        rdf = rdf.Define(
            "_jet_sel_mask",
            f"({jet_sel})",
        )
        jet_pt_col_use = "_jet_pt_sel"
        jet_eta_col_use = "_jet_eta_sel"
        jet_flavor_col_use = "_jet_flavor_sel"
        jet_btag_col_use = "_jet_btag_sel"
        rdf = (
            rdf.Define(jet_pt_col_use, f"{jet_pt_col}[_jet_sel_mask]")
               .Define(jet_eta_col_use, f"{jet_eta_col}[_jet_sel_mask]")
               .Define(jet_flavor_col_use, f"{jet_flavor_col}[_jet_sel_mask]")
               .Define(jet_btag_col_use, f"{jet_btag_col}[_jet_sel_mask]")
        )
    else:
        jet_pt_col_use = jet_pt_col
        jet_eta_col_use = jet_eta_col
        jet_flavor_col_use = jet_flavor_col
        jet_btag_col_use = jet_btag_col

    histograms: Dict[str, Dict[str, Tuple]] = {}
    accumulators: Dict[str, Dict[str, WeightedEfficiency]] = {}

    for wp_name, wp_cut in wps.items():
        histograms[wp_name] = {}
        accumulators[wp_name] = {}
        for flavor_name, flavor_id in flavors.items():
            log.debug(
                f"  Filling histograms: process={process}, WP={wp_name}, "
                f"flavor={flavor_name} (hadronFlavour=={flavor_id})"
            )
            # Masked RVec columns for this flavour
            flavor_mask = f"({jet_flavor_col_use} == {flavor_id})"
            pass_mask = f"({jet_flavor_col_use} == {flavor_id} && {jet_btag_col_use} >= {wp_cut})"

            col_pt_all = f"_pt_{flavor_name}_all"
            col_eta_all = f"_eta_{flavor_name}_all"
            col_pt_pass = f"_pt_{flavor_name}_pass_{wp_name}"
            col_eta_pass = f"_eta_{flavor_name}_pass_{wp_name}"
            # Per-jet weight RVecs: broadcast the scalar event weight onto every
            # selected jet so Histo2D fills each jet with the event weight.
            col_w_all = f"_w_{flavor_name}_all"
            col_w_pass = f"_w_{flavor_name}_pass_{wp_name}"

            rdf_loc = (
                rdf.Define(col_pt_all, f"{jet_pt_col_use}[{flavor_mask}]")
                   .Define(col_eta_all, f"abs({jet_eta_col_use}[{flavor_mask}])")
                   .Define(col_pt_pass, f"{jet_pt_col_use}[{pass_mask}]")
                   .Define(col_eta_pass, f"abs({jet_eta_col_use}[{pass_mask}])")
                   .Define(
                       col_w_all,
                       f"ROOT::VecOps::RVec<double>({col_pt_all}.size(), (double){weight_col})",
                   )
                   .Define(
                       col_w_pass,
                       f"ROOT::VecOps::RVec<double>({col_pt_pass}.size(), (double){weight_col})",
                   )
            )

            h_total = rdf_loc.Histo2D(
                ROOT.RDF.TH2DModel(
                    f"h_total_{process}_{flavor_name}_{wp_name}",
                    f"{process}: all {flavor_name} jets; jet p_{{T}} [GeV]; jet |#eta|",
                    n_pt, pt_arr, n_eta, eta_arr,
                ),
                col_pt_all,
                col_eta_all,
            )
            h_pass = rdf_loc.Histo2D(
                ROOT.RDF.TH2DModel(
                    f"h_pass_{process}_{flavor_name}_{wp_name}",
                    f"{process}: {flavor_name} jets passing {wp_name}; jet p_{{T}} [GeV]; jet |#eta|",
                    n_pt, pt_arr, n_eta, eta_arr,
                ),
                col_pt_pass,
                col_eta_pass,
            )
            h_total_w = rdf_loc.Histo2D(
                ROOT.RDF.TH2DModel(
                    f"h_total_w_{process}_{flavor_name}_{wp_name}",
                    f"{process}: all {flavor_name} jets (weighted); jet p_{{T}} [GeV]; jet |#eta|",
                    n_pt, pt_arr, n_eta, eta_arr,
                ),
                col_pt_all,
                col_eta_all,
                col_w_all,
            )
            h_pass_w = rdf_loc.Histo2D(
                ROOT.RDF.TH2DModel(
                    f"h_pass_w_{process}_{flavor_name}_{wp_name}",
                    f"{process}: {flavor_name} jets passing {wp_name} (weighted); jet p_{{T}} [GeV]; jet |#eta|",
                    n_pt, pt_arr, n_eta, eta_arr,
                ),
                col_pt_pass,
                col_eta_pass,
                col_w_pass,
            )

            # Trigger the event loop and clone the count histograms so they
            # outlive the lazy RDataFrame action (used by the plot count panels).
            h_pass_v = h_pass.GetValue().Clone()
            h_total_v = h_total.GetValue().Clone()
            h_pass_w_v = h_pass_w.GetValue()
            h_total_w_v = h_total_w.GetValue()

            histograms[wp_name][flavor_name] = (h_pass_v, h_total_v)
            accumulators[wp_name][flavor_name] = _accumulator_from_histograms(
                h_pass_w_v, h_total_w_v, h_pass_v, h_total_v, n_pt, n_eta
            )

    return histograms, accumulators


def _accumulator_from_histograms(
    h_pass_w: ROOT.TH2D,
    h_total_w: ROOT.TH2D,
    h_pass_raw: ROOT.TH2D,
    h_total_raw: ROOT.TH2D,
    n_pt: int,
    n_eta: int,
) -> WeightedEfficiency:
    """Pack weighted (sumw/sumw2) and raw counts into a WeightedEfficiency.

    ``sumw`` is the weighted bin content and ``sumw2`` is ``GetBinError()**2``
    (valid because ``TH1::SetDefaultSumw2`` is enabled). Arrays are laid out as
    ``[i_pt][i_eta]`` to match the efficiency-grid convention.

    ``raw_pt_overflow`` is read from ROOT's pt-*overflow* bin (index
    ``n_pt + 1``) of the raw/unweighted total histogram, summed over every
    in-range eta row (``1..n_eta``). ``Histo2D`` always tracks the overflow
    bin for entries above the top edge even though it is not displayed or
    included in ``raw_total`` -- this is the only way to see real jets above
    the top pt edge (e.g. pt > 1000 GeV for the usual binnings).
    """
    sumw_total = np.zeros((n_pt, n_eta))
    sumw2_total = np.zeros((n_pt, n_eta))
    sumw_pass = np.zeros((n_pt, n_eta))
    sumw2_pass = np.zeros((n_pt, n_eta))
    raw_total = np.zeros((n_pt, n_eta))
    raw_pass = np.zeros((n_pt, n_eta))

    for i_pt in range(1, n_pt + 1):
        for i_eta in range(1, n_eta + 1):
            sumw_total[i_pt - 1][i_eta - 1] = h_total_w.GetBinContent(i_pt, i_eta)
            sumw2_total[i_pt - 1][i_eta - 1] = h_total_w.GetBinError(i_pt, i_eta) ** 2
            sumw_pass[i_pt - 1][i_eta - 1] = h_pass_w.GetBinContent(i_pt, i_eta)
            sumw2_pass[i_pt - 1][i_eta - 1] = h_pass_w.GetBinError(i_pt, i_eta) ** 2
            raw_total[i_pt - 1][i_eta - 1] = h_total_raw.GetBinContent(i_pt, i_eta)
            raw_pass[i_pt - 1][i_eta - 1] = h_pass_raw.GetBinContent(i_pt, i_eta)

    raw_pt_overflow = sum(
        h_total_raw.GetBinContent(n_pt + 1, i_eta) for i_eta in range(1, n_eta + 1)
    )

    return WeightedEfficiency(
        sumw_total=sumw_total,
        sumw2_total=sumw2_total,
        sumw_pass=sumw_pass,
        sumw2_pass=sumw2_pass,
        raw_total=raw_total,
        raw_pass=raw_pass,
        raw_pt_overflow=raw_pt_overflow,
    )


def accumulators_to_efficiencies(
    accumulators: Dict[str, Dict[str, WeightedEfficiency]],
) -> Tuple[
    Dict[str, Dict[str, List[List[float]]]],
    Dict[str, Dict[str, List[List[float]]]],
]:
    """Convert accumulators to production (weighted) efficiency + uncertainty grids.

    The production efficiency is the weighted value ``sumw_pass / sumw_total`` and
    the displayed uncertainty is ``sqrt(var)`` from :func:`weighted_efficiency`.
    Empty bins (``nan``) are serialised as ``0.0`` here so the correctionlib JSON
    and plots stay finite -- matching the previous empty-bin behaviour. The raw
    ``nan``-carrying arrays are preserved in the diagnostic sidecar. No clamping
    of out-of-``[0, 1]`` values is applied (that is Task 18's gating).

    Returns:
        Tuple of two dicts ``[wp][flavor] -> [i_pt][i_eta]`` for efficiencies and
        their (weighted) uncertainties.
    """
    efficiencies: Dict[str, Dict[str, List[List[float]]]] = {}
    uncertainties: Dict[str, Dict[str, List[List[float]]]] = {}
    for wp_name, flavor_accs in accumulators.items():
        efficiencies[wp_name] = {}
        uncertainties[wp_name] = {}
        for flavor_name, acc in flavor_accs.items():
            eff, var = weighted_efficiency(acc)
            unc = np.sqrt(var)
            eff = np.where(np.isfinite(eff), eff, 0.0)
            unc = np.where(np.isfinite(unc), unc, 0.0)
            efficiencies[wp_name][flavor_name] = eff.tolist()
            uncertainties[wp_name][flavor_name] = unc.tolist()

    return efficiencies, uncertainties


# ---------------------------------------------------------------------------
# correctionlib JSON construction
# ---------------------------------------------------------------------------
def _build_pt_binning(
    pt_bins: List[float],
    pt_efficiencies: List[float],
) -> cs.Binning:
    """Leaf Binning node over jet pT."""
    return cs.Binning(
        nodetype="binning",
        input="jet_pt",
        edges=list(pt_bins),
        content=pt_efficiencies,
        flow="clamp",
    )


def _build_eta_binning(
    eta_bins: List[float],
    pt_bins: List[float],
    eff_grid: List[List[float]],
) -> cs.Binning:
    """Outer Binning node over jet |eta|; one pt-Binning per eta bin."""
    n_eta = len(eta_bins) - 1
    n_pt = len(pt_bins) - 1
    return cs.Binning(
        nodetype="binning",
        input="jet_eta",
        edges=list(eta_bins),
        content=[
            _build_pt_binning(pt_bins, [eff_grid[i_pt][i_eta] for i_pt in range(n_pt)])
            for i_eta in range(n_eta)
        ],
        flow="error",
    )


def _build_flavor_category(
    flavors: Dict[str, int],
    flavor_bins: Dict[str, Tuple[List[float], List[float]]],
    wp_efficiencies: Dict[str, List[List[float]]],
) -> cs.Category:
    """Category node over jet flavour labels.

    Each flavour carries its own ``(pt_bins, eta_bins)`` (from
    ``flavor_bins[flavor]``): after Task-18 bin merging the binning is frozen
    per ``(sample_type, flavor)`` and may differ between flavours.
    """
    return cs.Category(
        nodetype="category",
        input="jet_flavor",
        content=[
            cs.CategoryItem(
                key=flavors[flavor_name],
                value=_build_eta_binning(
                    flavor_bins[flavor_name][1],
                    flavor_bins[flavor_name][0],
                    wp_efficiencies[flavor_name],
                ),
            )
            for flavor_name in flavors
        ],
    )


def _build_wp_category(
    wps: Dict[str, float],
    flavors: Dict[str, int],
    flavor_bins: Dict[str, Tuple[List[float], List[float]]],
    sample_efficiencies: Dict[str, Dict[str, List[List[float]]]],
) -> cs.Category:
    """Category node over b-tag working points.

    All working points of a ``(sample_type, flavor)`` share the same binning
    (``flavor_bins[flavor]``) so working-point comparability is preserved.
    """
    return cs.Category(
        nodetype="category",
        input="working_point",
        content=[
            cs.CategoryItem(
                key=wp_name,
                value=_build_flavor_category(
                    flavors, flavor_bins, sample_efficiencies[wp_name]
                ),
            )
            for wp_name in wps
        ],
    )


# Modern, professional color palette for eta bins (works well for 4+ bins)
_ETA_COLORS = [
    "#2ca02c",  # forest green
    "#9467bd",  # purple
    "#bcbd22",  # olive
    "#e377c2",  # pink
]
 
_FLAVOR_LABELS = {
    "b":     r"$b$-jets",
    "c":     r"$c$-jets",
    "light": r"light jets",
    "udsg":  r"light jets",
    "bc":    r"$b/c$-jets",
}

_CHANNEL_LABELS = {
    "em":    r"$e\mu$",
    "et":    r"$e\tau_h$",
    "mt":    r"\mu\tau_h",
    "tt":    r"\tau_h\tau_h",
    "ee":    r"ee",
    "mm":    r"$\mu\mu$",
} 

def plot_histograms_and_efficiencies(
    histograms: Dict[str, Dict[str, Tuple[ROOT.TH2D, ROOT.TH2D]]],
    efficiencies: Dict[str, Dict[str, List[List[float]]]],
    uncertainties: Dict[str, Dict[str, List[List[float]]]],
    pt_bins: List[float],
    eta_bins: List[float],
    process: str,
    config: Dict,
    output_path: str,
) -> None:
    """Produce one figure per working point for *process*.
 
    Each figure has one column per jet flavour. Within each column there are
    two stacked panels:
 
    - **Upper (70%)**: raw jet counts (h_total solid, h_pass dashed) vs jet pT,
      one colour per |η| bin.
    - **Lower (30%)**: b-tag efficiency ± stat. uncertainty vs jet pT,
      one colour per |η| bin.
 
    The CMS label is placed at the top of each upper panel. The process/era/channel
    annotation is placed inside the upper panel as a text box in the top-left corner.
    Flavour labels are placed in the top-left of each subplot.
    """
    log = logging.getLogger("btag_efficiency")
 
    wps:     Dict[str, float] = config["btag_working_points"]
    flavors: Dict[str, int]   = config["jet_flavor_categories"]
    flavor_list = list(flavors.keys())
    n_flavors   = len(flavor_list)
 
    n_eta = len(eta_bins) - 1
    n_pt  = len(pt_bins)  - 1
 
    eta_labels = [
        f"{eta_bins[i]:.1f}–{eta_bins[i+1]:.1f}"
        for i in range(n_eta)
    ]
    pt_centres = np.array(
        [(pt_bins[i] + pt_bins[i + 1]) / 2.0 for i in range(n_pt)]
    )
 
    era     = config.get("era", "")
    channel = get_channel_display(config)
    center_of_mass = config.get("center_of_mass", 13.6)

    plot_dir = os.path.join(output_path, "plots")
    os.makedirs(plot_dir, exist_ok=True)
 
    # ================================================================== #
    #  One figure per working point, one column per flavour             #
    # ================================================================== #
    for wp in wps:
        fig, axes = plt.subplots(
            2, n_flavors,
            figsize=(4.5 * n_flavors, 5.0),
            gridspec_kw={
                "height_ratios": [7, 3],  # 70% / 30% split
                "hspace": 0.05,           # minimal vertical spacing
                "wspace": 0.30,           # horizontal spacing between columns
            },
            sharex="col",
        )
        # Normalise axes indexing so axes[row][col] always works
        if n_flavors == 1:
            axes = np.array(axes).reshape(2, 1)
 
        for i_flav, flavor in enumerate(flavor_list):
            ax_top = axes[0, i_flav]
            ax_bot = axes[1, i_flav]
 
            h_pass, h_total = histograms[wp][flavor]
            eff_grid = efficiencies[wp][flavor]
            unc_grid = uncertainties[wp][flavor]
            max_count = 0.0
            eff_y_min = np.inf
            eff_y_max = -np.inf
 
            # ============================================================ #
            #  Fill upper and lower panels                                 #
            # ============================================================ #
            for i_eta in range(n_eta):
                root_bin = i_eta + 1
                proj_total = h_total.ProjectionX(
                    f"_proj_total_{process}_{flavor}_{wp}_{i_eta}",
                    root_bin, root_bin,
                )
                proj_pass = h_pass.ProjectionX(
                    f"_proj_pass_{process}_{flavor}_{wp}_{i_eta}",
                    root_bin, root_bin,
                )
 
                total_vals = np.array(
                    [proj_total.GetBinContent(b) for b in range(1, n_pt + 1)]
                )
                pass_vals = np.array(
                    [proj_pass.GetBinContent(b) for b in range(1, n_pt + 1)]
                )
                eff_vals = np.array([eff_grid[i_pt][i_eta] for i_pt in range(n_pt)])
                unc_vals = np.array([unc_grid[i_pt][i_eta] for i_pt in range(n_pt)])
                if total_vals.size:
                    max_count = max(max_count, float(np.max(total_vals)))
                if pass_vals.size:
                    max_count = max(max_count, float(np.max(pass_vals)))

                eff_low = eff_vals - unc_vals
                eff_high = eff_vals + unc_vals
                finite_eff = np.isfinite(eff_low) & np.isfinite(eff_high)
                if np.any(finite_eff):
                    eff_y_min = min(eff_y_min, float(np.min(eff_low[finite_eff])))
                    eff_y_max = max(eff_y_max, float(np.max(eff_high[finite_eff])))
 
                color     = _ETA_COLORS[i_eta % len(_ETA_COLORS)]
                eta_label = (
                    r"$|\eta| \in ["
                    + eta_labels[i_eta]
                    + r"]$"
                )
 
                # upper panel — counts (solid for total, dashed for pass)
                step_x = np.append(pt_bins[:-1], pt_bins[-1])
                ax_top.step(
                    step_x, np.append(total_vals, total_vals[-1]),
                    where="post", color=color, linewidth=2.0, linestyle="-",
                    label=rf"total  {eta_label}",
                )
                ax_top.step(
                    step_x, np.append(pass_vals, pass_vals[-1]),
                    where="post", color=color, linewidth=2.0, linestyle="--",
                    label=rf"pass  {eta_label}",
                )
 
                # lower panel — efficiency with error bars
                ax_bot.step(
                    step_x, np.append(eff_vals, eff_vals[-1]),
                    where="post", color=color, linewidth=1.8,
                    label=eta_label,
                )
                ax_bot.errorbar(
                    pt_centres, eff_vals, yerr=unc_vals,
                    fmt="o", color=color, markersize=4.0,
                    linewidth=1.0, capsize=2.5,
                )
 
            # ========================================================== #
            #  Upper panel (counts) — cosmetics                          #
            # ========================================================== #
            ax_top.set_xscale("log")
            ax_top.set_xlim(pt_bins[0], pt_bins[-1])
            ax_top.set_ylim(0.0, max(1.0, 1.23 * max_count))
            ax_top.tick_params(axis="both", labelsize=9)
            ax_top.apply_aspect()
            ax_top.xaxis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f"{v:g}")
            )
            ax_top.xaxis.set_minor_formatter(plt.NullFormatter())
 
            # flavour label centred in the upper panel
            ax_top.text(
                0.06, 0.95,
                r"$\bf{Private\; work}$" + "\nCMS data/simulation",
                transform=ax_top.transAxes,
                fontsize=9,
                va="top", ha="left",
                fontweight="normal",
            )
 
            # CMS label at top of each upper panel
            ax_top.text(
                0.00, 1.00,
                f"{_CHANNEL_LABELS[channel]}: {process}, {_FLAVOR_LABELS.get(flavor.lower(), flavor)}, WP {wp}",
                transform=ax_top.transAxes,
                fontsize=9,
                va="bottom", ha="left",
            )
            # Right side: channel and era (on all panels)
            ax_top.text(
                1.00, 1.01,
                f"({era}, {center_of_mass} TeV)",
                transform=ax_top.transAxes,
                fontsize=9,
                va="bottom", ha="right",
            )
 
            ax_top.set_ylabel("Jet count", fontsize=10)
 
            ax_top.legend(
                fontsize=8,
                loc="upper right",
                bbox_to_anchor=(0.95, 0.95),
                borderaxespad=0.0,
                ncol=1,
                framealpha=0.85,
                edgecolor="black",
            )
 
            # ========================================================== #
            #  Lower panel (efficiency) — cosmetics                      #
            # ========================================================== #
            ax_bot.set_xlim(pt_bins[0], pt_bins[-1])
            if np.isfinite(eff_y_min) and np.isfinite(eff_y_max):
                eff_span = eff_y_max - eff_y_min
                eff_padding = max(0.02, 0.15 * eff_span)
                if eff_span <= 0.0:
                    eff_padding = 0.05
                ax_bot.set_ylim(eff_y_min - eff_padding, eff_y_max + eff_padding + 0.3)
            else:
                ax_bot.set_ylim(-0.05, 0.05)
            ax_bot.tick_params(axis="both", labelsize=9)
            ax_bot.xaxis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f"{v:g}")
            )
            ax_bot.xaxis.set_minor_formatter(plt.NullFormatter())
            ax_bot.set_xlabel(
                r"Jet $p_{\mathrm{T}}$ [GeV]", fontsize=10,
            )
 
            ax_bot.set_ylabel("Efficiency", fontsize=10)
 
            ax_bot.legend(
                fontsize=9,
                loc="upper right",
                bbox_to_anchor=(0.95, 0.95),
                borderaxespad=0.0,
                framealpha=0.85,
                edgecolor="black",
            )
            ax_bot.yaxis.set_major_formatter(
                plt.FuncFormatter(lambda v, _: f"{v:.2f}")
            )
 
            # If upper panel y-axis is in scientific notation, move it to the side
            ax_top_formatter = ax_top.yaxis.get_major_formatter()
            ax_top.yaxis.set_label_position("left")
            
            # Check if we need to adjust y-axis range to prevent overlap with title
            ax_top.margins(y=0.0)
 
        for ext in ("png", "pdf"):
            fname = os.path.join(
                plot_dir, f"btag_{era}_{channel}_wp{wp}_{process}.{ext}"
            )
            fig.savefig(fname, bbox_inches="tight", dpi=150)
            log.info(f"Saved plot: {fname}")
 
        plt.close(fig)


# ---------------------------------------------------------------------------
# Per-process worker (module-level so it is picklable by ProcessPoolExecutor)
# ---------------------------------------------------------------------------

def _process_one(
    task: Tuple,
) -> Tuple[str, List[float], List[float], Dict, Dict, Dict]:
    """Process one MC sample; designed to run in a subprocess.

    Handles file discovery, RDataFrame histogram filling, efficiency
    conversion, and plotting, then returns plain-Python results that can
    be pickled back to the parent process. The :class:`WeightedEfficiency`
    accumulators are plain numpy, so they pickle back cleanly across the spawn.

    Args:
        task: ``(process, proc_conf, config, output_path, n_threads, log_level)``

    Returns:
        ``(sample_type, pt_bins, eta_bins, efficiencies, uncertainties, accumulators)``

    Raises:
        RuntimeError: No ROOT files found for the process, a required input
                      branch is missing, the RDataFrame is empty, or the weight
                      column is missing. All are fatal (no partial payload).
    """
    process, proc_conf, config, output_path, n_threads, log_level = task

    # Configure ROOT implicit MT for this subprocess.
    if n_threads == 0:
        ROOT.EnableImplicitMT()
    else:
        ROOT.EnableImplicitMT(n_threads)
    ROOT.gROOT.SetBatch(True)

    # Set up a stderr logger for this worker (the main-process file handler is
    # not inherited by spawned subprocesses).  Always update the formatter so
    # that reused worker processes (ProcessPoolExecutor recycles OS processes)
    # show the correct process name rather than the one from the first task.
    log = logging.getLogger("btag_efficiency")
    if not log.handlers:
        log.addHandler(logging.StreamHandler())
    log.handlers[0].setFormatter(
        logging.Formatter(f"[%(levelname)s][{process}] %(message)s")
    )
    log.setLevel(log_level)

    sample_type: str = proc_conf["sample_type"]
    log.info(f"Processing process '{process}' -> sample type '{sample_type}'")

    pt_bins, eta_bins = get_process_bins(config, process)
    log.info(f"  pt bins:  {pt_bins}")
    log.info(f"  eta bins: {eta_bins}")

    files = get_process_files(config, process)
    if not files:
        raise RuntimeError(
            f"No ROOT files found for process '{process}' under "
            f"'{config['file_path']}'. A partial payload must never be produced, "
            f"so a missing process is fatal."
        )
    log.info(f"  Found {len(files)} file(s):")
    for f in files:
        log.info(f"    {f}")

    # Fail-fast branch validation before any graph is built.
    validate_process_branches(files, config)

    chain = ROOT.TChain(config["tree"])
    for f in files:
        chain.Add(f)

    histograms, accumulators = calculate_efficiency_histograms(
        chain, config, process, pt_bins, eta_bins
    )
    efficiencies, uncertainties = accumulators_to_efficiencies(accumulators)

    log.info(f"  Plotting histograms and efficiencies for '{process}' ...")
    plot_histograms_and_efficiencies(
        histograms, efficiencies, uncertainties,
        pt_bins, eta_bins, process, config, output_path,
    )
    log.info(f"  Finished processing '{process}' (sample type key: '{sample_type}')")
    return sample_type, pt_bins, eta_bins, efficiencies, uncertainties, accumulators

def build_correctionlib_json(
    all_efficiencies: Dict[str, Dict[str, Dict[str, List[List[float]]]]],
    all_flavor_bins: Dict[str, Dict[str, Tuple[List[float], List[float]]]],
    config: Dict,
    output_path: str,
) -> cs.CorrectionSet:
    """Assemble the correctionlib CorrectionSet and write it to disk.

    ``all_flavor_bins[sample_type][flavor] = (pt_bins, eta_bins)`` carries the
    post-merge binning frozen per ``(sample_type, flavor)`` (all working points
    of a flavour share it).
    """
    wps: Dict[str, float] = config["btag_working_points"]
    flavors: Dict[str, int] = config["jet_flavor_categories"]

    sample_items = [
        cs.CategoryItem(
            key=process,
            value=_build_wp_category(
                wps, flavors, all_flavor_bins[process], proc_eff
            ),
        )
        for process, proc_eff in all_efficiencies.items()
    ]

    pt_bins_global: List[float] = config["jet_pt_bins"]
    eta_bins_global: List[float] = config["jet_eta_bins"]

    correction = cs.Correction(
        name="btag_efficiency",
        description=(
            "B-tagging efficiency as a function of sample type, working point, "
            "jet flavour, jet pT [GeV] and jet |eta|. "
            f"Era: {config.get('era', 'unknown')}. "
            "Efficiency = N(jets passing WP) / N(total jets) per bin. "
            "Binning may differ per sample type."
        ),
        version=1,
        inputs=[
            cs.Variable(name="sample_type",   type="string", description="MC sample type / process name."),
            cs.Variable(name="working_point", type="string", description="B-tag working point label (e.g. L, M, T, XT, XXT)."),
            cs.Variable(name="jet_flavor",    type="int",    description="Jet hadron flavour (b=5, c=4, light/udsg=0)."),
            cs.Variable(name="jet_eta",       type="real",   description=f"Jet |eta|; range ({min(eta_bins_global)}, {max(eta_bins_global)})."),
            cs.Variable(name="jet_pt",        type="real",   description=f"Jet pT [GeV]; range ({min(pt_bins_global)}, {max(pt_bins_global)}) GeV."),
        ],
        output=cs.Variable(name="efficiency", type="real", description="B-tagging efficiency."),
        data=cs.Category(
            nodetype="category",
            input="sample_type",
            content=sample_items,
        ),
    )

    cset = cs.CorrectionSet(
        schema_version=2,
        description="B-tagging efficiency corrections",
        corrections=[correction],
        compound_corrections=None,
    )

    os.makedirs(output_path, exist_ok=True)
    out_base = os.path.join(output_path, "btag_efficiency")
    write_json(f"{out_base}.json", cset)
    write_json(f"{out_base}.json.gz", cset)

    return cset


def _json_safe(array: np.ndarray) -> List:
    """Return ``array.tolist()`` with non-finite entries replaced by ``None``.

    Empty bins carry ``nan`` in the derived (efficiency/variance/effective-
    population) arrays; JSON has no ``NaN`` literal, so those become ``null``.
    """
    return np.where(np.isfinite(array), array, None).tolist()


def write_accumulators_sidecar(
    all_accumulators: Dict[str, Dict[str, Dict[str, WeightedEfficiency]]],
    all_bins: Dict[str, Tuple[List[float], List[float]]],
    config: Dict,
    output_path: str,
) -> str:
    """Write the diagnostic accumulator sidecar next to the correctionlib JSON.

    The sidecar (``btag_efficiency_accumulators.json``) carries, per
    ``sample_type -> working_point -> jet_flavor``: the six raw accumulator
    arrays (``sumw_total``/``sumw2_total``/``sumw_pass``/``sumw2_pass``/
    ``raw_total``/``raw_pass``) plus the derived diagnostics -- weighted
    ``efficiency`` (the production value), its ``variance`` and
    ``uncertainty`` (``sqrt(var)``), the ``effective_population``, and the
    ``unweighted_efficiency`` / ``unweighted_uncertainty`` binomial diagnostic.
    All arrays are ``[i_pt][i_eta]`` plain lists; ``nan`` becomes ``null``.
    """
    payload: Dict = {
        "description": (
            "Diagnostic sidecar for btag_efficiency.json. Per sample_type / "
            "working_point / jet_flavor: weighted sumw/sumw2 accumulators, raw "
            "counts, and derived quantities. The 'efficiency' field is the "
            "weighted production value written to btag_efficiency.json; "
            "'unweighted_efficiency' is the binomial diagnostic. Arrays are "
            "indexed [pt_bin][eta_bin]."
        ),
        "era": config.get("era", "unknown"),
        "channel": get_channel_display(config),
        "working_points": dict(config["btag_working_points"]),
        "jet_flavor_categories": dict(config["jet_flavor_categories"]),
        "processes": {},
    }

    for sample_type, wp_accs in all_accumulators.items():
        pt_bins, eta_bins = all_bins[sample_type]
        proc_entry: Dict = {
            "pt_bins": list(pt_bins),
            "eta_bins": list(eta_bins),
            "working_points": {},
        }
        for wp_name, flavor_accs in wp_accs.items():
            proc_entry["working_points"][wp_name] = {}
            for flavor_name, acc in flavor_accs.items():
                eff, var = weighted_efficiency(acc)
                unc = np.sqrt(var)
                n_eff = effective_population(acc)
                unw_eff, unw_unc = unweighted_efficiency(acc)

                entry = acc.to_dict()
                entry.update(
                    {
                        "efficiency": _json_safe(eff),
                        "variance": _json_safe(var),
                        "uncertainty": _json_safe(unc),
                        "effective_population": _json_safe(n_eff),
                        "unweighted_efficiency": _json_safe(unw_eff),
                        "unweighted_uncertainty": _json_safe(unw_unc),
                    }
                )
                proc_entry["working_points"][wp_name][flavor_name] = entry
        payload["processes"][sample_type] = proc_entry

    os.makedirs(output_path, exist_ok=True)
    sidecar_path = os.path.join(output_path, "btag_efficiency_accumulators.json")
    with open(sidecar_path, "w") as fout:
        json.dump(payload, fout, indent=4, allow_nan=False)

    return sidecar_path


def run_for_channel(config: Dict, args: argparse.Namespace) -> None:
    """Execute the full efficiency workflow for one channel-specific config."""
    channel = get_channel_label(config)

    output_path = os.path.join(
        config["output_base"], config["workdir_name"], config["era"], f"{channel}",
    )
    func.check_path(output_path)

    log = logging_helper.setup_logging(
        output_file=os.path.join(output_path, "btag_efficiency.log"),
        logger=logging.getLogger("btag_efficiency"),
        level=logging_helper.LOG_LEVEL,
    )

    log.info("Starting b-tag efficiency calculation")
    log.info(f"Channel: {channel}")
    log.info(f"Output directory: {output_path}")

    all_efficiencies: Dict[str, Dict[str, Dict[str, List[List[float]]]]] = {}
    all_uncertainties: Dict[str, Dict[str, Dict[str, List[List[float]]]]] = {}
    all_accumulators: Dict[str, Dict[str, Dict[str, WeightedEfficiency]]] = {}
    all_bins: Dict[str, Tuple[List[float], List[float]]] = {}

    if args.workers == 1:
        # -------------------------------------------------------------------
        # Single-process path: enable ROOT implicit MT in the main process.
        # -------------------------------------------------------------------
        if args.threads == 0:
            ROOT.EnableImplicitMT()
        else:
            ROOT.EnableImplicitMT(args.threads)

        for process, proc_conf in config["processes"].items():
            sample_type: str = proc_conf["sample_type"]
            log.info(f"Processing process '{process}' -> sample type '{sample_type}'")

            pt_bins, eta_bins = get_process_bins(config, process)
            log.info(f"  pt bins:  {pt_bins}")
            log.info(f"  eta bins: {eta_bins}")

            files = get_process_files(config, process)
            if not files:
                raise RuntimeError(
                    f"No ROOT files found for process '{process}' under "
                    f"'{config['file_path']}'. A partial payload must never be "
                    f"produced, so a missing process is fatal."
                )

            log.info(f"  Found {len(files)} file(s):")
            for f in files:
                log.info(f"    {f}")

            # Fail-fast branch validation before any graph is built.
            validate_process_branches(files, config)

            chain = ROOT.TChain(config["tree"])
            for f in files:
                chain.Add(f)

            # A RuntimeError here (empty frame, missing weight column, probe
            # vector violation) is fatal -- no silent skip.
            histograms, accumulators = calculate_efficiency_histograms(
                chain, config, process, pt_bins, eta_bins
            )

            efficiencies, uncertainties = accumulators_to_efficiencies(accumulators)
            log.info(f"  Plotting histograms and efficiencies for '{process}' ...")
            plot_histograms_and_efficiencies(
                histograms, efficiencies, uncertainties,
                pt_bins, eta_bins, process, config, output_path,
            )
            all_efficiencies[sample_type] = efficiencies
            all_uncertainties[sample_type] = uncertainties
            all_accumulators[sample_type] = accumulators
            all_bins[sample_type] = (pt_bins, eta_bins)
            log.info(f"  Finished processing '{process}' (sample type key: '{sample_type}')")

    else:
        # -------------------------------------------------------------------
        # Multi-process path: one subprocess per MC sample.
        # Each worker enables ROOT implicit MT with --threads threads.
        # -------------------------------------------------------------------
        log.info(
            f"Using {args.workers} parallel worker(s) "
            f"(ROOT threads per worker: {'all' if args.threads == 0 else args.threads})"
        )
        tasks = [
            (process, proc_conf, config, output_path, args.threads, logging_helper.LOG_LEVEL)
            for process, proc_conf in config["processes"].items()
        ]
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {executor.submit(_process_one, task): task[0] for task in tasks}
            for future in concurrent.futures.as_completed(futures):
                process = futures[future]
                # A worker failure is fatal: re-raise so the whole run aborts
                # rather than silently dropping a process (no partial payload).
                (
                    sample_type,
                    pt_bins,
                    eta_bins,
                    efficiencies,
                    uncertainties,
                    accumulators,
                ) = future.result()
                log.info(f"Finished '{process}' (sample type: '{sample_type}')")
                all_efficiencies[sample_type] = efficiencies
                all_uncertainties[sample_type] = uncertainties
                all_accumulators[sample_type] = accumulators
                all_bins[sample_type] = (pt_bins, eta_bins)

    if not all_accumulators:
        # No processes were configured at all; a run that produces nothing is a
        # configuration error, not a success.
        raise RuntimeError(
            "No processes were configured for this channel; refusing to write an "
            "empty payload."
        )

    # -----------------------------------------------------------------------
    # Diagnostic accumulator sidecar (pre-merge, finest binning).
    # -----------------------------------------------------------------------
    log.info("Writing diagnostic accumulator sidecar ...")
    sidecar_path = write_accumulators_sidecar(
        all_accumulators, all_bins, config, output_path
    )
    log.info(f"Accumulator sidecar written to '{sidecar_path}'")

    # -----------------------------------------------------------------------
    # Strict validation: raw nesting (fatal), deterministic per-(process,
    # flavor) bin merging, and correctionlib round-trip. Produces the merged
    # payload plus the machine-readable validation report.
    # -----------------------------------------------------------------------
    provenance = validate_and_build_payload(
        all_accumulators, all_bins, config, output_path, channel, log
    )

    log.info("B-tag efficiency calculation finished.")
    return provenance


def validate_and_build_payload(
    all_accumulators: Dict[str, Dict[str, Dict[str, WeightedEfficiency]]],
    all_bins: Dict[str, Tuple[List[float], List[float]]],
    config: Dict,
    output_path: str,
    channel: str,
    log: logging.Logger,
) -> Dict:
    """Run the strict gates, merge bins, build the payload and the report.

    Raw pass-set nesting violations are fatal (:func:`bval.check_wp_nesting`).
    Weighted quality-gate violations are resolved by deterministic per-(process,
    flavor) bin merging (:func:`bval.merge_bins`), applied to all working points
    of the category together. The produced correctionlib JSON is round-tripped
    against the accumulator-derived efficiencies before the report is finalized;
    any mismatch fails the run. Returns a per-channel provenance dict.
    """
    wps: Dict[str, float] = config["btag_working_points"]
    flavors: Dict[str, int] = config["jet_flavor_categories"]
    thresholds = bval.get_validation_thresholds(config)
    log.info(f"Validation thresholds: {thresholds}")

    merged_efficiencies: Dict[str, Dict[str, Dict[str, List[List[float]]]]] = {}
    all_flavor_bins: Dict[str, Dict[str, Tuple[List[float], List[float]]]] = {}
    categories: Dict[str, Dict[str, Dict]] = {}

    for sample_type, wp_accs in all_accumulators.items():
        pt_bins, eta_bins = all_bins[sample_type]
        merged_efficiencies[sample_type] = {wp: {} for wp in wps}
        all_flavor_bins[sample_type] = {}
        categories[sample_type] = {}

        for flavor_name in flavors:
            acc_by_wp = {wp: wp_accs[wp][flavor_name] for wp in wps}

            # Raw exact-subset nesting: any violation is immediately fatal.
            nesting = bval.check_wp_nesting(acc_by_wp)
            if not nesting.raw_ok:
                raise RuntimeError(
                    f"Raw pass-set nesting violated for sample_type "
                    f"'{sample_type}', flavor '{flavor_name}' "
                    f"(L>=M>=T>=XT>=XXT must hold per bin); this indicates a "
                    f"selection/discriminator bug. Violations: "
                    f"{nesting.raw_violations}"
                )

            pre = bval.BinnedAccumulators(
                accumulators=acc_by_wp, pt_bins=pt_bins, eta_bins=eta_bins
            )
            merged, history = bval.merge_bins(pre, thresholds)
            if history:
                log.info(
                    f"  Merged bins for '{sample_type}'/'{flavor_name}': "
                    f"{len(history)} op(s); pt {pre.pt_bins} -> {merged.pt_bins}, "
                    f"eta {pre.eta_bins} -> {merged.eta_bins}"
                )

            all_flavor_bins[sample_type][flavor_name] = (
                list(merged.pt_bins),
                list(merged.eta_bins),
            )
            for wp in wps:
                eff, _var = weighted_efficiency(merged.accumulators[wp])
                eff = np.where(np.isfinite(eff), eff, 0.0)
                merged_efficiencies[sample_type][wp][flavor_name] = eff.tolist()

            categories[sample_type][flavor_name] = bval.summarize_category(
                pre, merged, history, thresholds, nesting
            )

    log.info("Building correctionlib JSON ...")
    build_correctionlib_json(merged_efficiencies, all_flavor_bins, config, output_path)
    out_base = os.path.join(output_path, "btag_efficiency")
    log.info(f"Correctionlib files written to '{out_base}.json[.gz]'")

    # Round-trip gate: the produced JSON must reproduce the stored efficiencies
    # exactly at every bin centre and (in-range) boundary.
    log.info("Running correctionlib round-trip gate ...")
    rt_ok, rt_n, rt_mismatches = bval.roundtrip_check(
        f"{out_base}.json", merged_efficiencies, all_flavor_bins, wps, flavors
    )
    log.info(f"  Round-trip evaluated {rt_n} lookup(s); ok={rt_ok}")

    status = "passed" if rt_ok else "failed"
    report_path = os.path.join(output_path, "btag_validation_report.json")
    bval.write_validation_report(
        report_path,
        status=status,
        thresholds=thresholds,
        merge_version=bval.MERGE_VERSION,
        categories=categories,
        roundtrip={
            "ok": rt_ok,
            "n_checked": rt_n,
            "mismatches": rt_mismatches[:50],
        },
        era=config.get("era", ""),
        channel=get_channel_display(config),
        channels_present=[channel],
        required_channels=_required_channels(config),
    )
    log.info(f"Validation report written to '{report_path}' (status: {status})")

    if not rt_ok:
        raise RuntimeError(
            f"correctionlib round-trip gate failed with {len(rt_mismatches)} "
            f"mismatch(es); refusing to publish. See '{report_path}'."
        )

    payload_files = [
        f"{out_base}.json",
        f"{out_base}.json.gz",
        sidecar_or_none(output_path),
        report_path,
    ]
    payload_files = [p for p in payload_files if p and os.path.isfile(p)]
    checksums = {
        os.path.basename(p): bval._sha256(p) for p in payload_files
    }

    return {
        "channel": channel,
        "status": status,
        "output_path": output_path,
        "report_path": report_path,
        "payload_files": payload_files,
        "checksums": checksums,
        "required_channels": _required_channels(config),
    }


def sidecar_or_none(output_path: str) -> str:
    """Return the accumulator sidecar path if it exists, else empty string."""
    path = os.path.join(output_path, "btag_efficiency_accumulators.json")
    return path if os.path.isfile(path) else ""


def _required_channels(config: Dict) -> List[str]:
    """Channels that must all be present before a payload may be published.

    Configurable via the ``required_channels`` key; defaults to the three
    fake-factor analysis channels ``et``, ``mt``, ``tt``.
    """
    required = config.get("required_channels")
    if required:
        return [str(c) for c in required]
    return ["et", "mt", "tt"]


def publish_all_channels(provenances: List[Dict], config: Dict) -> None:
    """Two-phase, all-or-nothing atomic install of every channel's payload.

    A partial multi-channel publish must never happen: every channel's
    payload is staged and checksum-verified into a tmp dir (``bval.stage_
    payload``) *before* any ``target_dir`` is touched. Only once every single
    channel has staged cleanly are the ``os.replace`` swaps performed
    (``bval.commit_staged``) for all of them.

    If any channel fails to stage (e.g. a checksum mismatch), every tmp dir
    staged so far -- including channels that staged fine before the failing
    one -- is cleaned up, nothing is published, and the exception is
    re-raised (aborting the whole run with a non-zero exit code under
    ``__main__``). This replaces a previous channel-by-channel loop that
    caught and merely logged a ``RuntimeError`` per channel, which could
    leave earlier channels published while a later one silently failed to
    install, and still exit 0.
    """
    log = logging.getLogger("btag_efficiency")
    channels_present = sorted({prov["channel"] for prov in provenances})
    required = _required_channels(config)

    staged: List[Tuple[str, str, str]] = []  # (tmp_dir, target_dir, channel)
    try:
        for prov in provenances:
            install_provenance = dict(prov)
            install_provenance["channels_present"] = channels_present
            install_provenance["required_channels"] = required
            target_dir = os.path.join(
                config["output_base"],
                config["workdir_name"],
                config["era"],
                "published",
                prov["channel"],
            )
            tmp_dir = bval.stage_payload(
                prov["payload_files"], install_provenance, target_dir
            )
            staged.append((tmp_dir, target_dir, prov["channel"]))
            log.info(f"Staged payload for channel '{prov['channel']}' -> '{tmp_dir}'")
    except Exception as exc:
        # Any channel failing to stage aborts the ENTIRE publish: clean up
        # every tmp dir staged so far and leave every target_dir untouched.
        for tmp_dir, _target_dir, channel in staged:
            shutil.rmtree(tmp_dir, ignore_errors=True)
            log.warning(f"Discarded staged payload for channel '{channel}'.")
        log.error(
            f"Not publishing any channel: at least one channel failed to "
            f"stage ({exc}). Nothing was published this run."
        )
        raise

    # Every channel staged and checksum-verified cleanly -- commit all of them.
    for tmp_dir, target_dir, channel in staged:
        bval.commit_staged(tmp_dir, target_dir)
        log.info(f"Published payload for channel '{channel}' to '{target_dir}'")


# ---------------------------------------------------------------------------
# Provenance manifest writer
# ---------------------------------------------------------------------------
# Marker string every consumer of the installed payload recognises as the
# validated production chain that produced it.
PRODUCED_BY = "sm_btag_efficiency_config -> TauFakeFactors -> install"


def _git_commit(repo_path: Optional[str]) -> Optional[str]:
    """Return the ``git rev-parse HEAD`` of ``repo_path`` or ``None``.

    Never raises: a missing path, a non-git directory, or a missing ``git``
    binary all yield ``None`` so provenance writing degrades gracefully.
    """
    if not repo_path or not os.path.isdir(repo_path):
        return None
    try:
        result = subprocess.run(
            ["git", "-C", repo_path, "rev-parse", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        )
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        return None
    commit = result.stdout.strip()
    return commit or None


def _read_digest_file(path: Optional[str]) -> Optional[Dict]:
    """Read the Task-14 ``--emit-digest`` JSON (nicks + per-sample filelist digests).

    Returns the parsed digest dict, or ``None`` if no ``digest_file`` is
    configured or the file is missing/unreadable.
    """
    if not path or not os.path.isfile(path):
        return None
    try:
        with open(path, "r") as handle:
            return json.load(handle)
    except (OSError, ValueError):
        return None


def _read_selection_contract(config_dir: Optional[str]) -> Optional[Dict]:
    """Return ``{path, sha256}`` for the selection contract, or ``None``.

    Reads ``selection_contract_2018_v1.yaml`` from the calculator config's
    directory and returns its stored ``contract_sha256`` (the value the bbtautau
    export script wrote); the parity test independently checks the preselection
    matches the contract.
    """
    if not config_dir:
        return None
    contract_path = os.path.join(config_dir, "selection_contract_2018_v1.yaml")
    if not os.path.isfile(contract_path):
        return None
    try:
        with open(contract_path, "r") as handle:
            contract = func.configured_yaml.load(handle)
    except Exception:
        return None
    return {
        "path": contract_path,
        "contract_version": contract.get("contract_version"),
        "sha256": contract.get("contract_sha256"),
    }


def _composite_weight_expression(mc_weights: Dict) -> str:
    """Build a human-readable composite weight expression from an mc_weights block.

    Mirrors preselection.py: the per-event ``weight`` is the product of every
    mc_weights term. Terms whose value is an empty string are computed
    internally (gen normalization, lumi) and are represented by their key name.
    """
    terms: List[str] = []
    for key, value in mc_weights.items():
        text = str(value).strip()
        if text == "":
            terms.append(f"<{key}>")
        else:
            terms.append(f"({text})")
    return " * ".join(terms)


def _channel_weight_expressions(config_dir: Optional[str], channels: List[str]) -> Dict[str, str]:
    """Return ``{channel: composite weight expression}`` from the preselection configs."""
    expressions: Dict[str, str] = {}
    if not config_dir:
        return expressions
    for channel in channels:
        presel_path = os.path.join(config_dir, f"preselection_{channel}.yaml")
        if not os.path.isfile(presel_path):
            continue
        try:
            with open(presel_path, "r") as handle:
                presel = func.configured_yaml.load(handle)
        except Exception:
            continue
        mc_weights = presel.get("mc_weights") or {}
        expressions[channel] = _composite_weight_expression(mc_weights)
    return expressions


def _merge_history_from_provenances(channel_provenances: List[Dict]) -> Dict[str, Dict]:
    """Extract binning + merge history per channel from the validation reports."""
    history: Dict[str, Dict] = {}
    for prov in channel_provenances or []:
        report_path = prov.get("report_path")
        channel = prov.get("channel")
        if not report_path or not os.path.isfile(report_path):
            continue
        try:
            with open(report_path, "r") as handle:
                report = json.load(handle)
        except (OSError, ValueError):
            continue
        per_category: Dict[str, Dict] = {}
        for sample_type, cats in (report.get("categories") or {}).items():
            per_category[sample_type] = {
                flavor: {
                    "pt_bins_pre": entry.get("pt_bins_pre"),
                    "pt_bins_post": entry.get("pt_bins_post"),
                    "eta_bins_pre": entry.get("eta_bins_pre"),
                    "eta_bins_post": entry.get("eta_bins_post"),
                    "merge_history": entry.get("merge_history"),
                }
                for flavor, entry in cats.items()
            }
        history[channel] = per_category
    return history


def build_provenance_meta(
    config: Dict,
    *,
    config_dir: Optional[str] = None,
    validation_status: str = "passed",
    channels: Optional[List[str]] = None,
    channel_provenances: Optional[List[Dict]] = None,
    working_points: Optional[Dict[str, float]] = None,
    bbtautau_repo: Optional[str] = None,
) -> Dict:
    """Assemble the provenance metadata dict (everything except the manifest).

    The returned dict carries the fields the bbtautau validated-payload gate
    (``btag_payloads.require_validated_payload``) reads -- ``validation_status``
    and ``produced_by`` -- plus the full production provenance: repo commits
    (TauFakeFactors / bbtautau / sample database), production tag, sample nick +
    filelist digests (Task-14 digest file), selection-contract checksum, the
    per-channel composite weight expression, the binning + merge history, the
    pinned BTV working-point payload path + sha256 and the resolved working
    points. ``install_provenance_payload`` fills in the ``manifest`` (per-scope
    payload sha256s) from the actually-installed files so it cannot drift.
    """
    if channels is None:
        channels = get_config_channels(config)

    btv_source = config.get("btag_working_points_source")
    btv_payload = None
    if btv_source:
        btv_payload = {
            "path": btv_source,
            "sha256": bval._sha256(btv_source) if os.path.isfile(btv_source) else None,
        }

    return {
        # -- fields the bbtautau require_validated_payload gate reads ----------
        "validation_status": validation_status,
        "produced_by": PRODUCED_BY,
        # -- full production provenance ---------------------------------------
        "era": config.get("era"),
        "production_tag": config.get("production_tag"),
        "repo_commits": {
            "TauFakeFactors": _git_commit(func.TAU_FAKE_FACTORS_DIR),
            "bbtautau": _git_commit(bbtautau_repo or config.get("bbtautau_repo")),
            "sample_database": _git_commit(config.get("sample_database")),
        },
        "selection_contract": _read_selection_contract(config_dir),
        "sample_digest": _read_digest_file(config.get("digest_file")),
        "weight_expression": _channel_weight_expressions(config_dir, channels),
        "binning": {
            "jet_pt_bins": list(config.get("jet_pt_bins", [])),
            "jet_eta_bins": list(config.get("jet_eta_bins", [])),
        },
        "merge_history": _merge_history_from_provenances(channel_provenances or []),
        "btv_working_points_payload": btv_payload,
        "working_points": dict(
            working_points or config.get("btag_working_points") or {}
        ),
    }


def install_provenance_payload(
    scope_payload_files: Dict[str, str],
    provenance_meta: Dict,
    target_dir: str,
) -> str:
    """Atomically install per-scope payloads + ``provenance.json`` into ``target_dir``.

    ``scope_payload_files`` maps each channel/scope to the correctionlib
    ``.json.gz`` produced for it. Each is copied into a fresh staging directory
    as ``btag_efficiency_<scope>.json.gz`` (the exact name the bbtautau gate
    globs), its SHA256 is recorded into the ``manifest`` of a copy of
    ``provenance_meta``, and once every file is staged the whole directory is
    ``os.replace``-d onto ``target_dir`` (via :func:`bval._atomic_replace_dir`).
    The manifest is computed from the staged files themselves, so a
    ``provenance.json`` written by this function always matches the payloads it
    ships next to -- exactly what ``require_validated_payload`` re-checks.

    Returns ``target_dir`` on success; on any failure the staging directory is
    removed and the exception propagates (no partial install).
    """
    tmp_dir = f"{target_dir}.tmp-{os.getpid()}"
    shutil.rmtree(tmp_dir, ignore_errors=True)
    os.makedirs(tmp_dir)
    try:
        manifest: Dict[str, Dict[str, str]] = {}
        for scope, src in sorted(scope_payload_files.items()):
            file_name = f"btag_efficiency_{scope}.json.gz"
            dst = os.path.join(tmp_dir, file_name)
            shutil.copy2(src, dst)
            manifest[file_name] = {"sha256": bval._sha256(dst)}

        provenance = dict(provenance_meta)
        provenance["manifest"] = manifest

        with open(os.path.join(tmp_dir, "provenance.json"), "w") as handle:
            json.dump(provenance, handle, indent=2, sort_keys=True)
            handle.write("\n")

        bval._atomic_replace_dir(tmp_dir, target_dir)
    except Exception:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise
    return target_dir


def _scope_payload_files(channel_provenances: List[Dict]) -> Dict[str, str]:
    """Map each channel to its correctionlib ``btag_efficiency.json.gz`` output."""
    scope_files: Dict[str, str] = {}
    for prov in channel_provenances:
        channel = prov["channel"]
        for payload_file in prov.get("payload_files", []):
            if os.path.basename(payload_file) == "btag_efficiency.json.gz":
                scope_files[channel] = payload_file
                break
    return scope_files


def write_provenance_payload(
    provenances: List[Dict],
    config: Dict,
    *,
    config_dir: Optional[str] = None,
    working_points: Optional[Dict[str, float]] = None,
) -> Optional[str]:
    """Assemble + install the validated-payload provenance directory.

    Gathers each channel's ``btag_efficiency.json.gz``, builds the provenance
    metadata (:func:`build_provenance_meta`) and installs the gate-compatible
    payload directory (per-scope ``btag_efficiency_<scope>.json.gz`` +
    ``provenance.json``) under ``<output_base>/<workdir_name>/<era>/payload``.
    Returns the installed directory, or ``None`` if no per-scope payloads were
    found.
    """
    channels = sorted({prov["channel"] for prov in provenances})
    scope_files = _scope_payload_files(provenances)
    if not scope_files:
        return None

    all_passed = all(prov.get("status") == "passed" for prov in provenances)
    meta = build_provenance_meta(
        config,
        config_dir=config_dir,
        validation_status="passed" if all_passed else "failed",
        channels=channels,
        channel_provenances=provenances,
        working_points=working_points,
    )
    target_dir = os.path.join(
        config["output_base"], config["workdir_name"], config["era"], "payload"
    )
    return install_provenance_payload(scope_files, meta, target_dir)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    args = parser.parse_args()

    if args.workers > 1:
        # Use 'spawn' so each worker gets a clean ROOT instance (safe on Linux).
        multiprocessing.set_start_method("spawn")

    logging_helper.LOG_LEVEL = getattr(logging, args.log_level.upper(), logging.INFO)

    config = func.load_config(args.config_file)

    # resolve path-like settings with CLI > env > config (> default)
    # precedence, once, right after loading the config, and write them back
    # into the config dict so that all downstream code (get_process_files(),
    # run_for_channel(), etc.) keeps reading config["file_path"] /
    # config["output_base"] unchanged. get_channel_configs() deep-copies this
    # config per channel, so the resolved values propagate to every channel.
    config["file_path"] = func.resolve_path_setting(
        config=config, key="file_path", cli_value=args.file_path, env_var="TFF_FILE_PATH"
    )
    config["output_base"] = func.resolve_path_setting(
        config=config,
        key="output_base",
        cli_value=args.output_path,
        env_var="TFF_OUTPUT_PATH",
        default="workdir",
    )

    # Validate the configured b-tag working points against the pinned BTV
    # payload (fatal on any disagreement) and write the validated set back into
    # the config so every channel run uses exactly the pinned cut values.
    resolved_working_points = resolve_working_points(config)
    if resolved_working_points is not None:
        config["btag_working_points"] = resolved_working_points

    config_dir = os.path.dirname(os.path.abspath(args.config_file))

    # Run every channel first (any fatal gate aborts the whole run before
    # anything is published -- a partial payload must never be produced).
    provenances = [
        run_for_channel(channel_config, args)
        for channel_config in get_channel_configs(config)
    ]

    # Publish atomically only once every channel has passed and all required
    # channels are present. channels_present is the union across channel runs.
    # This is two-phase and all-or-nothing across channels: a mid-loop
    # failure (e.g. a corrupted checksum on a later channel) must never leave
    # only some channels published (see publish_all_channels docstring).
    publish_all_channels(provenances, config)

    # Assemble + install the validated-payload provenance directory (per-scope
    # btag_efficiency_<scope>.json.gz + provenance.json) consumed by the
    # bbtautau require_validated_payload gate.
    payload_dir = write_provenance_payload(
        provenances,
        config,
        config_dir=config_dir,
        working_points=config.get("btag_working_points"),
    )
    if payload_dir:
        logging.getLogger("btag_efficiency").info(
            f"Installed validated b-tag efficiency payload + provenance to '{payload_dir}'"
        )
