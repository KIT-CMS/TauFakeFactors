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
import logging
import multiprocessing
import os
from typing import Dict, List, Tuple

import ROOT
import matplotlib.pyplot as plt
import mplhep as hep
import numpy as np
import correctionlib.schemav2 as cs

import CustomLogging as logging_helper
import helper.functions as func
from helper.correctionlib_json import write_json

hep.style.use(hep.style.CMS)
plt.rcParams["axes.linewidth"] = 1.0 # set non bold axes lines


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description="Calculate b-tagging efficiency and write a correctionlib JSON."
)
parser.add_argument(
    "--config-file",
    required=True,
    help="Path to the btag efficiency config YAML file.",
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
# Efficiency calculation
# ---------------------------------------------------------------------------
def calculate_efficiency_histograms(
    chain: ROOT.TChain,
    config: Dict,
    process: str,
    pt_bins: List[float],
    eta_bins: List[float],
) -> Dict[str, Dict[str, Tuple[ROOT.TH2D, ROOT.TH2D]]]:
    """Fill 2-D (pT, |eta|) histograms for every (WP, flavour) combination.

    For each combination the function returns a pair of histograms:
    ``(h_pass, h_total)`` where *h_pass* counts jets that pass the btag
    discriminator threshold and *h_total* counts all jets of the given
    flavour.

    Args:
        chain:   ROOT TChain with all files for this process already attached.
        config:  Loaded configuration dictionary.
        process: Process name (used only for histogram titles / logging).

    Returns:
        ``results[wp_name][flavor_name] = (h_pass, h_total)``
    """
    log = logging.getLogger("btag_efficiency")

    jet_pt_col = config["jet_pt_column"]
    jet_eta_col = config["jet_eta_column"]
    jet_flavor_col = config["jet_flavor_column"]
    jet_btag_col = config["jet_btag_column"]
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

    results: Dict[str, Dict[str, Tuple]] = {}

    for wp_name, wp_cut in wps.items():
        results[wp_name] = {}
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

            rdf_loc = (
                rdf.Define(col_pt_all, f"{jet_pt_col_use}[{flavor_mask}]")
                   .Define(col_eta_all, f"abs({jet_eta_col_use}[{flavor_mask}])")
                   .Define(col_pt_pass, f"{jet_pt_col_use}[{pass_mask}]")
                   .Define(col_eta_pass, f"abs({jet_eta_col_use}[{pass_mask}])")
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

            # Trigger the event loop and clone the histograms so they outlive
            # the lazy RDataFrame action.
            results[wp_name][flavor_name] = (
                h_pass.GetValue().Clone(),
                h_total.GetValue().Clone(),
            )

    return results


def histograms_to_efficiencies(
    histograms: Dict[str, Dict[str, Tuple[ROOT.TH2D, ROOT.TH2D]]],
    pt_bins: List[float],
    eta_bins: List[float],
) -> Tuple[Dict[str, Dict[str, List[List[float]]]], Dict[str, Dict[str, List[List[float]]]]]:
    """Convert (h_pass, h_total) histogram pairs to efficiency grids.

    Args:
        histograms: ``histograms[wp][flavor] = (h_pass, h_total)``
        pt_bins:    pT bin edges.
        eta_bins:   |eta| bin edges.

    Returns:
        Tuple of two dicts with identical structure ``[wp][flavor]``:
        - efficiencies: efficiency value per (pT, |eta|) bin.
        - stat_uncertainties: binomial statistical uncertainty per bin
          (sigma = sqrt(eff * (1 - eff) / total)); zero when total == 0.
    """
    log = logging.getLogger("btag_efficiency")

    n_pt = len(pt_bins) - 1
    n_eta = len(eta_bins) - 1

    efficiencies: Dict[str, Dict[str, List[List[float]]]] = {}
    stat_uncertainties: Dict[str, Dict[str, List[List[float]]]] = {}
    for wp_name, flavor_hists in histograms.items():
        efficiencies[wp_name] = {}
        stat_uncertainties[wp_name] = {}
        for flavor_name, (h_pass, h_total) in flavor_hists.items():
            eff_grid: List[List[float]] = []
            unc_grid: List[List[float]] = []
            for i_pt in range(1, n_pt + 1):
                eta_row: List[float] = []
                unc_row: List[float] = []
                for i_eta in range(1, n_eta + 1):
                    total = h_total.GetBinContent(i_pt, i_eta)
                    passed = h_pass.GetBinContent(i_pt, i_eta)
                    if total > 0.0:
                        eff = min(passed / total, 1.0)
                        unc = float(np.sqrt(eff * (1.0 - eff) / total))
                    else:
                        log.warning(
                            f"Zero denominator in bin (pT bin {i_pt}, eta bin {i_eta}) "
                            f"for WP={wp_name}, flavor={flavor_name}. Setting eff=0."
                        )
                        eff = 0.0
                        unc = 0.0
                    eta_row.append(float(eff))
                    unc_row.append(unc)
                eff_grid.append(eta_row)
                unc_grid.append(unc_row)
            efficiencies[wp_name][flavor_name] = eff_grid
            stat_uncertainties[wp_name][flavor_name] = unc_grid

    return efficiencies, stat_uncertainties


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
        flow="clamp",
    )


def _build_flavor_category(
    flavors: Dict[str, int],
    pt_bins: List[float],
    eta_bins: List[float],
    wp_efficiencies: Dict[str, List[List[float]]],
) -> cs.Category:
    """Category node over jet flavour labels."""
    return cs.Category(
        nodetype="category",
        input="jet_flavor",
        content=[
            cs.CategoryItem(
                key=flavors[flavor_name],
                value=_build_eta_binning(eta_bins, pt_bins, wp_efficiencies[flavor_name]),
            )
            for flavor_name in flavors
        ],
    )


def _build_wp_category(
    wps: Dict[str, float],
    flavors: Dict[str, int],
    pt_bins: List[float],
    eta_bins: List[float],
    sample_efficiencies: Dict[str, Dict[str, List[List[float]]]],
) -> cs.Category:
    """Category node over b-tag working points."""
    return cs.Category(
        nodetype="category",
        input="working_point",
        content=[
            cs.CategoryItem(
                key=wp_name,
                value=_build_flavor_category(
                    flavors, pt_bins, eta_bins, sample_efficiencies[wp_name]
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
                f"({era}, 13.6 TeV)",
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
) -> Tuple[str, List[float], List[float], Dict, Dict]:
    """Process one MC sample; designed to run in a subprocess.

    Handles file discovery, RDataFrame histogram filling, efficiency
    conversion, and plotting, then returns plain-Python results that can
    be pickled back to the parent process.

    Args:
        task: ``(process, proc_conf, config, output_path, n_threads, log_level)``

    Returns:
        ``(sample_type, pt_bins, eta_bins, efficiencies, stat_uncertainties)``

    Raises:
        FileNotFoundError: No ROOT files found for the process.
        RuntimeError:      The RDataFrame for the process is empty.
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
        raise FileNotFoundError(
            f"No ROOT files found for process '{process}' under "
            f"'{config['file_path']}'."
        )
    log.info(f"  Found {len(files)} file(s):")
    for f in files:
        log.info(f"    {f}")

    chain = ROOT.TChain(config["tree"])
    for f in files:
        chain.Add(f)

    histograms = calculate_efficiency_histograms(chain, config, process, pt_bins, eta_bins)
    efficiencies, stat_uncertainties = histograms_to_efficiencies(histograms, pt_bins, eta_bins)

    log.info(f"  Plotting histograms and efficiencies for '{process}' ...")
    plot_histograms_and_efficiencies(
        histograms, efficiencies, stat_uncertainties,
        pt_bins, eta_bins, process, config, output_path,
    )
    log.info(f"  Finished processing '{process}' (sample type key: '{sample_type}')")
    return sample_type, pt_bins, eta_bins, efficiencies, stat_uncertainties

def build_correctionlib_json(
    all_efficiencies: Dict[str, Dict[str, Dict[str, List[List[float]]]]],
    all_bins: Dict[str, Tuple[List[float], List[float]]],
    config: Dict,
    output_path: str,
) -> cs.CorrectionSet:
    """Assemble the correctionlib CorrectionSet and write it to disk."""
    wps: Dict[str, float] = config["btag_working_points"]
    flavors: Dict[str, int] = config["jet_flavor_categories"]

    sample_items = [
        cs.CategoryItem(
            key=process,
            value=_build_wp_category(
                wps, flavors, all_bins[process][0], all_bins[process][1], proc_eff
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


def run_for_channel(config: Dict, args: argparse.Namespace) -> None:
    """Execute the full efficiency workflow for one channel-specific config."""
    channel = get_channel_label(config)

    output_path = os.path.join(
        "workdir", config["workdir_name"], config["era"], f"{channel}",
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
                log.warning(
                    f"No ROOT files found for process '{process}' under "
                    f"'{config['file_path']}'. Skipping."
                )
                continue

            log.info(f"  Found {len(files)} file(s):")
            for f in files:
                log.info(f"    {f}")

            chain = ROOT.TChain(config["tree"])
            for f in files:
                chain.Add(f)

            try:
                histograms = calculate_efficiency_histograms(chain, config, process, pt_bins, eta_bins)
            except RuntimeError as exc:
                log.warning(f"Skipping process '{process}': {exc}")
                continue

            efficiencies, stat_uncertainties = histograms_to_efficiencies(
                histograms,
                pt_bins,
                eta_bins,
            )
            log.info(f"  Plotting histograms and efficiencies for '{process}' ...")
            plot_histograms_and_efficiencies(
                histograms, efficiencies, stat_uncertainties,
                pt_bins, eta_bins, process, config, output_path,
            )
            all_efficiencies[sample_type] = efficiencies
            all_uncertainties[sample_type] = stat_uncertainties
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
                try:
                    sample_type, pt_bins, eta_bins, efficiencies, stat_uncertainties = future.result()
                except (FileNotFoundError, RuntimeError) as exc:
                    log.warning(f"Skipping process '{process}': {exc}")
                    continue
                log.info(f"Finished '{process}' (sample type: '{sample_type}')")
                all_efficiencies[sample_type] = efficiencies
                all_uncertainties[sample_type] = stat_uncertainties
                all_bins[sample_type] = (pt_bins, eta_bins)

    if not all_efficiencies:
        log.error("No processes were successfully processed. Aborting.")
        raise SystemExit(1)

    log.info("Building correctionlib JSON ...")
    build_correctionlib_json(all_efficiencies, all_bins, config, output_path)
    log.info(
        f"Correctionlib files written to '{output_path}/btag_efficiency.json[.gz]'"
    )
    log.info("B-tag efficiency calculation finished.")


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
    for channel_config in get_channel_configs(config):
        run_for_channel(channel_config, args)
