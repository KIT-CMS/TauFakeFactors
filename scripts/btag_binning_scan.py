#!/usr/bin/env python3
"""Apply the at-least-N-jets binning rule to a b-tag efficiency config.

The strict CROWN consumer (btag_sf_strict.cxx) throws on an efficiency of 0 or 1
and on a tie between two working points, so every cell of the (pT, |eta|) grid
must hold jets in each disjoint score region of the configured working points:
below the loosest WP, between each pair of adjacent WPs and above the tightest
WP (three regions for [T, M]). The binning rule of the efficiency configs is:
per pool and jet flavour, take the finest grid whose thinnest region still holds
at least --min-jets jets (default 5) in every channel the pool is measured in.

For one era directory (configs/btag_efficiency/<era>) the script counts the
probe jets of every pool in those regions, on the preselection outputs that
btag_efficiency.py reads ({file_path}/preselection/{era}/{channel}/{process}*.root,
unweighted, pT above the last edge not counted, as in the measurement), and
reports
  * the margin of the configured grid: the thinnest region per pool, flavour
    and channel, and every cell below --min-jets;
  * a proposed grid: starting from the candidate edges (--pt-edges, one or two
    |eta| bins), the bin with the highest pT that fails is merged into its lower
    neighbour until all bins pass; the grid with more cells wins, a tie goes to
    two |eta| bins.

Before a production exists the rule can be projected from another era's
preselection: --preselection-dir points at that era's outputs and --scale
multiplies the counts of each input process by its expected statistics ratio
(e.g. nevents of the era over nevents of the reference). The working points
are always the configured ones, so the region boundaries are the era's.

Usage (from the TauFakeFactors directory):
    python3 scripts/btag_binning_scan.py --config-dir configs/btag_efficiency/2017
    python3 scripts/btag_binning_scan.py --config-dir configs/btag_efficiency/2016preVFP \\
        --preselection-dir /ceph/sdaigler/CMS_FF/btag_eff_v3/preselection/2018 \\
        --scale-json scale_2016preVFP.json --allow-missing-inputs
"""
import argparse
import glob
import json
import os
import sys
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from helper.config_loading import load_config  # noqa: E402

DEFAULT_PT_EDGES = [20.0, 30.0, 50.0, 70.0, 100.0, 140.0, 200.0, 300.0, 1000.0]
DEFAULT_ETA_EDGES = [0.0, 1.4, 2.4]


def channel_configs(config_dir: str, channels: Optional[Sequence[str]]) -> Dict[str, Dict]:
    """Load btag_efficiency_<channel>.yaml (or the shared config pinned to the channel) per channel."""
    shared_path = os.path.join(config_dir, "btag_efficiency.yaml")
    shared = load_config(shared_path)
    names = list(channels) if channels else list(shared.get("channels") or [shared["channel"]])
    configs = {}
    for channel in names:
        path = os.path.join(config_dir, f"btag_efficiency_{channel}.yaml")
        config = load_config(path) if os.path.exists(path) else dict(shared)
        config["channel"] = channel
        config.pop("channels", None)
        configs[channel] = config
    return configs


def wps_loose_to_tight(config: Dict) -> List[Tuple[str, float]]:
    return sorted(config["btag_working_points"].items(), key=lambda item: item[1])


def region_names(wps: List[Tuple[str, float]]) -> List[str]:
    """Disjoint score regions, loosest first: fail <WP0>, <WP0> not <WP1>, ..., pass <WPn>."""
    names = [f"fail {wps[0][0]}"]
    names += [f"{low} not {high}" for (low, _), (high, _) in zip(wps[:-1], wps[1:])]
    names.append(f"pass {wps[-1][0]}")
    return names


def count_jets(files: List[str], config: Dict, pt_edges: List[float], eta_edges: List[float], tree: str) -> Dict[str, np.ndarray]:
    """Unweighted probe-jet counts per flavour: array [region, pt bin, eta bin] on the fine grid."""
    import ROOT

    wps = wps_loose_to_tight(config)
    flavors = config["jet_flavor_categories"]
    pt_col, eta_col = config["jet_pt_column"], config["jet_eta_column"]
    flavor_col, score_col = config["jet_flavor_column"], config["jet_btag_column"]
    n_regions = len(wps) + 1
    result = {name: np.zeros((n_regions, len(pt_edges) - 1, len(eta_edges) - 1)) for name in flavors}
    if not files:
        return result

    chain = ROOT.TChain(tree)
    for path in files:
        chain.Add(path)
    rdf = ROOT.RDataFrame(chain)
    pt_arr = np.array(pt_edges, dtype=float)
    eta_arr = np.array(eta_edges, dtype=float)
    # region index of each jet = number of WPs it passes (score >= cut, as in the measurement)
    passed = " + ".join(f"({score_col} >= {cut})" for _, cut in wps)
    rdf = rdf.Define("_scan_region", f"ROOT::VecOps::RVec<int>({passed})")
    booked = {}
    for name, flavor_id in flavors.items():
        for region in range(n_regions):
            mask = f"(({flavor_col} == {flavor_id}) && (_scan_region == {region}))"
            pt_name, eta_name = f"_scan_pt_{name}_{region}", f"_scan_eta_{name}_{region}"
            node = rdf.Define(pt_name, f"{pt_col}[{mask}]").Define(eta_name, f"abs({eta_col}[{mask}])")
            model = ROOT.RDF.TH2DModel(f"h_{name}_{region}", "", len(pt_arr) - 1, pt_arr, len(eta_arr) - 1, eta_arr)
            booked[(name, region)] = node.Histo2D(model, pt_name, eta_name)
    for (name, region), histogram in booked.items():
        h = histogram.GetValue()
        for i_pt in range(len(pt_edges) - 1):
            for i_eta in range(len(eta_edges) - 1):
                result[name][region, i_pt, i_eta] = h.GetBinContent(i_pt + 1, i_eta + 1)
    return result


def rebin(counts: np.ndarray, fine_pt: List[float], fine_eta: List[float], pt: List[float], eta: List[float]) -> np.ndarray:
    """Sum fine-grid counts [region, pt, eta] into a coarser grid whose edges are a subset of the fine ones."""
    for edge in pt:
        if not any(np.isclose(edge, fine_pt)):
            raise ValueError(f"pT edge {edge} is not on the fine grid {fine_pt}")
    for edge in eta:
        if not any(np.isclose(edge, fine_eta)):
            raise ValueError(f"|eta| edge {edge} is not on the fine grid {fine_eta}")
    pt_index = [int(np.argmin(np.abs(np.array(fine_pt) - e))) for e in pt]
    eta_index = [int(np.argmin(np.abs(np.array(fine_eta) - e))) for e in eta]
    out = np.zeros((counts.shape[0], len(pt) - 1, len(eta) - 1))
    for i in range(len(pt) - 1):
        for j in range(len(eta) - 1):
            out[:, i, j] = counts[:, pt_index[i]:pt_index[i + 1], eta_index[j]:eta_index[j + 1]].sum(axis=(1, 2))
    return out


def thinnest(grids: Dict[str, np.ndarray]) -> float:
    """Smallest region count over all channels and cells."""
    return min(float(grid.min()) for grid in grids.values()) if grids else float("inf")


def propose(counts: Dict[str, np.ndarray], fine_pt: List[float], fine_eta: List[float], pt_edges: List[float],
            eta_options: List[List[float]], min_jets: float) -> Optional[Tuple[List[float], List[float]]]:
    """Finest grid passing the rule in every channel (greedy merge from high pT), or None."""
    best = None
    for eta in eta_options:
        pt = list(pt_edges)
        while True:
            per_channel = {ch: rebin(c, fine_pt, fine_eta, pt, eta) for ch, c in counts.items()}
            per_bin = [min(float(grid[:, i, :].min()) for grid in per_channel.values()) for i in range(len(pt) - 1)]
            failing = [i for i, value in enumerate(per_bin) if value < min_jets]
            if not failing:
                candidate = (pt, list(eta))
                break
            if len(pt) == 2:
                candidate = None
                break
            i = max(failing)
            del pt[i if i > 0 else 1]
        if candidate is None:
            continue
        cells = (len(candidate[0]) - 1) * (len(candidate[1]) - 1)
        key = (cells, len(candidate[1]))
        if best is None or key > best[0]:
            best = (key, candidate)
    return best[1] if best else None


def fmt_edges(edges: Sequence[float]) -> str:
    return "[" + ", ".join(f"{e:g}" for e in edges) + "]"


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--config-dir", required=True, help="era directory, e.g. configs/btag_efficiency/2018")
    parser.add_argument("--channels", nargs="*", default=None, help="channels (default: the 'channels' list of btag_efficiency.yaml)")
    parser.add_argument("--preselection-dir", default=None,
                        help="directory holding <channel>/<process>*.root (default: {file_path}/preselection/{era})")
    parser.add_argument("--scale", action="append", default=[], metavar="PROCESS=FACTOR",
                        help="multiply the counts of an input process (projection to another era's statistics)")
    parser.add_argument("--scale-json", default=None, help="JSON {process: factor}, applied before --scale")
    parser.add_argument("--min-jets", type=float, default=5.0, help="minimum jets per score region and cell (default 5)")
    parser.add_argument("--pt-edges", type=float, nargs="+", default=DEFAULT_PT_EDGES, help="candidate pT edges")
    parser.add_argument("--eta-edges", type=float, nargs="+", default=DEFAULT_ETA_EDGES,
                        help="candidate |eta| edges; the proposal tries these and a single |eta| bin")
    parser.add_argument("--allow-missing-inputs", action="store_true", help="count a missing input process as empty")
    parser.add_argument("--threads", type=int, default=4, help="ROOT implicit-MT threads (default 4)")
    parser.add_argument("--json", default=None, help="write the margins and proposals to this JSON file")
    args = parser.parse_args(argv)

    configs = channel_configs(args.config_dir, args.channels)
    first = next(iter(configs.values()))
    era = first["era"]
    base_dir = args.preselection_dir or os.path.join(first["file_path"], "preselection", era)

    scale: Dict[str, float] = {}
    if args.scale_json:
        with open(args.scale_json) as handle:
            scale.update({k: float(v) for k, v in json.load(handle).items()})
    for item in args.scale:
        name, _, value = item.partition("=")
        scale[name] = float(value)

    # fine grid: the candidate edges plus every configured edge, so that both can be summed from it
    fine_pt, fine_eta = set(args.pt_edges), set(args.eta_edges)
    for config in configs.values():
        fine_pt.update(config["jet_pt_bins"])
        fine_eta.update(config["jet_eta_bins"])
        for pool in config["processes"].values():
            for key in ("jet_pt_bins", "jet_eta_bins"):
                (fine_pt if key == "jet_pt_bins" else fine_eta).update((pool or {}).get(key, []))
            for bins in ((pool or {}).get("flavor_bins") or {}).values():
                fine_pt.update((bins or {}).get("jet_pt_bins", []))
                fine_eta.update((bins or {}).get("jet_eta_bins", []))
    fine_pt, fine_eta = sorted(fine_pt), sorted(fine_eta)
    eta_options = [sorted(args.eta_edges)]
    if len(eta_options[0]) > 2:
        eta_options.append([eta_options[0][0], eta_options[0][-1]])

    import ROOT
    ROOT.gROOT.SetBatch(True)
    if args.threads > 0:
        ROOT.EnableImplicitMT(args.threads)

    # count once per (channel, input process)
    process_counts: Dict[Tuple[str, str], Dict[str, np.ndarray]] = {}
    for channel, config in configs.items():
        members = sorted({m for name, pool in config["processes"].items() for m in (pool or {}).get("input_processes", [name])})
        for member in members:
            files = sorted(glob.glob(os.path.join(base_dir, channel, f"{member}*.root")))
            if not files and not args.allow_missing_inputs:
                raise SystemExit(f"no preselection file for '{member}' in {os.path.join(base_dir, channel)}")
            counts = count_jets(files, config, fine_pt, fine_eta, config.get("tree", "ntuple"))
            factor = scale.get(member, 1.0)
            process_counts[(channel, member)] = {flavor: grid * factor for flavor, grid in counts.items()}
            print(f"counted {channel}/{member}: {len(files)} file(s), scale {factor:g}", file=sys.stderr)

    regions = region_names(wps_loose_to_tight(first))
    report = {"era": era, "preselection_dir": base_dir, "min_jets": args.min_jets, "regions": regions, "pools": {}}
    pool_names = sorted({name for config in configs.values() for name in config["processes"]})
    total_below = 0
    for pool in pool_names:
        for flavor in first["jet_flavor_categories"]:
            per_channel_fine, configured, margins, below = {}, {}, {}, []
            for channel, config in configs.items():
                pool_conf = config["processes"].get(pool)
                if pool_conf is None:
                    continue
                members = pool_conf.get("input_processes", [pool])
                per_channel_fine[channel] = sum(process_counts[(channel, m)][flavor] for m in members)
                bins = (pool_conf.get("flavor_bins") or {}).get(flavor) or {}
                pt = list(bins.get("jet_pt_bins", pool_conf.get("jet_pt_bins", config["jet_pt_bins"])))
                eta = list(bins.get("jet_eta_bins", pool_conf.get("jet_eta_bins", config["jet_eta_bins"])))
                configured[channel] = (pt, eta)
                grid = rebin(per_channel_fine[channel], fine_pt, fine_eta, pt, eta)
                margins[channel] = float(grid.min())
                for r, i, j in zip(*np.where(grid < args.min_jets)):
                    below.append(f"{channel} {regions[r]} pt=[{pt[i]:g},{pt[i + 1]:g}) eta=[{eta[j]:g},{eta[j + 1]:g}): {grid[r, i, j]:.1f}")
            if not per_channel_fine:
                continue
            total_below += len(below)
            proposal = propose(per_channel_fine, fine_pt, fine_eta, sorted(args.pt_edges), eta_options, args.min_jets)
            proposal_margin = None
            if proposal:
                proposal_margin = thinnest({ch: rebin(c, fine_pt, fine_eta, *proposal) for ch, c in per_channel_fine.items()})
            grids = {tuple(map(tuple, g)) for g in configured.values()}
            print(f"\n== {pool} / {flavor}")
            for grid in grids:
                chans = [ch for ch, g in configured.items() if tuple(map(tuple, g)) == grid]
                print(f"   configured pt {fmt_edges(grid[0])} eta {fmt_edges(grid[1])} ({', '.join(chans)})")
            print("   margin  " + "  ".join(f"{ch} {m:.1f}" for ch, m in margins.items()))
            for line in below:
                print(f"   BELOW {args.min_jets:g}: {line}")
            if proposal:
                print(f"   proposal pt {fmt_edges(proposal[0])} eta {fmt_edges(proposal[1])}  (margin {proposal_margin:.1f})")
            else:
                print(f"   proposal: none, even one cell holds fewer than {args.min_jets:g} jets in some region")
            report["pools"].setdefault(pool, {})[flavor] = {
                "configured": {ch: {"pt": g[0], "eta": g[1]} for ch, g in configured.items()},
                "margin": margins,
                "below": below,
                "proposal": {"pt": proposal[0], "eta": proposal[1], "margin": proposal_margin} if proposal else None,
            }

    print(f"\n{total_below} configured cell-region(s) below {args.min_jets:g} jets")
    if args.json:
        with open(args.json, "w") as handle:
            json.dump(report, handle, indent=1)
    return 1 if total_below else 0


if __name__ == "__main__":
    sys.exit(main())
