#!/usr/bin/env python3
"""
Regenerate the per-channel preselection_{ch}.yaml of the fake factor measurement from the unified
config YAML files, for each (era, channel) pair:
  fake_factors/configs/smhtt/{era}/preselection_{ch}.yaml

Sections written:
  processes       -- sample lists and tau_gen_modes (common/datasets.yaml); the combined eras
                     list every sample under era_overrides with its sub-era, 2025 data also takes
                     the 2026 samples
  event_selection -- CROWN's presel_mask (weights, region masks and output columns come from the
                     CROWN friend and common/variables.yaml's selection_masks, nothing is built here)

Normally invoked via `run.py GEN-TFF-CONFIGS`.
"""

import argparse
import logging
import sys
from pathlib import Path

import yaml

_REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_REPO_ROOT))

from CustomLogging import setup_logging
from common.analysis import DATA_ERA_MERGE, PRESELECTION_MASK_BRANCH

logger = setup_logging(logger=logging.getLogger(__name__))

TFF_CONFIG_ROOT = _REPO_ROOT / "fake_factors" / "configs" / "smhtt"
_DATASETS_YAML = _REPO_ROOT / "common" / "datasets.yaml"

# FF config era -> the sub-eras of datasets.yaml it combines
FF_ERAS = {
    "2022": ["2022preEE", "2022postEE"],
    "2023": ["2023preBPix", "2023postBPix"],
    "2024": ["2024"],
    "2025": ["2025"],
}
CHANNELS = ["et", "mt", "tt"]

_PROCESS_ORDER = ["DYjets", "ST", "Wjets", "data", "diboson", "ttbar"]


def build_processes(era: str, channel: str, backgrounds: dict) -> dict:
    sub_eras = FF_ERAS[era]
    combined = len(sub_eras) > 1
    processes = {}
    for proc in _PROCESS_ORDER:
        by_era = {}
        for sub_era in [*sub_eras, *(DATA_ERA_MERGE.get(era, ()) if proc == "data" else ())]:
            info = backgrounds[sub_era].get(proc)
            if info is not None:
                samples = info["by_channel"].get(channel, []) if proc == "data" else info.get("samples", [])
                by_era[sub_era] = list(samples)
        if not by_era:
            continue
        info = backgrounds[sub_eras[0]][proc]
        entry = {} if combined else {"samples": by_era[era]}
        overrides = {s: sub_era for sub_era, samples in by_era.items() if combined or sub_era != era for s in samples}
        if overrides:
            entry["era_overrides"] = overrides
        entry["tau_gen_modes"] = info.get("tau_gen_modes", ["all"])
        processes[proc] = entry
    return processes


def main() -> None:
    parser = argparse.ArgumentParser(description="Regenerate the preselection YAMLs of the fake factor measurement from the unified config.")
    parser.add_argument("--era", default="all", help=f"Comma-separated eras of {list(FF_ERAS)}, or 'all' (default: all)")
    parser.add_argument("--channels", default="all", help="Comma-separated channels, or 'all' (default: all)")
    parser.add_argument("--dry-run", action="store_true", help="Print what would be written without writing files")
    args = parser.parse_args()

    eras = list(FF_ERAS) if args.era == "all" else [e.strip() for e in args.era.split(",")]
    channels = CHANNELS if args.channels == "all" else [c.strip() for c in args.channels.split(",")]
    for values, valid, kind in [(eras, FF_ERAS, "era"), (channels, CHANNELS, "channel")]:
        if invalid := [v for v in values if v not in valid]:
            parser.error(f"Unknown {kind}(s): {invalid}. Valid: {list(valid)}")

    datasets = yaml.safe_load(_DATASETS_YAML.read_text())
    backgrounds = {era: {**procs, **datasets["data"].get(era, {})} for era, procs in datasets["backgrounds"].items()}
    for era in eras:
        for channel in channels:
            path = TFF_CONFIG_ROOT / era / f"preselection_{channel}.yaml"
            text = yaml.dump(
                {
                    "channel": channel,
                    "processes": build_processes(era, channel, backgrounds),
                    "event_selection": {PRESELECTION_MASK_BRANCH: f"({PRESELECTION_MASK_BRANCH} > 0.5)"},
                },
                default_flow_style=False, sort_keys=False, allow_unicode=True,
            )
            if args.dry_run:
                logger.info(f"[dry-run] would write {path}")
                print(text)
            else:
                path.write_text(text)
                logger.info(f"updated {path.relative_to(_REPO_ROOT)}")


if __name__ == "__main__":
    main()
