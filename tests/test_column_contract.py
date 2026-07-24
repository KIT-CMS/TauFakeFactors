"""Static contract regression test for the b-tag efficiency calculator configs.

The btag efficiency calculator (``btag_efficiency.py``) reads jet columns
whose names are configured via ``jet_pt_column``/``jet_eta_column``/
``jet_flavor_column``/``jet_btag_column`` (and, optionally, an extra
``jet_selection`` C++ expression) out of the ROOT files produced by
``preselection.py``. ``preselection.py`` only ever writes out the columns
listed under ``output_features`` in the paired ``preselection_<channel>.yaml``
(see ``preselection.py`` around the ``Snapshot(..., output_features)`` call) --
so any column requested by a calculator config that isn't also present in
``output_features`` will blow up at runtime with a missing-branch error.

This test parses the YAML configs directly (no ROOT/RDataFrame involved) and
asserts that contract holds for every era/channel combination that currently
exists, plus pins down the two known column-naming conventions:

- Run-3 eras (2022preEE, 2022postEE, 2023preBPix, 2023postBPix, 2024, 2025)
  write lowercase ``jet_*_vec`` columns.
- 2018 (added later, in Task 19) is expected to use the ``btag_probe_jet_*``
  naming instead; this is asserted only once the ``2018`` config directory
  actually exists, and Run-3 configs are asserted to never drift onto that
  naming scheme.
"""

import re
from pathlib import Path
from typing import Dict, Iterator, List, Set, Tuple

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[1]
BTAG_EFFICIENCY_DIR = REPO_ROOT / "configs" / "btag_efficiency"

# Column keys read directly by btag_efficiency.py (see
# calculate_efficiency_histograms() in btag_efficiency.py).
JET_COLUMN_KEYS = (
    "jet_pt_column",
    "jet_eta_column",
    "jet_flavor_column",
    "jet_btag_column",
)

# The mc_weights sum is always snapshotted under this name by preselection.py
# and is required by any downstream calculator config.
WEIGHT_COLUMN = "weight"

# Identifier regex used to pull candidate column names out of a
# `jet_selection` C++ boolean expression.
IDENTIFIER_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")

# C++ keywords/functions that can legally appear in a `jet_selection`
# expression but never denote a branch/column name.
CPP_NON_COLUMN_TOKENS = {
    "abs",
    "fabs",
    "std",
    "max",
    "min",
    "true",
    "false",
}

RUN3_ERAS = (
    "2022preEE",
    "2022postEE",
    "2023preBPix",
    "2023postBPix",
    "2024",
    "2025",
)

# The lowercase column names Run-3 preselections actually snapshot.
RUN3_JET_COLUMNS = {
    "jet_pt_column": "jet_pt_vec",
    "jet_eta_column": "jet_eta_vec",
    "jet_flavor_column": "jet_hadronflavour_vec",
    "jet_btag_column": "jet_btag_value_vec",
}

# The naming contract 2018 (Task 19) is expected to use instead.
CONTRACT_2018_JET_COLUMNS = {
    "jet_pt_column": "btag_probe_jet_pt",
    "jet_eta_column": "btag_probe_jet_eta",
    "jet_flavor_column": "btag_probe_jet_hadron_flavour",
    "jet_btag_column": "btag_probe_jet_upart",
}


def _load_yaml(path: Path) -> Dict:
    with open(path, "r") as config_file:
        return yaml.safe_load(config_file)


def _identifiers_in_jet_selection(jet_selection: str) -> Set[str]:
    """Extract candidate column-name identifiers from a jet_selection expression.

    ``jet_selection`` is a raw C++ boolean expression (e.g.
    ``"jet_pt_vec > 20 && abs(jet_eta_vec) < 2.5"``). Numeric literals never
    match the identifier regex (it requires a leading letter/underscore), and
    C++ keywords/operators like ``abs``/``&&`` are filtered out separately.
    """
    tokens = set(IDENTIFIER_RE.findall(jet_selection))
    return {token for token in tokens if token not in CPP_NON_COLUMN_TOKENS}


def get_config_channels(config: Dict) -> List[str]:
    """Return the list of channels a calculator config declares.

    Mirrors ``btag_efficiency.get_config_channels``: supports both the
    shared multi-channel ``channels`` list and the single-channel ``channel``
    key.
    """
    channels = config.get("channels")
    if channels is not None:
        return [str(channel) for channel in channels]

    channel = config.get("channel")
    if channel is None:
        raise ValueError("Config must specify either 'channel' or 'channels'.")
    return [str(channel)]


def iter_era_configs() -> Iterator[Tuple[str, Path, Path]]:
    """Yield ``(era, presel_cfg_dir, calc_cfg_path)`` for every calculator config.

    ``calc_cfg_path`` is every ``btag_efficiency*.yaml`` file found directly
    under an era directory (the shared ``btag_efficiency.yaml`` plus any
    per-channel ``btag_efficiency_<channel>.yaml`` files). ``presel_cfg_dir``
    is that same era directory, since ``preselection_<channel>.yaml`` files
    are always colocated with the calculator configs -- callers resolve the
    preselection config for a given channel via
    ``presel_cfg_dir / f"preselection_{channel}.yaml"``.
    """
    if not BTAG_EFFICIENCY_DIR.is_dir():
        return
    for era_dir in sorted(BTAG_EFFICIENCY_DIR.iterdir()):
        if not era_dir.is_dir():
            continue
        era = era_dir.name
        for calc_cfg_path in sorted(era_dir.glob("btag_efficiency*.yaml")):
            yield era, era_dir, calc_cfg_path


def _required_columns(calc_config: Dict) -> Set[str]:
    """Columns a calculator config needs present in output_features."""
    required: Set[str] = {WEIGHT_COLUMN}
    for key in JET_COLUMN_KEYS:
        if key in calc_config:
            required.add(str(calc_config[key]))

    jet_selection = str(calc_config.get("jet_selection", "") or "").strip()
    if jet_selection:
        required |= _identifiers_in_jet_selection(jet_selection)

    return required


ERA_CONFIG_CASES = list(iter_era_configs())
ERA_CONFIG_CASE_IDS = [
    f"{era}/{calc_cfg_path.name}" for era, _, calc_cfg_path in ERA_CONFIG_CASES
]


@pytest.mark.parametrize(
    "era,presel_cfg_dir,calc_cfg_path",
    ERA_CONFIG_CASES,
    ids=ERA_CONFIG_CASE_IDS,
)
def test_calculator_columns_are_written_by_preselection(
    era: str, presel_cfg_dir: Path, calc_cfg_path: Path
) -> None:
    """Every column a calculator config asks for must be in output_features.

    Covers the four jet_*_column keys, every identifier referenced in
    jet_selection, and the "weight" column, checked against every channel
    the calculator config declares.
    """
    calc_config = _load_yaml(calc_cfg_path)
    required_columns = _required_columns(calc_config)

    for channel in get_config_channels(calc_config):
        presel_cfg_path = presel_cfg_dir / f"preselection_{channel}.yaml"
        assert presel_cfg_path.is_file(), (
            f"{calc_cfg_path} declares channel '{channel}' but "
            f"{presel_cfg_path} does not exist"
        )
        presel_config = _load_yaml(presel_cfg_path)
        output_features = set(presel_config.get("output_features") or [])

        missing = sorted(required_columns - output_features)
        assert not missing, (
            f"{calc_cfg_path} (channel={channel}) requires columns {missing} "
            f"that are not in {presel_cfg_path}'s output_features "
            f"({sorted(output_features)})"
        )


def test_run3_eras_keep_lowercase_jet_vec_columns() -> None:
    """Run-3 calculator configs must use the lowercase jet_*_vec names.

    Also guards against an accidental future rename onto the 2018
    btag_probe_jet_* naming contract.
    """
    run3_cases = [
        (era, calc_cfg_path)
        for era, _, calc_cfg_path in ERA_CONFIG_CASES
        if era in RUN3_ERAS
    ]
    assert run3_cases, "expected to find Run-3 btag_efficiency configs"

    for era, calc_cfg_path in run3_cases:
        calc_config = _load_yaml(calc_cfg_path)
        for key, expected in RUN3_JET_COLUMNS.items():
            actual = calc_config.get(key)
            assert actual == expected, (
                f"{calc_cfg_path}: expected {key}={expected!r}, got {actual!r}"
            )
            assert actual != CONTRACT_2018_JET_COLUMNS[key], (
                f"{calc_cfg_path}: {key} must not be renamed to the 2018 "
                f"btag_probe_jet_* contract"
            )


def test_2018_calculator_uses_btag_probe_jet_columns_once_added() -> None:
    """2018 (Task 19) must use the btag_probe_jet_* naming contract.

    The 2018 config directory does not exist yet, so this is a no-op today
    and will start enforcing the contract the moment Task 19 adds it.
    """
    era_dir = BTAG_EFFICIENCY_DIR / "2018"
    if not era_dir.is_dir():
        pytest.skip("configs/btag_efficiency/2018 does not exist yet (Task 19)")

    calc_cfg_paths = sorted(era_dir.glob("btag_efficiency*.yaml"))
    assert calc_cfg_paths, f"no btag_efficiency*.yaml configs found under {era_dir}"

    for calc_cfg_path in calc_cfg_paths:
        calc_config = _load_yaml(calc_cfg_path)
        for key, expected in CONTRACT_2018_JET_COLUMNS.items():
            actual = calc_config.get(key)
            assert actual == expected, (
                f"{calc_cfg_path}: expected {key}={expected!r}, got {actual!r}"
            )
