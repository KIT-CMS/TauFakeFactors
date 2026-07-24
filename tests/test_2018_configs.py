"""Tests for the 2018 UL (NanoAOD v15, UParT) b-tag efficiency configs (Task 19).

Covers, for ``configs/btag_efficiency/2018/``:

* every preselection / calculator config loads via ``helper.functions.load_config``
  for all three channels (et, mt, tt);
* the probe-jet column contract (btag_probe_jet_*) and the "no literal user
  paths" rule (these also live in the shared contract / path-override tests --
  re-asserted here so this module is a self-contained 2018 gate);
* no weight block references Zpt (ZPtReweighting is Run-3-only) and none applies
  a b-tag weight (this IS the b-tag efficiency measurement);
* the pinned BTV working-point validation in ``btag_efficiency.py`` (agree ->
  ok, disagree -> fatal), including that the committed 2018 WPs match the
  pinned payload;
* selection-contract parity: every threshold in ``selection_contract_2018_v1.yaml``
  appears with the same value in the matching ``preselection_<ch>.yaml``
  event_selection, and the parity check has teeth (a tampered threshold fails);
* cross-repo provenance compatibility: a synthetic install produced by
  ``btag_efficiency.build_provenance_meta`` + ``install_provenance_payload`` is
  accepted by the bbtautau CROWN gate ``btag_payloads.require_validated_payload``.
"""

import gzip
import importlib.util
import json
import os
import re
import sys
from pathlib import Path

import pytest
import yaml

import helper.functions as func

REPO_ROOT = Path(__file__).resolve().parents[1]
ERA_DIR = REPO_ROOT / "configs" / "btag_efficiency" / "2018"
CONTRACT_PATH = ERA_DIR / "selection_contract_2018_v1.yaml"
CHANNELS = ("et", "mt", "tt")

# Sibling CROWN checkout holding the bbtautau analysis config (the payload gate
# lives there). Overridable via env for non-standard layouts.
CROWN_ROOT = Path(
    os.environ.get(
        "TFF_BBTAUTAU_CROWN_ROOT", str(REPO_ROOT.parent / "KingMaker" / "CROWN")
    )
)

USER_PATH_RE = re.compile(r"(/ceph|/store/user)/[A-Za-z0-9]+")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def _load_raw(path: Path) -> dict:
    with open(path, "r") as handle:
        return yaml.safe_load(handle)


def _presel_path(channel: str) -> Path:
    return ERA_DIR / f"preselection_{channel}.yaml"


def _calc_paths():
    return sorted(ERA_DIR.glob("btag_efficiency*.yaml"))


def _weight_blocks(config: dict):
    """Yield every weight-block value string of a preselection config."""
    for block_key in ("mc_weights", "emb_weights"):
        block = config.get(block_key) or {}
        for key, value in block.items():
            yield block_key, key, str(value)


# ---------------------------------------------------------------------------
# config loading
# ---------------------------------------------------------------------------
def test_era_dir_exists():
    assert ERA_DIR.is_dir(), f"{ERA_DIR} must exist (Task 19)"


@pytest.mark.parametrize("channel", CHANNELS)
def test_preselection_loads_via_load_config(channel):
    config = func.load_config(str(_presel_path(channel)))
    assert config["era"] == "2018"
    assert config["channel"] == channel
    assert config["tree"] == "ntuple"
    assert config.get("event_selection")
    assert config.get("mc_weights")
    assert config.get("output_features")


@pytest.mark.parametrize("channel", CHANNELS)
def test_calculator_loads_via_load_config(channel):
    config = func.load_config(str(ERA_DIR / f"btag_efficiency_{channel}.yaml"))
    assert config["era"] == "2018"
    assert config["channel"] == channel


def test_shared_calculator_loads_via_load_config():
    config = func.load_config(str(ERA_DIR / "btag_efficiency.yaml"))
    assert list(config["channels"]) == list(CHANNELS)


# ---------------------------------------------------------------------------
# probe-jet column contract + output_features
# ---------------------------------------------------------------------------
CONTRACT_2018_JET_COLUMNS = {
    "jet_pt_column": "btag_probe_jet_pt",
    "jet_eta_column": "btag_probe_jet_eta",
    "jet_flavor_column": "btag_probe_jet_hadron_flavour",
    "jet_btag_column": "btag_probe_jet_upart",
}


@pytest.mark.parametrize("calc_path", _calc_paths(), ids=lambda p: p.name)
def test_calculator_uses_probe_jet_contract(calc_path):
    config = _load_raw(calc_path)
    for key, expected in CONTRACT_2018_JET_COLUMNS.items():
        assert config.get(key) == expected, f"{calc_path.name}: {key}"
    assert config.get("weight_column") == "weight"


@pytest.mark.parametrize("channel", CHANNELS)
def test_output_features_carry_probe_vectors_and_weight(channel):
    config = _load_raw(_presel_path(channel))
    features = set(config.get("output_features") or [])
    required = {"weight"} | set(CONTRACT_2018_JET_COLUMNS.values())
    assert required <= features, f"missing {sorted(required - features)}"


# ---------------------------------------------------------------------------
# no literal user paths
# ---------------------------------------------------------------------------
def test_no_literal_user_paths_in_2018_configs():
    violations = []
    for path in sorted(ERA_DIR.rglob("*.yaml")):
        text = path.read_text()
        for match in USER_PATH_RE.finditer(text):
            violations.append((path.name, match.group(0)))
    assert not violations, f"literal user paths found: {violations}"


def test_common_settings_paths_are_relative_placeholders():
    common = _load_raw(ERA_DIR / "common_settings.yaml")
    for key in ("ntuple_path", "output_path", "file_path"):
        value = str(common[key])
        assert not value.startswith("/"), f"{key} must be a relative placeholder"
        assert not value.startswith("root://"), f"{key} must be a relative placeholder"


# ---------------------------------------------------------------------------
# no Zpt reweighting, no b-tag weight in any 2018 weight block
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("channel", CHANNELS)
def test_no_zpt_in_weight_blocks(channel):
    config = _load_raw(_presel_path(channel))
    zpt_re = re.compile(r"zpt|z_pt|zptreweight", re.IGNORECASE)
    for block_key, key, value in _weight_blocks(config):
        assert not zpt_re.search(key), f"{channel}: {block_key}.{key} references zpt"
        assert not zpt_re.search(
            value
        ), f"{channel}: {block_key}.{key} value references zpt"


@pytest.mark.parametrize("channel", CHANNELS)
def test_no_btag_weight_in_weight_blocks(channel):
    config = _load_raw(_presel_path(channel))
    btag_re = re.compile(r"btag", re.IGNORECASE)
    for block_key, key, value in _weight_blocks(config):
        assert not btag_re.search(key), f"{channel}: {block_key}.{key} references btag"
        assert not btag_re.search(
            value
        ), f"{channel}: {block_key}.{key} value references btag"


# ---------------------------------------------------------------------------
# pinned BTV working-point validation
# ---------------------------------------------------------------------------
def _write_synthetic_btv_payload(path: str, wps: dict) -> None:
    payload = {
        "schema_version": 2,
        "corrections": [
            {
                "name": "UParTAK4_wp_values",
                "data": {
                    "nodetype": "category",
                    "input": "working_point",
                    "content": [{"key": k, "value": v} for k, v in wps.items()],
                },
            }
        ],
    }
    with gzip.open(path, "wt") as handle:
        json.dump(payload, handle)


PINNED_WPS = {"L": 0.0308, "M": 0.161, "T": 0.5405, "XT": 0.6992, "XXT": 0.9655}


def test_resolve_working_points_agrees(tmp_path):
    import btag_efficiency

    source = str(tmp_path / "btagging.json.gz")
    _write_synthetic_btv_payload(source, PINNED_WPS)
    config = {
        "btag_working_points_source": source,
        "btag_working_points": dict(PINNED_WPS),
    }
    assert btag_efficiency.resolve_working_points(config) == PINNED_WPS


def test_resolve_working_points_mismatch_is_fatal(tmp_path):
    import btag_efficiency

    source = str(tmp_path / "btagging.json.gz")
    _write_synthetic_btv_payload(source, PINNED_WPS)
    tampered = dict(PINNED_WPS)
    tampered["M"] = 0.5  # drift the medium WP
    config = {"btag_working_points_source": source, "btag_working_points": tampered}
    with pytest.raises(RuntimeError, match="working points"):
        btag_efficiency.resolve_working_points(config)


def test_resolve_working_points_missing_block_is_fatal(tmp_path):
    import btag_efficiency

    source = str(tmp_path / "btagging.json.gz")
    _write_synthetic_btv_payload(source, PINNED_WPS)
    with pytest.raises(RuntimeError, match="both must be present"):
        btag_efficiency.resolve_working_points({"btag_working_points_source": source})


def test_committed_2018_wps_match_configured_source():
    """The committed 2018 WP block agrees with the pinned BTV payload."""
    import btag_efficiency

    common = _load_raw(ERA_DIR / "common_settings.yaml")
    source = common["btag_working_points_source"]
    if not os.path.isfile(source):
        pytest.skip(f"pinned BTV payload not available at {source}")
    calc = _load_raw(ERA_DIR / "btag_efficiency.yaml")
    config = {
        "btag_working_points_source": source,
        "btag_working_points": calc["btag_working_points"],
    }
    assert btag_efficiency.resolve_working_points(config) == calc["btag_working_points"]


# ---------------------------------------------------------------------------
# selection-contract parity
# ---------------------------------------------------------------------------
def _num_present(value, text: str) -> bool:
    """True if the numeric ``value`` appears as a standalone number in ``text``."""
    candidates = set()
    fvalue = float(value)
    if fvalue.is_integer():
        candidates.add(str(int(fvalue)))
    candidates.add(repr(fvalue))
    candidates.add(str(value))
    for cand in candidates:
        pattern = r"(?<![\d.])" + re.escape(cand) + r"(?![\d.])"
        if re.search(pattern, text):
            return True
    return False


def parity_violations(contract: dict, presel: dict, channel: str):
    """Return the list of contract thresholds NOT reproduced by the preselection.

    Checks the event_selection block against the contract's per-channel
    thresholds: trigger flag(s) + trigger pt, hadronic-tau pt (single-lepton
    channels), decay modes, the three tau-ID working-point branches, the
    charge requirement, the extra-lepton vetoes, and -- for the single-lepton
    channels -- the light-lepton isolation and transverse-mass windows.
    """
    ch = contract["channels"][channel]
    selection = presel.get("event_selection") or {}
    text = " ".join(str(v) for v in selection.values())
    violations = []

    trigger = ch["trigger"]
    for flag in trigger["flags"]:
        if flag not in text:
            violations.append(f"trigger flag {flag}")
    if "lepton_min_pt" in trigger and not _num_present(trigger["lepton_min_pt"], text):
        violations.append(f"trigger lepton_min_pt {trigger['lepton_min_pt']}")
    if "tau_min_pt" in trigger and not _num_present(trigger["tau_min_pt"], text):
        violations.append(f"trigger tau_min_pt {trigger['tau_min_pt']}")

    tau = ch["tau"]
    if channel in ("et", "mt") and not _num_present(tau["min_pt"], text):
        violations.append(f"tau min_pt {tau['min_pt']}")
    for dm in tau["decay_modes"]:
        if not _num_present(dm, text):
            violations.append(f"decay mode {dm}")
    for disc, wp in (
        ("vsJet", tau["vs_jet_wp"]),
        ("vsEle", tau["vs_ele_wp"]),
        ("vsMu", tau["vs_mu_wp"]),
    ):
        if f"{disc}_{wp}" not in text:
            violations.append(f"tau {disc} WP {wp}")

    if "transverse_mass_max" in ch and not _num_present(
        ch["transverse_mass_max"], text
    ):
        violations.append(f"transverse_mass_max {ch['transverse_mass_max']}")

    lepton = ch.get("lepton")
    if lepton and "max_iso" in lepton and not _num_present(lepton["max_iso"], text):
        violations.append(f"lepton max_iso {lepton['max_iso']}")

    if ch.get("charge_requirement") == "opposite_sign":
        if "q_1*q_2" not in text.replace(" ", ""):
            violations.append("charge opposite_sign (q_1*q_2 < 0)")

    for veto in ch["vetoes"]:
        if veto not in text:
            violations.append(f"veto {veto}")

    return violations


def test_contract_exists_and_has_checksum():
    contract = _load_raw(CONTRACT_PATH)
    assert contract["era"] == "2018"
    assert re.fullmatch(r"[0-9a-f]{64}", contract["contract_sha256"])
    assert set(contract["channels"]) == set(CHANNELS)


@pytest.mark.parametrize("channel", CHANNELS)
def test_preselection_matches_contract(channel):
    contract = _load_raw(CONTRACT_PATH)
    presel = _load_raw(_presel_path(channel))
    violations = parity_violations(contract, presel, channel)
    assert not violations, f"{channel}: contract not reproduced: {violations}"


def test_parity_check_fails_on_threshold_edit(tmp_path):
    """The parity check has teeth: edit a threshold and it must report a violation."""
    contract = _load_raw(CONTRACT_PATH)
    presel = _load_raw(_presel_path("mt"))
    # sanity: unmodified config passes
    assert not parity_violations(contract, presel, "mt")

    # tamper the transverse-mass window (50 -> 60) in a copy
    tampered = json.loads(json.dumps(presel))
    tampered["event_selection"]["lep_mt"] = "(mt_1 < 60)"
    violations = parity_violations(contract, tampered, "mt")
    assert any("transverse_mass_max" in v for v in violations), violations


def test_parity_check_fails_on_working_point_edit():
    """Editing a tau-ID working point in the preselection breaks parity."""
    contract = _load_raw(CONTRACT_PATH)
    presel = _load_raw(_presel_path("mt"))
    tampered = json.loads(json.dumps(presel))
    # mt cuts vsMu at Tight; drop it to VLoose -> contract (Tight) no longer met
    tampered["event_selection"]["had_tau_id_vs_mu"] = "id_tau_vsMu_VLoose_2 > 0.5"
    violations = parity_violations(contract, tampered, "mt")
    assert any("vsMu" in v for v in violations), violations


# ---------------------------------------------------------------------------
# ground-truth reality check (mt): every column the mt preselection references
# must exist in the real 2018-v15 CROWN probe ntuple.
# ---------------------------------------------------------------------------
GROUND_TRUTH_MT = os.environ.get(
    "TFF_GROUND_TRUTH_MT",
    "/work/sdaigler/bbtautau/test_inputs/output_ttbar_probe_mt.root",
)

_IDENTIFIER_RE = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")


def test_mt_preselection_columns_exist_in_ground_truth_ntuple():
    """Reality check: the mt preselection reads only branches CROWN actually writes."""
    if not os.path.isfile(GROUND_TRUTH_MT):
        pytest.skip(f"ground-truth probe ntuple not available at {GROUND_TRUTH_MT}")

    import btag_efficiency

    branches = set(btag_efficiency._tree_columns(GROUND_TRUTH_MT, "ntuple"))

    presel = _load_raw(_presel_path("mt"))
    common = _load_raw(ERA_DIR / "common_settings.yaml")

    referenced = set()
    for block in ("event_selection", "mc_weights"):
        for value in (presel.get(block) or {}).values():
            for token in _IDENTIFIER_RE.findall(str(value)):
                if not token.isdigit():
                    referenced.add(token)
    referenced |= set(presel["output_features"])
    for wp in common["tau_vs_jet_wps"]:
        referenced.add(f"id_tau_vsJet_{wp}_2")
    for wp in common["tau_vs_jet_wgt_wps"]:
        referenced.add(f"id_wgt_tau_vsJet_{wp}_2")

    # "weight" is built at runtime by preselection.py, not read from the ntuple.
    referenced.discard("weight")

    missing = sorted(col for col in referenced if col not in branches)
    assert (
        not missing
    ), f"mt preselection references columns absent from CROWN ntuple: {missing}"


# ---------------------------------------------------------------------------
# cross-repo provenance compatibility with the bbtautau gate
# ---------------------------------------------------------------------------
def _import_bbtautau_gate():
    """Import btag_payloads from the sibling CROWN checkout, or skip."""
    gate_path = CROWN_ROOT / "analysis_configurations" / "bbtautau" / "btag_payloads.py"
    if not gate_path.is_file():
        pytest.skip(f"bbtautau btag_payloads.py not found at {gate_path}")
    if str(CROWN_ROOT) not in sys.path:
        sys.path.insert(0, str(CROWN_ROOT))
    spec = importlib.util.spec_from_file_location(
        "bbtautau_btag_payloads_under_test", str(gate_path)
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_synthetic_efficiency_payload(path: str, categories) -> None:
    """Write a minimal correctionlib efficiency payload with a sample_type axis."""

    def flavor_cat(value):
        return {
            "nodetype": "category",
            "input": "jet_flavor",
            "content": [
                {
                    "key": flavor,
                    "value": {
                        "nodetype": "binning",
                        "input": "jet_eta",
                        "edges": [0.0, 2.4],
                        "content": [
                            {
                                "nodetype": "binning",
                                "input": "jet_pt",
                                "edges": [20.0, 1000.0],
                                "content": [value],
                                "flow": "clamp",
                            }
                        ],
                        "flow": "error",
                    },
                }
                for flavor in (0, 4, 5)
            ],
        }

    payload = {
        "schema_version": 2,
        "corrections": [
            {
                "name": "btag_efficiency",
                "version": 1,
                "inputs": [
                    {"name": "sample_type", "type": "string"},
                    {"name": "working_point", "type": "string"},
                    {"name": "jet_flavor", "type": "int"},
                    {"name": "jet_eta", "type": "real"},
                    {"name": "jet_pt", "type": "real"},
                ],
                "output": {"name": "efficiency", "type": "real"},
                "data": {
                    "nodetype": "category",
                    "input": "sample_type",
                    "content": [
                        {
                            "key": sample_type,
                            "value": {
                                "nodetype": "category",
                                "input": "working_point",
                                "content": [
                                    {"key": wp, "value": flavor_cat(0.1)}
                                    for wp in ("L", "M", "T", "XT", "XXT")
                                ],
                            },
                        }
                        for sample_type in categories
                    ],
                },
            }
        ],
    }
    with gzip.open(path, "wt") as handle:
        json.dump(payload, handle)


def test_provenance_writer_output_accepted_by_bbtautau_gate(tmp_path):
    """The provenance my writer produces is accepted by the bbtautau gate.

    Builds a synthetic per-scope payload set (with the 13 SM sample_type
    categories), runs btag_efficiency.build_provenance_meta +
    install_provenance_payload to produce a gate-compatible payload directory
    (btag_efficiency_<scope>.json.gz + provenance.json), then runs the real
    bbtautau require_validated_payload against it -- the true cross-repo
    compatibility proof.
    """
    import btag_efficiency

    gate = _import_bbtautau_gate()

    # synthetic per-scope payloads
    source_dir = tmp_path / "source"
    source_dir.mkdir()
    scope_files = {}
    for scope in CHANNELS:
        p = str(source_dir / f"eff_{scope}.json.gz")
        _write_synthetic_efficiency_payload(p, gate.SM_BTAG_EFFICIENCY_CATEGORIES)
        scope_files[scope] = p

    config = func.load_config(str(ERA_DIR / "btag_efficiency.yaml"))
    config["production_tag"] = "unit_test_tag"
    meta = btag_efficiency.build_provenance_meta(
        config,
        config_dir=str(ERA_DIR),
        validation_status="passed",
        channels=list(CHANNELS),
    )

    target_dir = str(tmp_path / "payload")
    btag_efficiency.install_provenance_payload(scope_files, meta, target_dir)

    # provenance.json + per-scope payloads must be installed
    assert os.path.isfile(os.path.join(target_dir, "provenance.json"))
    for scope in CHANNELS:
        assert os.path.isfile(
            os.path.join(target_dir, f"btag_efficiency_{scope}.json.gz")
        )

    # the real bbtautau gate must accept it
    provenance = gate.require_validated_payload(
        target_dir, list(CHANNELS), gate.SM_BTAG_EFFICIENCY_CATEGORIES
    )
    assert provenance["validation_status"] == "passed"
    assert set(provenance["manifest"]) == {
        f"btag_efficiency_{s}.json.gz" for s in CHANNELS
    }
    # the extra provenance fields the brief requires are present
    assert provenance["production_tag"] == "unit_test_tag"
    assert (
        provenance["selection_contract"]["sha256"]
        == _load_raw(CONTRACT_PATH)["contract_sha256"]
    )
    assert set(provenance["weight_expression"]) == set(CHANNELS)
    assert "repo_commits" in provenance
    assert provenance["binning"]["jet_eta_bins"] == [0.0, 1.4, 2.4]


def test_gate_rejects_tampered_payload(tmp_path):
    """Sanity: the gate rejects a payload whose file was modified post-install."""
    import btag_efficiency

    gate = _import_bbtautau_gate()

    source_dir = tmp_path / "source"
    source_dir.mkdir()
    scope_files = {}
    for scope in CHANNELS:
        p = str(source_dir / f"eff_{scope}.json.gz")
        _write_synthetic_efficiency_payload(p, gate.SM_BTAG_EFFICIENCY_CATEGORIES)
        scope_files[scope] = p

    config = func.load_config(str(ERA_DIR / "btag_efficiency.yaml"))
    meta = btag_efficiency.build_provenance_meta(
        config,
        config_dir=str(ERA_DIR),
        validation_status="passed",
        channels=list(CHANNELS),
    )
    target_dir = str(tmp_path / "payload")
    btag_efficiency.install_provenance_payload(scope_files, meta, target_dir)

    with open(os.path.join(target_dir, "btag_efficiency_mt.json.gz"), "ab") as handle:
        handle.write(b"tampered")

    with pytest.raises(ValueError, match="checksum"):
        gate.require_validated_payload(
            target_dir, list(CHANNELS), gate.SM_BTAG_EFFICIENCY_CATEGORIES
        )
