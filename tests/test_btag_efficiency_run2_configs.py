"""Offline checks of the Run-2 b-tag efficiency configs (v7: 2016preVFP, 2016postVFP, 2017, 2018).

Runs without the probe-jet production: config loading, the sample database, pool coverage, the working points
against the pinned BTV payload, and every selection and weight expression JIT-compiled by ROOT.

    python3 -m pytest tests/test_btag_efficiency_run2_configs.py
"""
import gzip
import json
import os
import re

import pytest

from helper.config_loading import load_config

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_DIR = os.path.join(ROOT_DIR, "configs", "btag_efficiency")
DATABASE = os.path.join(ROOT_DIR, "datasets", "nanoAOD_v15", "datasets.json")
ERAS = ["2016preVFP", "2016postVFP", "2017", "2018"]
CHANNELS = ["et", "mt", "tt", "ee", "em", "mm"]
SM_SAMPLE_TYPES = {
    "hh2b2tau", "dyjets", "wjets", "ttbar", "singletop", "diboson", "electroweak_boson", "ggh_htautau",
    "vbf_htautau", "vbf_hbb", "rem_hbb", "rem_higgs", "rem_htautau", "rem_ttbar",
}
# v7 spec 7.2 / E14: [T, M] with the era values of the pinned BTV payload
WORKING_POINTS = {
    "2016preVFP": {"M": 0.1847, "T": 0.5467},
    "2016postVFP": {"M": 0.1898, "T": 0.5538},
    "2017": {"M": 0.1776, "T": 0.5755},
    "2018": {"M": 0.161, "T": 0.5405},
}
BTV_PIN = "/cvmfs/cms-griddata.cern.ch/cat/metadata/BTV/Run2-{era}-UL-NanoAODv15/2026-06-18/btagging.json.gz"
# The production smrun2_btag_eff_v1 runs the CROWN _mc lists (v7 spec 8, P2); a nick that is not in the list has
# no directory under CROWNRun and check_inputfiles aborts the whole preselection. Searched in this order, the
# ShapeSmith inventories are byte-identical copies; CROWN_SAMPLE_LIST_DIR overrides.
SAMPLE_LIST_DIRS = [
    os.environ.get("CROWN_SAMPLE_LIST_DIR", ""),
    os.path.join(ROOT_DIR, "..", "KingMaker", "CROWN", "analysis_configurations", "bbtautau", "sample_list"),
    os.path.join(ROOT_DIR, "..", "ShapeSmith", "BBTauTauAnalysis-ShapeSmith", "inventory"),
]
# in the list but deliberately not measured: the madgraph HH cross-check of the powheg signal (sample type hh2b2tau)
NOT_MEASURED = "_node_SM_"
# v7 spec 6.2: the trigger flag and offline threshold of the triggering lepton per year and channel
TRIGGER = {
    "2016": {"mt": ("trg_single_mu22", "pt_1 > 24"), "et": ("trg_single_ele25", "pt_1 > 27"),
             "tt": ("trg_double_tau35_mediumiso", "pt_1 > 40"), "ee": ("trg_single_ele25", "pt_1 > 27"),
             "em": ("trg_single_mu22", "pt_2 > 24"), "mm": ("trg_single_mu22", "pt_1 > 24")},
    "2017": {"mt": ("trg_single_mu27", "pt_1 > 29"), "et": ("trg_single_ele35", "pt_1 >= 37"),
             "tt": ("trg_double_tau40_tightiso", "pt_1 > 40"), "ee": ("trg_single_ele35", "pt_1 >= 34"),
             "em": ("trg_single_mu27", "pt_2 > 29"), "mm": ("trg_single_mu27", "pt_1 > 29")},
    "2018": {"mt": ("trg_single_mu24", "pt_1 > 25"), "et": ("trg_single_ele32", "pt_1 > 33"),
             "tt": ("trg_double_tau40_tightiso", "pt_1 > 40"), "ee": ("trg_single_ele32", "pt_1 > 34"),
             "em": ("trg_single_mu24", "pt_2 > 26"), "mm": ("trg_single_mu24", "pt_1 > 26")},
}


def year(era):
    return "2016" if era.startswith("2016") else era


def presel(era, channel):
    return load_config(os.path.join(CONFIG_DIR, era, f"preselection_{channel}.yaml"))


def efficiency(era, channel):
    return load_config(os.path.join(CONFIG_DIR, era, f"btag_efficiency_{channel}.yaml"))


@pytest.fixture(scope="module")
def database():
    with open(DATABASE) as handle:
        return json.load(handle)


@pytest.mark.parametrize("era", ERAS)
def test_common_settings(era):
    config = presel(era, "mt")
    assert config["era"] == era
    assert config["ntuple_path"].endswith("/smrun2_btag_eff_v1/CROWNRun")
    assert config["output_path"] == config["file_path"] == f"/ceph/sdaigler/CMS_FF/btag_eff_{era}"
    assert config["workdir_name"] == f"btag_eff_{era}"
    assert config["nanoAOD_version"] == "nanoAOD_v15"


@pytest.mark.parametrize("era", ERAS)
def test_samples_in_database(era, database):
    """Every nick exists in the DB with the config's era and its process as sample_type; one inclusive DY cut."""
    config = presel(era, "mt")
    expected = SM_SAMPLE_TYPES - ({"rem_htautau"} if era == "2016postVFP" else set())
    assert set(config["processes"]) == expected
    for process, entry in config["processes"].items():
        assert entry["samples"], process
        for nick in entry["samples"]:
            assert nick in database, nick
            assert database[nick]["era"] == era, nick
            assert database[nick]["sample_type"] == process, nick
        assert set(entry.get("sample_cuts", {})) <= set(entry["samples"])
    cuts = config["processes"]["dyjets"]["sample_cuts"]
    assert len(cuts) == 1
    (nick, cut), = cuts.items()
    assert nick.startswith("DYJetsToLL_M-50_") and cut == "npartons == 0"
    nicks = [n for entry in config["processes"].values() for n in entry["samples"]]
    assert len(nicks) == len(set(nicks))
    # never both ZZ -> 2nu2q samples in one era (v7 spec 3.2)
    assert sum(n.startswith(("ZZTo2Nu2Q_5f", "ZZTo2Q2Nu_")) for n in nicks) == 1


@pytest.mark.parametrize("era", ERAS)
def test_nicks_match_the_crown_mc_list(era):
    """The nicks are exactly those of sm_<era>_binned_mc.txt without the node_SM cross-check."""
    name = f"sm_{era}_binned_mc.txt"
    paths = [os.path.join(d, name) for d in SAMPLE_LIST_DIRS if d and os.path.isfile(os.path.join(d, name))]
    if not paths:
        pytest.skip(f"{name} not found in {SAMPLE_LIST_DIRS} (set CROWN_SAMPLE_LIST_DIR)")
    with open(paths[0]) as handle:
        listed = {line.strip() for line in handle if line.strip() and not line.lstrip().startswith("#")}
    not_measured = {nick for nick in listed if NOT_MEASURED in nick}
    assert len(not_measured) == 1, not_measured
    config = presel(era, "mt")
    nicks = {nick for entry in config["processes"].values() for nick in entry["samples"]}
    expected = listed - not_measured
    assert nicks == expected, f"not in {paths[0]}: {sorted(nicks - expected)}; missing: {sorted(expected - nicks)}"


@pytest.mark.parametrize("era", ERAS)
@pytest.mark.parametrize("channel", CHANNELS)
def test_selection_per_era(era, channel):
    """Trigger, threshold, lepton |eta| < 2.1 (E31) and the 2018 HEM veto (E29) after v7 spec 6.2."""
    config = presel(era, channel)
    assert config["channel"] == channel
    selection = " && ".join(config["event_selection"].values())
    flag, threshold = TRIGGER[year(era)][channel]
    assert flag in selection and threshold in selection
    assert ("hem_veto < 0.5" in selection) == (era == "2018")
    if channel in ("et", "mt") or (year(era) == "2016" and channel in ("ee", "mm")):
        assert "abs(eta_1) < 2.1" in selection
    if year(era) == "2016" and channel == "em":
        assert "abs(eta_2) < 2.1" in selection
    if year(era) == "2016" and channel == "tt":
        assert "mediumcombinediso" not in selection  # the MC takes the Iso path only (E7, D9)
    weights = " ".join(str(v) for v in config["mc_weights"].values())
    assert ("0.991" in weights) == (year(era) == "2017" and channel in ("et", "ee"))
    used_trigger_weights = set(re.findall(r"trg_wgt\w+", weights))
    for name in used_trigger_weights:
        if channel != "tt":
            assert name.replace("trg_wgt_", "trg_") in selection, name


@pytest.mark.parametrize("era", ERAS)
@pytest.mark.parametrize("channel", CHANNELS)
def test_efficiency_config(era, channel):
    config = efficiency(era, channel)
    assert config["era"] == era and config["channel"] == channel
    assert dict(config["btag_working_points"]) == WORKING_POINTS[era]
    assert config["validate_payload"] == "error"
    assert set(config["required_sample_types"]) == SM_SAMPLE_TYPES
    keys, inputs = [], []
    for name, pool in config["processes"].items():
        keys += list(pool["sample_types"])
        inputs += list(pool["input_processes"])
        for flavor in ("b", "c", "light"):
            bins = pool["flavor_bins"][flavor]
            pt, eta = list(bins["jet_pt_bins"]), list(bins["jet_eta_bins"])
            assert pt == sorted(set(pt)) and pt[0] == 20 and pt[-1] == 1000, (name, flavor, pt)
            assert eta == sorted(set(eta)) and eta[0] == 0.0 and eta[-1] == 2.4, (name, flavor, eta)
    assert sorted(keys) == sorted(SM_SAMPLE_TYPES)  # every key exactly once
    assert len(inputs) == len(set(inputs))  # disjoint pools
    produced = set(presel(era, channel)["processes"])
    assert set(inputs) == produced  # every measured process is preselected and vice versa
    pools = set(config["processes"])
    if channel in ("ee", "mm"):
        assert "higgs_tautau" not in pools and {"ggh_htautau", "vbf_htautau"} <= set(config["processes"]["vjets"]["input_processes"])
    else:
        assert pools == {"vjets", "top", "diboson", "higgs_bb", "higgs_tautau"}


def test_efficiency_config_passes_the_measurement_checks():
    """btag_efficiency.validate_process_config (duplicate keys, shared inputs, required sample types)."""
    be = pytest.importorskip("btag_efficiency")
    for era in ERAS:
        for channel in CHANNELS:
            be.validate_process_config(efficiency(era, channel))


@pytest.mark.parametrize("era", ERAS)
def test_working_points_match_btv_payload(era):
    path = BTV_PIN.format(era=era)
    if not os.path.exists(path):
        pytest.skip("cvmfs not available")
    with gzip.open(path, "rt") as handle:
        payload = json.load(handle)
    (wp_values,) = [c for c in payload["corrections"] if c["name"] == "UParTAK4_wp_values"]
    values = {item["key"]: item["value"] for item in wp_values["data"]["content"]}
    for wp, cut in WORKING_POINTS[era].items():
        assert values[wp] == pytest.approx(cut)


def _columns(expressions):
    names = set()
    for expression in expressions:
        names.update(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", expression))
    return names - {"abs", "float", "int", "double", "bool"}


@pytest.mark.parametrize("era", ERAS)
@pytest.mark.parametrize("channel", CHANNELS)
def test_expressions_compile(era, channel):
    """Every selection and weight expression is valid C++ for RDataFrame (dummy scalar columns)."""
    ROOT = pytest.importorskip("ROOT")
    config = presel(era, channel)
    expressions = [str(v) for v in config["event_selection"].values()]
    expressions += [str(v) for k, v in config["mc_weights"].items() if k not in ("generator", "lumi") and str(v)]
    rdf = ROOT.RDataFrame(1)
    for name in sorted(_columns(expressions)):
        rdf = rdf.Define(name, "1.f")
    for index, expression in enumerate(expressions):
        rdf = rdf.Define(f"_expr_{index}", f"(double)({expression})")
    assert rdf.Count().GetValue() == 1
