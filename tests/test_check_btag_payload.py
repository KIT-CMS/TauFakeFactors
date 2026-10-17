"""Tests of scripts/check_btag_payload.py on small synthetic correctionlib payloads.

The scale-factor payload mimics the layout of the BTV UParTAK4 payload (systematic -> working_point -> flavor ->
abseta -> pt; the b/c pT axis starts at 30 GeV with a constant flow value), so the consumer emulation of
btag_sf_strict.cxx is exercised without cvmfs. Needs correctionlib (LCG_108 view or conda env fake_factors):

    python3 -m pytest tests/test_check_btag_payload.py
"""
import gzip
import json
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

import check_btag_payload as chk  # noqa: E402

THRESHOLDS = {"L": 0.05, "M": 0.161, "T": 0.5405}
EFF_PT = [20.0, 30.0, 1000.0]
# per (wp, flavour) the efficiency in pT [20, 30) and [30, 1000)
EFF = {
    ("M", 5): (0.85, 0.5), ("T", 5): (0.6, 0.3),
    ("M", 4): (0.3, 0.3), ("T", 4): (0.1, 0.1),
    ("M", 0): (0.1, 0.1), ("T", 0): (0.02, 0.02),
}


def write(path, corrections):
    with gzip.open(path, "wt") as handle:
        json.dump({"schema_version": 2, "corrections": corrections}, handle)
    return str(path)


def eff_payload(path, sample_types=chk.SM_SAMPLE_TYPES, effs=EFF):
    wps = sorted({wp for wp, _ in effs})
    flavours = sorted({flavour for _, flavour in effs})

    def pt_node(wp, flavour):
        return {"nodetype": "binning", "input": "jet_pt", "edges": EFF_PT, "flow": "clamp",
                "content": list(effs[(wp, flavour)])}

    def wp_node(wp):
        return {"nodetype": "category", "input": "jet_flavor", "content": [
            {"key": flavour, "value": {"nodetype": "binning", "input": "jet_eta", "edges": [0.0, 2.4],
                                       "flow": "clamp", "content": [pt_node(wp, flavour)]}}
            for flavour in flavours if (wp, flavour) in effs]}

    data = {"nodetype": "category", "input": "sample_type", "content": [
        {"key": sample_type, "value": {"nodetype": "category", "input": "working_point",
                                       "content": [{"key": wp, "value": wp_node(wp)} for wp in wps]}}
        for sample_type in sample_types]}
    inputs = [{"name": "sample_type", "type": "string"}, {"name": "working_point", "type": "string"},
              {"name": "jet_flavor", "type": "int"}, {"name": "jet_eta", "type": "real"},
              {"name": "jet_pt", "type": "real"}]
    return write(path, [{"name": "btag_efficiency", "version": 1, "inputs": inputs,
                         "output": {"name": "efficiency", "type": "real"}, "data": data}])


def sf_payload(path, comb, comb_flow=None, light=None):
    """*comb*: {systematic: {wp: SF in pT [30, 300)}}; above 300 GeV (and, if evaluated there, below 30 GeV)
    the flow value comb_flow[systematic][wp] applies, 1.0 by default. *light*: {systematic: {wp: SF}}."""
    light = light or {"central": {"M": 1.0, "T": 1.0}, "up": {"M": 1.0, "T": 1.0}}
    inputs = [{"name": "systematic", "type": "string"}, {"name": "working_point", "type": "string"},
              {"name": "flavor", "type": "int"}, {"name": "abseta", "type": "real"}, {"name": "pt", "type": "real"}]

    def correction(name, table, flavours, pt_edges, flow):
        def leaf(systematic, wp):
            return {"nodetype": "binning", "input": "abseta", "edges": [0.0, 2.4 if name == chk.COMB_SF else 2.5],
                    "flow": "error", "content": [{
                        "nodetype": "binning", "input": "pt", "edges": pt_edges,
                        "flow": flow(systematic, wp), "content": [table[systematic][wp]]}]}
        data = {"nodetype": "category", "input": "systematic", "content": [
            {"key": systematic, "value": {"nodetype": "category", "input": "working_point", "content": [
                {"key": wp, "value": {"nodetype": "category", "input": "flavor", "content": [
                    {"key": flavour, "value": leaf(systematic, wp)} for flavour in flavours]}}
                for wp in wps]}}
            for systematic, wps in table.items()]}
        return {"name": name, "version": 1, "inputs": inputs, "output": {"name": "sf", "type": "real"}, "data": data}

    flows = comb_flow or {}
    wp_values = {"name": chk.WP_VALUES, "version": 1, "inputs": [{"name": "working_point", "type": "string"}],
                 "output": {"name": "cut", "type": "real"},
                 "data": {"nodetype": "category", "input": "working_point",
                          "content": [{"key": wp, "value": cut} for wp, cut in THRESHOLDS.items()]}}
    return write(path, [
        correction(chk.COMB_SF, comb, [4, 5], [30.0, 300.0], lambda s, wp: flows.get(s, {}).get(wp, 1.0)),
        correction(chk.LIGHT_SF, light, [0], [20.0, 1000.0], lambda s, wp: "clamp"),
        wp_values,
    ])


UNITY = {"central": {"M": 1.0, "T": 1.0}, "up": {"M": 1.0, "T": 1.0}}


def kinds(problems):
    return sorted({kind for kind, _ in problems})


# --- the consumer's jet weight -------------------------------------------------------------------------------------

def test_k1_sets_a_negative_intermediate_numerator_to_unity():
    """eff M 0.8 / T 0.6 with SF M 0.9 / T 1.3: 0.72 - 0.78 < 0 -> factor 1, flagged as K1."""
    weights = chk.jet_weights({"M": 0.8, "T": 0.6}, {"M": 0.9, "T": 1.3}, ["M", "T"])
    assert [region for region, _, _ in weights] == ["pass T", "M not T", "fails M"]
    assert weights[0] == ("pass T", 1.3, False)
    assert weights[1] == ("M not T", 1.0, True)
    assert weights[2][1] == pytest.approx((1 - 0.72) / 0.2) and not weights[2][2]


def test_zero_numerator_is_no_k1_case():
    """SF M 1.0 / T 1.5 on eff M 0.75 / T 0.5 gives num = 0 exactly: weight 0, not unity."""
    region, jet_w, unity = chk.jet_weights({"M": 0.75, "T": 0.5}, {"M": 1.0, "T": 1.5}, ["M", "T"])[1]
    assert (region, jet_w, unity) == ("M not T", 0.0, False)


def test_sf_check_counts_k1_without_a_violation(tmp_path):
    effs = {**EFF, ("M", 5): (0.8, 0.8), ("T", 5): (0.6, 0.6)}
    comb = {"central": {"M": 0.9, "T": 1.3}}
    problems, stats = chk.check_sf_weights(
        eff_payload(tmp_path / "eff.json.gz", ["ttbar"], effs),
        sf_payload(tmp_path / "sf.json.gz", comb, comb_flow={"central": {"M": 0.9, "T": 1.3}},
                   light={"central": {"M": 1.0, "T": 1.0}}),
        wps=["T", "M"],
    )
    assert problems == []
    assert stats["wps"] == ["M", "T"]  # ordered by the BTV thresholds, not as passed
    assert stats["k1_unity"] > 0


def test_sf_check_reports_a_zero_numerator(tmp_path):
    effs = {**EFF, ("M", 5): (0.75, 0.75), ("T", 5): (0.5, 0.5)}
    comb = {"central": {"M": 1.0, "T": 1.5}}
    problems, stats = chk.check_sf_weights(
        eff_payload(tmp_path / "eff.json.gz", ["ttbar"], effs),
        sf_payload(tmp_path / "sf.json.gz", comb, comb_flow={"central": {"M": 1.0, "T": 1.5}},
                   light={"central": {"M": 1.0, "T": 1.0}}),
    )
    assert problems and kinds(problems) == ["JET WEIGHT <= 0"]
    assert all("[M not T]: jet_w=0" in detail for _, detail in problems)
    assert stats["k1_unity"] == 0


def test_b_jet_below_30_gev_uses_the_first_bin_with_twice_the_deviation(tmp_path):
    """A 25 GeV b jet takes the SF at 30.001 GeV with 2x(up - central): 1 + 2 x 0.1 = 1.2, and 1.2 x 0.85 > 1.

    Evaluated at 25 GeV it would get the flow value (1.0) and pass; with the plain deviation (1.1 x 0.85 < 1)
    it would pass as well."""
    comb = {"central": {"M": 1.0, "T": 1.0}, "up": {"M": 1.1, "T": 1.1}}
    problems, stats = chk.check_sf_weights(
        eff_payload(tmp_path / "eff.json.gz", ["ttbar", "dyjets"]),
        sf_payload(tmp_path / "sf.json.gz", comb),
    )
    assert stats["bc_min_pt"] == 30.0 and 25.0 in stats["pt_probes"]
    assert kinds(problems) == ["JET WEIGHT <= 0"]
    assert sorted(detail.split(":")[0] for _, detail in problems) == [
        f"{sample_type}/up/flavor=5/eta=1.2/pt=25/[fails M]" for sample_type in ("dyjets", "ttbar")]
    # the margin report names the same cell, with eff(M) and the ceiling 1/SF(M) = 1/1.2
    jet_w, where, eff_value, ceiling = stats["fails_loosest"]["ttbar"]
    assert jet_w == pytest.approx((1 - 1.2 * 0.85) / 0.15)
    assert where == "up/flavor=5/eta=1.2/pt=25" and eff_value == 0.85 and ceiling == pytest.approx(1 / 1.2)
    lines = chk.fails_loosest_report(stats["fails_loosest"], "M", 0.05)
    assert len(lines) == 1 and lines[0].lstrip().startswith("WARNING dyjets,ttbar: jet_w=")


# --- working points --------------------------------------------------------------------------------------------------

def test_wps_are_ordered_by_the_btv_thresholds():
    assert chk.order_wps({"T", "M"}, THRESHOLDS) == ["M", "T"]
    assert chk.order_wps(["A", "B"], {"A": 0.6, "B": 0.2}) == ["B", "A"]
    assert chk.order_wps(["XT", "T", "L", "M"]) == ["L", "M", "T", "XT"]  # by name without a BTV payload
    assert chk.order_wps(["Q", "M"]) == ["M", "Q"]  # unknown names last


def test_era_requires_exactly_the_consumer_wps(tmp_path, capsys):
    """CROWN evaluates T and M only; a leftover L would be emulated as a 'fails L' region the consumer never has."""
    sf = sf_payload(tmp_path / "sf.json.gz", UNITY)
    extra = {**EFF, ("L", 5): (0.95, 0.9), ("L", 4): (0.6, 0.6), ("L", 0): (0.3, 0.3)}
    assert chk.main(["--era", "2018", "--sf-check", "--sf-payload", sf,
                     eff_payload(tmp_path / "extra.json.gz", effs=extra)]) == 1
    out = capsys.readouterr().out
    assert "WP SET" in out and "consumer WPs: T,M" in out
    assert "WPs M,T, b/c below" in out  # the emulation still uses the consumer's list, not L,M,T

    only_m = {key: value for key, value in EFF.items() if key[0] == "M"}
    problems, _ = chk.check_file(eff_payload(tmp_path / "only_m.json.gz", effs=only_m), THRESHOLDS, ("T", "M"))
    assert kinds(problems) == ["WP SET"] and len(problems) == len(chk.SM_SAMPLE_TYPES)
    problems, _ = chk.check_sf_weights(str(tmp_path / "only_m.json.gz"), sf, wps=["T", "M"])
    assert len(problems) == 1 and problems[0][0] == "WP SET" and "no T" in problems[0][1]


def test_without_era_the_payload_wps_are_emulated(tmp_path, capsys):
    sf = sf_payload(tmp_path / "sf.json.gz", UNITY)
    assert chk.main(["--sf-check", "--sf-payload", sf, eff_payload(tmp_path / "eff.json.gz")]) == 0
    assert "consumer WPs: those of the payload" in capsys.readouterr().out


# --- sample types under --era ----------------------------------------------------------------------------------------

def test_era_reports_a_missing_rem_htautau(tmp_path, capsys):
    sf = sf_payload(tmp_path / "sf.json.gz", UNITY)
    without = [name for name in chk.SM_SAMPLE_TYPES if name != "rem_htautau"]
    assert chk.main(["--era", "2016postVFP", "--sf-payload", sf, eff_payload(tmp_path / "eff.json.gz", without)]) == 1
    out = capsys.readouterr().out
    assert "MISSING SAMPLE TYPE 'rem_htautau' is not in the payload" in out
    assert "by class: MISSING SAMPLE TYPE: 1\n" in out


def test_complete_payload_is_safe_to_ship(tmp_path, capsys):
    sf = sf_payload(tmp_path / "sf.json.gz", UNITY)
    for era in chk.RUN2_ERAS:
        assert chk.main(["--era", era, "--sf-check", "--sf-payload", sf, eff_payload(tmp_path / "eff.json.gz")]) == 0
        assert "PAYLOAD IS SAFE TO SHIP" in capsys.readouterr().out
