"""Tests for the weighted b-tag efficiency accumulators (Task 17).

Two layers:

1. Pure-python unit tests of ``helper.btag_accumulators`` against a
   hand-computed, mixed-sign-weight 2-bin example (all arithmetic is worked out
   in the comments), including the ``sumw2_fail`` identity and the empty-bin
   guard.
2. An integration test that writes a small synthetic ROOT file (Run-3 lowercase
   ``jet_*_vec`` branch contract + a scalar ``weight`` including negatives), runs
   ``btag_efficiency.calculate_efficiency_histograms`` end to end, and asserts the
   extracted accumulators / weighted efficiencies match a numpy reference
   computed directly from the generated events, while the unweighted diagnostic
   reproduces the old binomial behaviour. Skipped when ROOT is unavailable.
"""

import math
import os
import unittest

import numpy as np

from helper.btag_accumulators import (
    WeightedEfficiency,
    effective_population,
    unweighted_efficiency,
    weighted_efficiency,
)

try:
    import ROOT  # noqa: F401

    HAS_ROOT = True
except Exception:  # pragma: no cover - environment dependent
    HAS_ROOT = False


class TestWeightedEfficiencyPurePython(unittest.TestCase):
    """Hand-computed mixed-sign-weight arithmetic for the accumulator math.

    Two pt bins, one eta bin (arrays are shape (2, 1), indexed [i_pt][i_eta]).

    Bin 0 -- all-positive weights (well behaved):
        total jets: weights [2, 3, 5]        -> sumw_total  = 10, sumw2_total  = 4+9+25 = 38
        passing jets: weights [2, 3]         -> sumw_pass   = 5,  sumw2_pass   = 4+9    = 13
        (the failing jet has weight 5)       -> sumw2_fail  = 38 - 13 = 25 = 5^2  (identity)
        eff = 5 / 10                         = 0.5
        var = ((1-0.5)^2*13 + 0.5^2*25)/10^2 = (0.25*13 + 0.25*25)/100
                                             = (3.25 + 6.25)/100 = 9.5/100 = 0.095
        n_eff = 10^2 / 38                    = 100/38 = 2.631578947368421
        raw: 2/3 jets pass -> unweighted eff = 0.666..., <= 1

    Bin 1 -- mixed sign (a negative-weight failing jet drives eff > 1):
        total jets: weights [10, -4]         -> sumw_total  = 6,   sumw2_total  = 100+16 = 116
        passing jets: weights [10]           -> sumw_pass   = 10,  sumw2_pass   = 100
        (failing jet has weight -4)          -> sumw2_fail  = 116 - 100 = 16 = (-4)^2  (identity)
        eff = 10 / 6                         = 1.6666666666666667  (> 1, NOT clamped)
        var = ((1-5/3)^2*100 + (5/3)^2*16)/6^2
            = ((4/9)*100 + (25/9)*16)/36 = ((400+400)/9)/36 = (800/9)/36
            = 800/324                        = 2.4691358024691357
        n_eff = 6^2 / 116                    = 36/116 = 0.3103448275862069
        raw: 1/2 jets pass -> unweighted eff = 0.5, still <= 1
    """

    def setUp(self) -> None:
        self.acc = WeightedEfficiency(
            sumw_total=[[10.0], [6.0]],
            sumw2_total=[[38.0], [116.0]],
            sumw_pass=[[5.0], [10.0]],
            sumw2_pass=[[13.0], [100.0]],
            raw_total=[[3.0], [2.0]],
            raw_pass=[[2.0], [1.0]],
        )

    def test_weighted_efficiency_values(self) -> None:
        eff, var = weighted_efficiency(self.acc)
        np.testing.assert_allclose(eff, [[0.5], [10.0 / 6.0]])
        np.testing.assert_allclose(var, [[0.095], [800.0 / 324.0]])

    def test_negative_weight_bin_exceeds_one_while_unweighted_does_not(self) -> None:
        # The mixed-sign bin (bin 1) has a weighted efficiency > 1 ...
        eff, _ = weighted_efficiency(self.acc)
        self.assertGreater(eff[1][0], 1.0)
        self.assertAlmostEqual(eff[1][0], 1.6666666666666667)
        # ... while the unweighted, clamped diagnostic stays <= 1.
        unw_eff, _ = unweighted_efficiency(self.acc)
        self.assertLessEqual(unw_eff[1][0], 1.0)
        self.assertAlmostEqual(unw_eff[1][0], 0.5)

    def test_effective_population(self) -> None:
        n_eff = effective_population(self.acc)
        np.testing.assert_allclose(n_eff, [[100.0 / 38.0], [36.0 / 116.0]])

    def test_unweighted_diagnostic_values(self) -> None:
        eff, unc = unweighted_efficiency(self.acc)
        np.testing.assert_allclose(eff, [[2.0 / 3.0], [0.5]])
        expected_unc = [
            [math.sqrt((2.0 / 3.0) * (1.0 - 2.0 / 3.0) / 3.0)],
            [math.sqrt(0.5 * 0.5 / 2.0)],
        ]
        np.testing.assert_allclose(unc, expected_unc)

    def test_sumw2_fail_identity(self) -> None:
        # sumw2_fail = sumw2_total - sumw2_pass must equal the sum of squared
        # weights of the failing jets: 5^2 = 25 in bin 0, (-4)^2 = 16 in bin 1.
        sumw2_fail = self.acc.sumw2_total - self.acc.sumw2_pass
        np.testing.assert_allclose(sumw2_fail, [[25.0], [16.0]])

    def test_empty_bin_guards(self) -> None:
        empty = WeightedEfficiency(
            sumw_total=[[0.0]],
            sumw2_total=[[0.0]],
            sumw_pass=[[0.0]],
            sumw2_pass=[[0.0]],
            raw_total=[[0.0]],
            raw_pass=[[0.0]],
        )
        eff, var = weighted_efficiency(empty)
        self.assertTrue(np.isnan(eff[0][0]))
        self.assertTrue(np.isnan(var[0][0]))
        self.assertTrue(np.isnan(effective_population(empty)[0][0]))
        # The unweighted diagnostic keeps the old empty-bin behaviour: 0.0/0.0.
        unw_eff, unw_unc = unweighted_efficiency(empty)
        self.assertEqual(unw_eff[0][0], 0.0)
        self.assertEqual(unw_unc[0][0], 0.0)


# ---------------------------------------------------------------------------
# Integration test with a synthetic ROOT file
# ---------------------------------------------------------------------------

# One WP, three flavours, 2 pt bins x 1 eta bin (Run-3-style config view).
_WP_NAME = "M"
_WP_CUT = 0.5
_FLAVORS = {"b": 5, "c": 4, "light": 0}
_PT_BINS = [20.0, 100.0, 1000.0]
_ETA_BINS = [0.0, 2.5]


def _make_test_config() -> dict:
    return {
        "jet_pt_column": "jet_pt_vec",
        "jet_eta_column": "jet_eta_vec",
        "jet_flavor_column": "jet_hadronflavour_vec",
        "jet_btag_column": "jet_btag_value_vec",
        "weight_column": "weight",
        "jet_selection": "",
        "btag_working_points": {_WP_NAME: _WP_CUT},
        "jet_flavor_categories": dict(_FLAVORS),
        "jet_pt_bins": list(_PT_BINS),
        "jet_eta_bins": list(_ETA_BINS),
        "tree": "ntuple",
    }


def _generate_events(n_events: int, seed: int = 1234):
    """Return per-event jagged jet arrays + scalar weights (numpy reference)."""
    rng = np.random.default_rng(seed)
    weight_choices = np.array([1.0, 2.0, 0.5, -1.0, -0.5, 3.0])

    events = []
    for _ in range(n_events):
        n_jets = int(rng.integers(1, 5))  # 1..4 jets
        # pt strictly inside [20, 1000) and |eta| strictly inside [0, 2.5)
        # so every jet lands in a real (non-overflow) histogram bin.
        pt = rng.uniform(20.0, 999.0, size=n_jets).astype(np.float32)
        eta = rng.uniform(-2.4, 2.4, size=n_jets).astype(np.float32)
        flavor = rng.choice([0, 4, 5], size=n_jets).astype(np.int32)
        btag = rng.uniform(0.0, 1.0, size=n_jets).astype(np.float32)
        weight = float(rng.choice(weight_choices))
        events.append((pt, eta, flavor, btag, weight))
    return events


def _write_root_file(path: str, events, include_weight: bool = True) -> None:
    import ROOT

    fout = ROOT.TFile(path, "RECREATE")
    tree = ROOT.TTree("ntuple", "ntuple")

    pt_vec = ROOT.std.vector("float")()
    eta_vec = ROOT.std.vector("float")()
    flav_vec = ROOT.std.vector("int")()
    btag_vec = ROOT.std.vector("float")()
    tree.Branch("jet_pt_vec", pt_vec)
    tree.Branch("jet_eta_vec", eta_vec)
    tree.Branch("jet_hadronflavour_vec", flav_vec)
    tree.Branch("jet_btag_value_vec", btag_vec)

    import array as _array

    weight_holder = _array.array("f", [0.0])
    if include_weight:
        tree.Branch("weight", weight_holder, "weight/F")

    for pt, eta, flavor, btag, weight in events:
        pt_vec.clear()
        eta_vec.clear()
        flav_vec.clear()
        btag_vec.clear()
        for value in pt:
            pt_vec.push_back(float(value))
        for value in eta:
            eta_vec.push_back(float(value))
        for value in flavor:
            flav_vec.push_back(int(value))
        for value in btag:
            btag_vec.push_back(float(value))
        weight_holder[0] = float(weight)
        tree.Fill()

    tree.Write()
    fout.Close()


def _numpy_reference(events):
    """Compute reference sumw/sumw2/raw arrays per flavour, shape (n_pt, n_eta)."""
    n_pt = len(_PT_BINS) - 1
    n_eta = len(_ETA_BINS) - 1
    pt_edges = np.array(_PT_BINS)
    eta_edges = np.array(_ETA_BINS)

    ref = {
        flavor: {
            "sumw_total": np.zeros((n_pt, n_eta)),
            "sumw2_total": np.zeros((n_pt, n_eta)),
            "sumw_pass": np.zeros((n_pt, n_eta)),
            "sumw2_pass": np.zeros((n_pt, n_eta)),
            "raw_total": np.zeros((n_pt, n_eta)),
            "raw_pass": np.zeros((n_pt, n_eta)),
        }
        for flavor in _FLAVORS
    }

    for pt, eta, flavor, btag, weight in events:
        for j in range(len(pt)):
            # np.digitize with the bin edges; subtract 1 for 0-based index.
            i_pt = int(np.digitize(pt[j], pt_edges)) - 1
            i_eta = int(np.digitize(abs(eta[j]), eta_edges)) - 1
            if not (0 <= i_pt < n_pt and 0 <= i_eta < n_eta):
                continue
            for flavor_name, flavor_id in _FLAVORS.items():
                if int(flavor[j]) != flavor_id:
                    continue
                r = ref[flavor_name]
                r["sumw_total"][i_pt][i_eta] += weight
                r["sumw2_total"][i_pt][i_eta] += weight * weight
                r["raw_total"][i_pt][i_eta] += 1.0
                if btag[j] >= _WP_CUT:
                    r["sumw_pass"][i_pt][i_eta] += weight
                    r["sumw2_pass"][i_pt][i_eta] += weight * weight
                    r["raw_pass"][i_pt][i_eta] += 1.0
    return ref


@unittest.skipUnless(HAS_ROOT, "ROOT/RDataFrame not importable")
class TestWeightedEfficiencyIntegration(unittest.TestCase):
    """End-to-end: synthetic ROOT file -> calculate_efficiency_histograms."""

    @classmethod
    def setUpClass(cls) -> None:
        import ROOT

        ROOT.gROOT.SetBatch(True)
        cls.tmpdir = os.environ.get("TMPDIR", "/tmp")
        cls.events = _generate_events(200, seed=20240723)
        cls.root_path = os.path.join(cls.tmpdir, "btag_weighted_eff_test.root")
        _write_root_file(cls.root_path, cls.events, include_weight=True)
        cls.ref = _numpy_reference(cls.events)

    def _run_calculator(self, root_path: str):
        import ROOT

        import btag_efficiency

        chain = ROOT.TChain("ntuple")
        chain.Add(root_path)
        config = _make_test_config()
        _histos, accumulators = btag_efficiency.calculate_efficiency_histograms(
            chain, config, "test_proc", list(_PT_BINS), list(_ETA_BINS)
        )
        return accumulators

    def test_accumulators_match_numpy_reference(self) -> None:
        accumulators = self._run_calculator(self.root_path)
        for flavor in _FLAVORS:
            acc = accumulators[_WP_NAME][flavor]
            ref = self.ref[flavor]
            for key in (
                "sumw_total",
                "sumw2_total",
                "sumw_pass",
                "sumw2_pass",
                "raw_total",
                "raw_pass",
            ):
                np.testing.assert_allclose(
                    getattr(acc, key),
                    ref[key],
                    rtol=1e-5,
                    atol=1e-4,
                    err_msg=f"{key} mismatch for flavor {flavor}",
                )

    def test_weighted_efficiency_matches_reference(self) -> None:
        accumulators = self._run_calculator(self.root_path)
        for flavor in _FLAVORS:
            acc = accumulators[_WP_NAME][flavor]
            eff, var = weighted_efficiency(acc)

            ref = self.ref[flavor]
            with np.errstate(divide="ignore", invalid="ignore"):
                ref_eff = np.where(
                    ref["sumw_total"] != 0.0,
                    ref["sumw_pass"] / ref["sumw_total"],
                    np.nan,
                )
                ref_fail = ref["sumw2_total"] - ref["sumw2_pass"]
                ref_var = np.where(
                    ref["sumw_total"] != 0.0,
                    ((1.0 - ref_eff) ** 2 * ref["sumw2_pass"] + ref_eff**2 * ref_fail)
                    / ref["sumw_total"] ** 2,
                    np.nan,
                )
            np.testing.assert_allclose(eff, ref_eff, rtol=1e-5, atol=1e-6)
            np.testing.assert_allclose(var, ref_var, rtol=1e-5, atol=1e-6)

    def test_unweighted_diagnostic_matches_old_behavior(self) -> None:
        accumulators = self._run_calculator(self.root_path)
        for flavor in _FLAVORS:
            acc = accumulators[_WP_NAME][flavor]
            eff, unc = unweighted_efficiency(acc)

            ref = self.ref[flavor]
            # Reproduce the exact previous binomial formula bin-by-bin.
            n_pt, n_eta = ref["raw_total"].shape
            for i_pt in range(n_pt):
                for i_eta in range(n_eta):
                    total = ref["raw_total"][i_pt][i_eta]
                    passed = ref["raw_pass"][i_pt][i_eta]
                    if total > 0.0:
                        exp_eff = min(passed / total, 1.0)
                        exp_unc = math.sqrt(exp_eff * (1.0 - exp_eff) / total)
                    else:
                        exp_eff = 0.0
                        exp_unc = 0.0
                    self.assertAlmostEqual(eff[i_pt][i_eta], exp_eff)
                    self.assertAlmostEqual(unc[i_pt][i_eta], exp_unc)

    def test_missing_weight_column_raises(self) -> None:
        import ROOT

        import btag_efficiency

        no_weight_path = os.path.join(self.tmpdir, "btag_no_weight_test.root")
        _write_root_file(no_weight_path, self.events, include_weight=False)

        chain = ROOT.TChain("ntuple")
        chain.Add(no_weight_path)
        config = _make_test_config()
        with self.assertRaises(RuntimeError) as ctx:
            btag_efficiency.calculate_efficiency_histograms(
                chain, config, "test_proc", list(_PT_BINS), list(_ETA_BINS)
            )
        self.assertIn("weight", str(ctx.exception))


if __name__ == "__main__":
    unittest.main()
