"""Tests for the strict b-tag validation gates, bin merging and atomic install (Task 18).

All tests here are pure Python (synthetic ``WeightedEfficiency`` accumulators,
filesystem fixtures, correctionlib round-trip) and do not require ROOT, except
where explicitly noted / skipped.

The scenarios mirror the Task-18 spec:

- an exact-subset (raw) nesting violation is *immediately fatal*;
- a signed-weight (weighted) nesting violation is *not* fatal and is resolved
  by one deterministic pt-bin merge;
- a violation surviving the coarsest 1x1 binning is fatal;
- a missing process is fatal (tightened calculator code path);
- a correctionlib round-trip mismatch is fatal;
- ``atomic_install`` refuses to publish unless the report status is ``passed``
  and all required channels are present, and it is strictly all-or-nothing.
"""

import hashlib
import json
import os
import types
import unittest

import numpy as np

from helper.btag_accumulators import WeightedEfficiency
from helper.btag_validation import (
    BinnedAccumulators,
    atomic_install,
    check_wp_nesting,
    evaluate_gates,
    get_validation_thresholds,
    merge_bins,
    raw_flow_counts,
    roundtrip_check,
    validate_input_branches,
    write_validation_report,
)


def _acc(sumw_total, sumw2_total, sumw_pass, sumw2_pass, raw_total, raw_pass):
    """Build a WeightedEfficiency from column lists (each entry is one pt bin)."""
    return WeightedEfficiency(
        sumw_total=sumw_total,
        sumw2_total=sumw2_total,
        sumw_pass=sumw_pass,
        sumw2_pass=sumw2_pass,
        raw_total=raw_total,
        raw_pass=raw_pass,
    )


# ---------------------------------------------------------------------------
# validate_input_branches
# ---------------------------------------------------------------------------
class TestValidateInputBranches(unittest.TestCase):
    REQUIRED = [
        "jet_pt_vec",
        "jet_eta_vec",
        "jet_hadronflavour_vec",
        "jet_btag_value_vec",
        "weight",
    ]

    def test_all_present_single_tree(self):
        # Returns None (no raise) when every required column is present.
        self.assertIsNone(
            validate_input_branches(list(self.REQUIRED) + ["extra"], self.REQUIRED)
        )

    def test_missing_single_tree_raises_listing_names(self):
        cols = ["jet_pt_vec", "jet_eta_vec"]
        with self.assertRaises(RuntimeError) as ctx:
            validate_input_branches(cols, self.REQUIRED)
        msg = str(ctx.exception)
        self.assertIn("weight", msg)
        self.assertIn("jet_btag_value_vec", msg)

    def test_missing_per_file_mapping_lists_each_file(self):
        columns_by_file = {
            "/data/a.root": list(self.REQUIRED),  # complete
            "/data/b.root": ["jet_pt_vec"],  # missing four
        }
        with self.assertRaises(RuntimeError) as ctx:
            validate_input_branches(columns_by_file, self.REQUIRED)
        msg = str(ctx.exception)
        self.assertIn("/data/b.root", msg)
        self.assertNotIn("/data/a.root", msg)  # complete file not reported
        self.assertIn("weight", msg)


# ---------------------------------------------------------------------------
# check_wp_nesting
# ---------------------------------------------------------------------------
class TestCheckWpNesting(unittest.TestCase):
    def test_raw_exact_subset_violation_is_flagged_fatal(self):
        # 1 pt x 1 eta, two WPs ordered loose -> tight. The tighter WP has MORE
        # raw passes than the looser one -> impossible for nested pass sets.
        loose = _acc([[10.0]], [[10.0]], [[3.0]], [[3.0]], [[10.0]], [[3.0]])
        tight = _acc([[10.0]], [[10.0]], [[5.0]], [[5.0]], [[10.0]], [[5.0]])
        result = check_wp_nesting({"L": loose, "T": tight})
        self.assertFalse(result.raw_ok)
        self.assertTrue(result.raw_violations)
        v = result.raw_violations[0]
        self.assertEqual((v["loose"], v["tight"]), ("L", "T"))

    def test_raw_nesting_ok_but_weighted_range_violation_not_fatal(self):
        # Raw counts nest correctly (loose 5 >= tight 2), but a negative-weight
        # failing jet drives the weighted efficiency above 1 -> weighted only.
        loose = _acc([[6.0]], [[116.0]], [[10.0]], [[100.0]], [[3.0]], [[3.0]])
        tight = _acc([[6.0]], [[116.0]], [[10.0]], [[100.0]], [[3.0]], [[2.0]])
        result = check_wp_nesting({"L": loose, "T": tight})
        self.assertTrue(result.raw_ok)  # raw subset intact
        self.assertFalse(result.weighted_ok)  # weighted eff out of [0, 1]
        self.assertTrue(result.weighted_violations)


# ---------------------------------------------------------------------------
# merge_bins
# ---------------------------------------------------------------------------
class TestMergeBins(unittest.TestCase):
    THRESHOLDS = {"min_effective_denominator": 1.0, "max_stat_uncertainty": 100.0}

    def _signed_weight_case(self):
        # 2 pt bins x 1 eta, WPs L (loose) and T (tight). Bin 0 is well behaved;
        # bin 1 has a negative L-only weight so weighted eff is non-nesting
        # (eff_T > eff_L) there while raw counts still nest. Numbers worked out
        # so that merging bin 1 into bin 0 restores weighted nesting.
        #
        # bin 0: total sumw 17 (sumw2 129); L pass 15 (sumw2 125); T pass 5 (25)
        # bin 1: total sumw 14 (sumw2 180); L pass  6 (sumw2 116); T pass 10 (100)
        loose = _acc(
            sumw_total=[[17.0], [14.0]],
            sumw2_total=[[129.0], [180.0]],
            sumw_pass=[[15.0], [6.0]],
            sumw2_pass=[[125.0], [116.0]],
            raw_total=[[3.0], [3.0]],
            raw_pass=[[2.0], [2.0]],
        )
        tight = _acc(
            sumw_total=[[17.0], [14.0]],
            sumw2_total=[[129.0], [180.0]],
            sumw_pass=[[5.0], [10.0]],
            sumw2_pass=[[25.0], [100.0]],
            raw_total=[[3.0], [3.0]],
            raw_pass=[[1.0], [1.0]],
        )
        return BinnedAccumulators(
            accumulators={"L": loose, "T": tight},
            pt_bins=[20.0, 100.0, 1000.0],
            eta_bins=[0.0, 2.5],
        )

    def test_signed_weight_nesting_resolved_by_one_pt_merge(self):
        binned = self._signed_weight_case()

        # Pre-merge: bin 1 is a weighted non-nesting failure.
        pre = evaluate_gates(binned, self.THRESHOLDS)
        self.assertTrue(pre.any_fail)
        self.assertTrue(bool(pre.fail_mask[1][0]))

        merged, history = merge_bins(binned, self.THRESHOLDS)

        # Exactly one pt merge was needed.
        self.assertEqual(len(history), 1)
        self.assertEqual(history[0]["op"], "merge_pt")
        self.assertEqual(merged.pt_bins, [20.0, 1000.0])

        # Arrays were summed (never efficiencies): combined total sumw = 31.
        np.testing.assert_allclose(merged.accumulators["L"].sumw_total, [[31.0]])
        np.testing.assert_allclose(merged.accumulators["L"].sumw_pass, [[21.0]])
        np.testing.assert_allclose(merged.accumulators["T"].sumw_pass, [[15.0]])

        # Re-running every gate post-merge: no bin fails, weighted nesting holds.
        post = evaluate_gates(merged, self.THRESHOLDS)
        self.assertFalse(post.any_fail)
        nres = check_wp_nesting(merged.accumulators)
        self.assertTrue(nres.raw_ok)
        self.assertTrue(nres.weighted_ok)

    def test_violation_surviving_coarsest_binning_is_fatal(self):
        # 1 pt x 1 eta with effective denominator 1 < min 100 -> cannot be
        # merged any coarser -> fatal.
        acc = _acc([[5.0]], [[25.0]], [[2.0]], [[10.0]], [[5.0]], [[2.0]])
        binned = BinnedAccumulators(
            accumulators={"M": acc}, pt_bins=[20.0, 1000.0], eta_bins=[0.0, 2.5]
        )
        with self.assertRaises(RuntimeError) as ctx:
            merge_bins(
                binned,
                {"min_effective_denominator": 100.0, "max_stat_uncertainty": 0.25},
            )
        self.assertIn("1", str(ctx.exception))  # references the 1x1 coarsest grid

    def test_bin_failing_only_for_tightest_wp_merges_all_wps(self):
        # 2 pt bins x 1 eta, WPs L and XXT sharing totals. Only XXT fails (its
        # weighted eff goes negative in bin 1); L stays in range. The merge must
        # collapse the pt bin for BOTH WPs (merged-together semantics).
        loose = _acc(
            sumw_total=[[20.0], [30.0]],
            sumw2_total=[[40.0], [200.0]],
            sumw_pass=[[10.0], [9.0]],
            sumw2_pass=[[20.0], [125.0]],
            raw_total=[[5.0], [6.0]],
            raw_pass=[[3.0], [3.0]],
        )
        xxt = _acc(
            sumw_total=[[20.0], [30.0]],
            sumw2_total=[[40.0], [200.0]],
            sumw_pass=[[6.0], [-3.0]],  # negative -> eff_XXT < 0 in bin 1
            sumw2_pass=[[12.0], [25.0]],
            raw_total=[[5.0], [6.0]],
            raw_pass=[[2.0], [2.0]],
        )
        binned = BinnedAccumulators(
            accumulators={"L": loose, "XXT": xxt},
            pt_bins=[20.0, 100.0, 1000.0],
            eta_bins=[0.0, 2.5],
        )
        thresholds = {"min_effective_denominator": 1.0, "max_stat_uncertainty": 100.0}

        # Pre-merge the failure is confined to the XXT working point.
        pre = evaluate_gates(binned, thresholds)
        self.assertTrue(bool(pre.fail_mask[1][0]))
        self.assertIn("XXT", pre.failing_working_points[(1, 0)])
        self.assertNotIn("L", pre.failing_working_points[(1, 0)])

        merged, history = merge_bins(binned, thresholds)
        self.assertEqual(len(history), 1)
        self.assertEqual(history[0]["op"], "merge_pt")

        # BOTH WPs collapsed to a single pt bin with summed arrays.
        self.assertEqual(merged.accumulators["L"].shape, (1, 1))
        self.assertEqual(merged.accumulators["XXT"].shape, (1, 1))
        np.testing.assert_allclose(merged.accumulators["L"].sumw_total, [[50.0]])
        np.testing.assert_allclose(merged.accumulators["XXT"].sumw_total, [[50.0]])
        np.testing.assert_allclose(merged.accumulators["XXT"].sumw_pass, [[3.0]])
        np.testing.assert_allclose(merged.accumulators["L"].sumw_pass, [[19.0]])

    def test_no_failing_bins_returns_empty_history(self):
        acc = _acc([[500.0]], [[500.0]], [[250.0]], [[250.0]], [[500.0]], [[250.0]])
        binned = BinnedAccumulators(
            accumulators={"M": acc}, pt_bins=[20.0, 1000.0], eta_bins=[0.0, 2.5]
        )
        merged, history = merge_bins(
            binned, {"min_effective_denominator": 100.0, "max_stat_uncertainty": 0.25}
        )
        self.assertEqual(history, [])
        self.assertEqual(merged.pt_bins, [20.0, 1000.0])


# ---------------------------------------------------------------------------
# thresholds
# ---------------------------------------------------------------------------
class TestThresholds(unittest.TestCase):
    def test_defaults_when_absent(self):
        th = get_validation_thresholds({})
        self.assertEqual(th["min_effective_denominator"], 100.0)
        self.assertEqual(th["max_stat_uncertainty"], 0.25)

    def test_explicit_config_values_win(self):
        th = get_validation_thresholds(
            {
                "validation": {
                    "min_effective_denominator": 250,
                    "max_stat_uncertainty": 0.1,
                }
            }
        )
        self.assertEqual(th["min_effective_denominator"], 250.0)
        self.assertEqual(th["max_stat_uncertainty"], 0.1)


# ---------------------------------------------------------------------------
# correctionlib round-trip gate
# ---------------------------------------------------------------------------
class TestRoundTrip(unittest.TestCase):
    def setUp(self):
        import btag_efficiency

        self.btag_efficiency = btag_efficiency
        self.tmpdir = os.environ.get("TMPDIR", "/tmp")
        self.wps = {"L": 0.1, "M": 0.5}
        self.flavors = {"b": 5, "c": 4, "light": 0}
        # Per (process, flavor) binning; here shared for simplicity.
        pt = [20.0, 100.0, 1000.0]
        eta = [0.0, 2.5]
        self.flavor_bins = {"b": (pt, eta), "c": (pt, eta), "light": (pt, eta)}
        self.all_flavor_bins = {"ttbar": self.flavor_bins}

        def grid(scale):
            return [[0.1 * scale], [0.2 * scale]]

        self.all_efficiencies = {
            "ttbar": {
                wp: {flav: grid(i + 1) for i, flav in enumerate(self.flavors)}
                for wp in self.wps
            }
        }
        self.config = {
            "btag_working_points": dict(self.wps),
            "jet_flavor_categories": dict(self.flavors),
            "jet_pt_bins": pt,
            "jet_eta_bins": eta,
            "era": "test",
        }
        self.json_path = os.path.join(self.tmpdir, "btag_rt_test.json")
        self.btag_efficiency.build_correctionlib_json(
            self.all_efficiencies, self.all_flavor_bins, self.config, self.tmpdir
        )
        # build_correctionlib_json writes <output>/btag_efficiency.json
        self.json_path = os.path.join(self.tmpdir, "btag_efficiency.json")

    def test_roundtrip_matches_before_tamper(self):
        ok, n_checked, mismatches = roundtrip_check(
            self.json_path,
            self.all_efficiencies,
            self.all_flavor_bins,
            self.wps,
            self.flavors,
        )
        self.assertTrue(ok, msg=f"unexpected mismatches: {mismatches}")
        self.assertGreater(n_checked, 0)

    def test_roundtrip_mismatch_is_detected_after_tamper(self):
        with open(self.json_path) as f:
            payload = json.load(f)
        # Tamper: bump the first efficiency float we can find in the tree.
        self._bump_first_float(payload)
        with open(self.json_path, "w") as f:
            json.dump(payload, f)

        ok, _n, mismatches = roundtrip_check(
            self.json_path,
            self.all_efficiencies,
            self.all_flavor_bins,
            self.wps,
            self.flavors,
        )
        self.assertFalse(ok)
        self.assertTrue(mismatches)

    def test_roundtrip_tolerates_ulp_serialization_roundoff(self):
        """A ~1-ULP difference between the stored efficiency and the value
        correctionlib returns after JSON (de)serialization must NOT be flagged.

        This reproduces the real failure seen in the Task-24 integration run:
        efficiencies such as 0.9210030090094377 come back from correctionlib as
        0.9210030090094375 (last-bit round-off), which an exact ``==`` compare
        wrongly rejected. Here we perturb every expected efficiency by one ULP
        and require the round-trip gate to still pass.
        """
        import copy
        import math

        perturbed = copy.deepcopy(self.all_efficiencies)
        for wp in perturbed["ttbar"]:
            for flav in perturbed["ttbar"][wp]:
                grid = perturbed["ttbar"][wp][flav]
                for r in range(len(grid)):
                    for c in range(len(grid[r])):
                        grid[r][c] = math.nextafter(grid[r][c], math.inf)
        ok, _n, mismatches = roundtrip_check(
            self.json_path,
            perturbed,
            self.all_flavor_bins,
            self.wps,
            self.flavors,
        )
        self.assertTrue(ok, msg=f"ULP round-off must not fail: {mismatches[:3]}")

    def test_roundtrip_still_rejects_physically_meaningful_difference(self):
        """The tolerance must not be so loose that a real (0.05) efficiency
        error slips through -- the gate still has teeth."""
        import copy

        perturbed = copy.deepcopy(self.all_efficiencies)
        perturbed["ttbar"]["L"]["b"][0][0] += 0.05
        ok, _n, mismatches = roundtrip_check(
            self.json_path,
            perturbed,
            self.all_flavor_bins,
            self.wps,
            self.flavors,
        )
        self.assertFalse(ok)
        self.assertTrue(mismatches)

    def _bump_first_float(self, node):
        # Walk the correctionlib content tree and mutate the first leaf float.
        corr = node["corrections"][0]["data"]

        def walk(n):
            if isinstance(n, dict):
                content = n.get("content")
                if isinstance(content, list):
                    for i, c in enumerate(content):
                        if isinstance(c, (int, float)) and not isinstance(c, bool):
                            content[i] = float(c) + 0.5
                            return True
                        if walk(c):
                            return True
                if "value" in n and walk(n["value"]):
                    return True
            elif isinstance(n, list):
                for c in n:
                    if walk(c):
                        return True
            return False

        self.assertTrue(walk(corr), "could not find a float leaf to tamper")


# ---------------------------------------------------------------------------
# missing process is fatal (tightened calculator code path)
# ---------------------------------------------------------------------------
class TestMissingProcessFatal(unittest.TestCase):
    def test_missing_process_files_raises(self):
        import btag_efficiency

        tmpdir = os.environ.get("TMPDIR", "/tmp")
        out_base = os.path.join(tmpdir, "btag_missing_proc_test")

        config = {
            "file_path": "/nonexistent",
            "output_base": out_base,
            "workdir_name": "wd",
            "era": "test",
            "tree": "ntuple",
            "channel": "mt",
            "jet_pt_column": "jet_pt_vec",
            "jet_eta_column": "jet_eta_vec",
            "jet_flavor_column": "jet_hadronflavour_vec",
            "jet_btag_column": "jet_btag_value_vec",
            "weight_column": "weight",
            "jet_selection": "",
            "btag_working_points": {"M": 0.5},
            "jet_flavor_categories": {"b": 5},
            "jet_pt_bins": [20.0, 1000.0],
            "jet_eta_bins": [0.0, 2.5],
            "processes": {"ttbar": {"sample_type": "ttbar"}},
        }
        args = types.SimpleNamespace(workers=1, threads=1)

        # Force "no files found" without touching the filesystem.
        original = btag_efficiency.get_process_files
        btag_efficiency.get_process_files = lambda cfg, proc: []
        try:
            with self.assertRaises(RuntimeError) as ctx:
                btag_efficiency.run_for_channel(config, args)
            self.assertIn("ttbar", str(ctx.exception))
        finally:
            btag_efficiency.get_process_files = original


# ---------------------------------------------------------------------------
# atomic_install
# ---------------------------------------------------------------------------
def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


class TestAtomicInstall(unittest.TestCase):
    def setUp(self):
        self.tmpdir = os.path.join(
            os.environ.get("TMPDIR", "/tmp"), f"btag_install_{os.getpid()}"
        )
        os.makedirs(self.tmpdir, exist_ok=True)
        self.src_dir = os.path.join(self.tmpdir, "staging")
        os.makedirs(self.src_dir, exist_ok=True)
        self.payload = []
        for name, text in (
            ("btag_efficiency.json", '{"a": 1}'),
            ("btag_efficiency.json.gz", "binary-ish"),
        ):
            p = os.path.join(self.src_dir, name)
            with open(p, "w") as f:
                f.write(text)
            self.payload.append(p)
        self.checksums = {os.path.basename(p): _sha256(p) for p in self.payload}
        self.target = os.path.join(self.tmpdir, "published")

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _provenance(self, **overrides):
        prov = {
            "status": "passed",
            "channels_present": ["et", "mt", "tt"],
            "required_channels": ["et", "mt", "tt"],
            "checksums": dict(self.checksums),
        }
        prov.update(overrides)
        return prov

    def test_installs_when_passed_and_channels_present(self):
        atomic_install(self.payload, self._provenance(), self.target)
        self.assertTrue(os.path.isdir(self.target))
        self.assertTrue(
            os.path.isfile(os.path.join(self.target, "btag_efficiency.json"))
        )
        self.assertTrue(
            os.path.isfile(os.path.join(self.target, "btag_efficiency.json.gz"))
        )

    def test_refuses_when_status_not_passed(self):
        with self.assertRaises(RuntimeError):
            atomic_install(self.payload, self._provenance(status="failed"), self.target)
        self.assertFalse(os.path.exists(self.target))
        # No stray temp dir left behind.
        self.assertEqual(
            [d for d in os.listdir(self.tmpdir) if d.startswith("published.tmp")], []
        )

    def test_refuses_when_a_required_channel_is_missing(self):
        with self.assertRaises(RuntimeError):
            atomic_install(
                self.payload,
                self._provenance(channels_present=["et", "mt"]),
                self.target,
            )
        self.assertFalse(os.path.exists(self.target))

    def test_all_or_nothing_on_checksum_error_leaves_target_untouched(self):
        # Pre-existing published payload must survive a failed install verbatim.
        os.makedirs(self.target)
        sentinel = os.path.join(self.target, "OLD")
        with open(sentinel, "w") as f:
            f.write("previous")

        bad = dict(self.checksums)
        bad["btag_efficiency.json"] = "0" * 64  # wrong checksum -> corruption
        with self.assertRaises(RuntimeError):
            atomic_install(self.payload, self._provenance(checksums=bad), self.target)

        # Old payload untouched; the new files never landed.
        self.assertTrue(os.path.isfile(sentinel))
        self.assertFalse(
            os.path.isfile(os.path.join(self.target, "btag_efficiency.json"))
        )
        self.assertEqual(
            [d for d in os.listdir(self.tmpdir) if d.startswith("published.tmp")], []
        )


# ---------------------------------------------------------------------------
# publish_all_channels: cross-channel two-phase all-or-nothing install
# ---------------------------------------------------------------------------
class TestPublishAllChannelsAtomicity(unittest.TestCase):
    """A mid-loop staging failure (e.g. a later channel's checksum mismatch)
    must publish NOTHING -- not even channels that staged fine before it."""

    def setUp(self):
        import btag_efficiency

        self.btag_efficiency = btag_efficiency
        self.tmpdir = os.path.join(
            os.environ.get("TMPDIR", "/tmp"), f"btag_publish_all_{os.getpid()}"
        )
        os.makedirs(self.tmpdir, exist_ok=True)
        self.config = {
            "output_base": self.tmpdir,
            "workdir_name": "wd",
            "era": "test",
            # Match the two synthetic channels below so the "required
            # channels present" gate passes and the failure under test comes
            # specifically from the checksum step.
            "required_channels": ["et", "mt"],
        }

    def tearDown(self):
        import shutil

        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def _make_channel_payload(self, channel):
        src_dir = os.path.join(self.tmpdir, f"staging_{channel}")
        os.makedirs(src_dir, exist_ok=True)
        payload = []
        for name, text in (
            ("btag_efficiency.json", f'{{"channel": "{channel}"}}'),
            ("btag_efficiency.json.gz", f"binary-ish-{channel}"),
        ):
            p = os.path.join(src_dir, name)
            with open(p, "w") as f:
                f.write(text)
            payload.append(p)
        checksums = {os.path.basename(p): _sha256(p) for p in payload}
        return payload, checksums

    def _provenance(self, channel, payload, checksums):
        return {
            "channel": channel,
            "status": "passed",
            "payload_files": payload,
            "checksums": checksums,
        }

    def _target(self, channel):
        return os.path.join(self.tmpdir, "wd", "test", "published", channel)

    def test_second_channel_checksum_failure_leaves_first_channel_unpublished(self):
        payload_et, checksums_et = self._make_channel_payload("et")
        payload_mt, checksums_mt = self._make_channel_payload("mt")
        # Simulate a corrupted checksum for the second channel's payload.
        checksums_mt["btag_efficiency.json"] = "0" * 64

        provenances = [
            self._provenance("et", payload_et, checksums_et),
            self._provenance("mt", payload_mt, checksums_mt),
        ]

        with self.assertRaises(RuntimeError):
            self.btag_efficiency.publish_all_channels(provenances, self.config)

        # Neither channel was published -- not even 'et', which staged fine
        # before 'mt' failed.
        self.assertFalse(os.path.exists(self._target("et")))
        self.assertFalse(os.path.exists(self._target("mt")))

        # No stray tmp dirs left behind for either channel.
        published_dir = os.path.join(self.tmpdir, "wd", "test", "published")
        leftover = []
        if os.path.isdir(published_dir):
            leftover = [d for d in os.listdir(published_dir) if ".tmp-" in d]
        self.assertEqual(leftover, [])

    def test_all_channels_publish_when_all_stage_cleanly(self):
        payload_et, checksums_et = self._make_channel_payload("et")
        payload_mt, checksums_mt = self._make_channel_payload("mt")
        provenances = [
            self._provenance("et", payload_et, checksums_et),
            self._provenance("mt", payload_mt, checksums_mt),
        ]

        self.btag_efficiency.publish_all_channels(provenances, self.config)

        self.assertTrue(
            os.path.isfile(os.path.join(self._target("et"), "btag_efficiency.json"))
        )
        self.assertTrue(
            os.path.isfile(os.path.join(self._target("mt"), "btag_efficiency.json"))
        )


# ---------------------------------------------------------------------------
# write_validation_report
# ---------------------------------------------------------------------------
class TestValidationReport(unittest.TestCase):
    def test_report_is_machine_readable_json_with_status_and_thresholds(self):
        tmpdir = os.environ.get("TMPDIR", "/tmp")
        path = os.path.join(tmpdir, f"btag_report_{os.getpid()}.json")
        thresholds = {"min_effective_denominator": 100.0, "max_stat_uncertainty": 0.25}
        report = write_validation_report(
            path,
            status="passed",
            thresholds=thresholds,
            merge_version="merge_v1",
            categories={"ttbar": {"b": {"merge_history": []}}},
            roundtrip={"ok": True, "n_checked": 12, "mismatches": []},
            era="test",
            channel="mt",
            channels_present=["mt"],
            required_channels=["mt"],
        )
        self.assertEqual(report["status"], "passed")
        self.assertEqual(report["thresholds"], thresholds)
        with open(path) as f:
            reloaded = json.load(f)
        self.assertEqual(reloaded["merge_version"], "merge_v1")
        self.assertEqual(reloaded["roundtrip"]["n_checked"], 12)
        self.assertIn("ttbar", reloaded["categories"])
        os.remove(path)


try:
    import ROOT  # noqa: F401

    HAS_ROOT = True
except Exception:  # pragma: no cover - environment dependent
    HAS_ROOT = False


def _write_probe_tree(path, pt, eta, flav, btag, weight):
    """Write a single-event ntuple with the four probe jet vectors + weight.

    Shared by the ROOT-gated tests below (probe-vector validation and the
    real pt-overflow diagnostic).
    """
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

    w = _array.array("f", [float(weight)])
    tree.Branch("weight", w, "weight/F")
    for vec, vals in (
        (pt_vec, pt),
        (eta_vec, eta),
        (flav_vec, flav),
        (btag_vec, btag),
    ):
        vec.clear()
        for v in vals:
            vec.push_back(v)
    tree.Fill()
    tree.Write()
    fout.Close()


@unittest.skipUnless(HAS_ROOT, "ROOT/RDataFrame not importable")
class TestProbeVectorFilter(unittest.TestCase):
    """The compiled equal-length assertion aborts the event loop on bad input."""

    def _write(self, path, pt, eta, flav, btag, weight):
        _write_probe_tree(path, pt, eta, flav, btag, weight)

    def test_unequal_length_vectors_raise(self):
        import ROOT

        import btag_efficiency

        ROOT.gROOT.SetBatch(True)
        tmpdir = os.environ.get("TMPDIR", "/tmp")
        path = os.path.join(tmpdir, f"btag_probe_bad_{os.getpid()}.root")
        # jet_eta_vec deliberately shorter than the other three vectors.
        self._write(path, [30.0, 40.0], [0.5], [5, 5], [0.9, 0.9], 1.0)

        chain = ROOT.TChain("ntuple")
        chain.Add(path)
        config = {
            "jet_pt_column": "jet_pt_vec",
            "jet_eta_column": "jet_eta_vec",
            "jet_flavor_column": "jet_hadronflavour_vec",
            "jet_btag_column": "jet_btag_value_vec",
            "weight_column": "weight",
            "jet_selection": "",
            "btag_working_points": {"M": 0.5},
            "jet_flavor_categories": {"b": 5},
            "jet_pt_bins": [20.0, 1000.0],
            "jet_eta_bins": [0.0, 2.5],
            "tree": "ntuple",
        }
        with self.assertRaises(Exception) as ctx:
            btag_efficiency.calculate_efficiency_histograms(
                chain, config, "bad_proc", [20.0, 1000.0], [0.0, 2.5]
            )
        self.assertIn("btag_validation", str(ctx.exception))


@unittest.skipUnless(HAS_ROOT, "ROOT/RDataFrame not importable")
class TestRawPtOverflowDiagnostic(unittest.TestCase):
    """``raw_flow_counts["pt_gt_1000"]`` must count real jets above the top pt
    edge (ROOT's Histo2D overflow bin), never be structurally zero."""

    def test_jet_above_top_edge_is_counted_via_root_overflow(self):
        import ROOT

        import btag_efficiency

        ROOT.gROOT.SetBatch(True)
        tmpdir = os.environ.get("TMPDIR", "/tmp")
        path = os.path.join(tmpdir, f"btag_overflow_{os.getpid()}.root")
        # One event with two jets: one in-range (pt=500) and one above the
        # top pt edge (pt=1500, edge at 1000) -- the latter must land in
        # ROOT's pt-overflow bin.
        _write_probe_tree(
            path,
            pt=[500.0, 1500.0],
            eta=[0.5, 0.5],
            flav=[5, 5],
            btag=[0.9, 0.9],
            weight=1.0,
        )

        chain = ROOT.TChain("ntuple")
        chain.Add(path)
        pt_bins = [20.0, 1000.0]
        eta_bins = [0.0, 2.5]
        config = {
            "jet_pt_column": "jet_pt_vec",
            "jet_eta_column": "jet_eta_vec",
            "jet_flavor_column": "jet_hadronflavour_vec",
            "jet_btag_column": "jet_btag_value_vec",
            "weight_column": "weight",
            "jet_selection": "",
            "btag_working_points": {"M": 0.5},
            "jet_flavor_categories": {"b": 5},
            "jet_pt_bins": pt_bins,
            "jet_eta_bins": eta_bins,
            "tree": "ntuple",
        }

        _histograms, accumulators = btag_efficiency.calculate_efficiency_histograms(
            chain, config, "overflow_proc", pt_bins, eta_bins
        )

        # The scalar is real and threaded onto the accumulator.
        acc = accumulators["M"]["b"]
        self.assertEqual(acc.raw_pt_overflow, 1.0)

        # ... and it reaches the report-facing raw_flow_counts() diagnostic.
        binned = BinnedAccumulators(
            accumulators={"M": acc}, pt_bins=pt_bins, eta_bins=eta_bins
        )
        counts = raw_flow_counts(binned)
        self.assertEqual(counts["pt_gt_1000"], 1.0)


if __name__ == "__main__":
    unittest.main()
