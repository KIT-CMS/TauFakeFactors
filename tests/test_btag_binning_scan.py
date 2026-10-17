"""Unit tests of the grid logic of scripts/btag_binning_scan.py (ROOT only for the empty-process case)."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "scripts"))

import btag_binning_scan as scan  # noqa: E402

PT = [20.0, 30.0, 50.0, 100.0, 1000.0]
ETA = [0.0, 1.4, 2.4]


def counts(per_bin):
    """[region, pt, eta] with the same count in every region and |eta| bin, per pT bin."""
    grid = np.zeros((3, len(PT) - 1, len(ETA) - 1))
    for i, value in enumerate(per_bin):
        grid[:, i, :] = value
    return grid


def test_rebin_sums_fine_cells():
    grid = counts([1, 2, 3, 4])
    coarse = scan.rebin(grid, PT, ETA, [20.0, 50.0, 1000.0], [0.0, 2.4])
    assert coarse.shape == (3, 2, 1)
    assert coarse[0, 0, 0] == (1 + 2) * 2 and coarse[0, 1, 0] == (3 + 4) * 2


def test_rebin_rejects_edges_off_the_fine_grid():
    with pytest.raises(ValueError):
        scan.rebin(counts([1, 1, 1, 1]), PT, ETA, [20.0, 40.0, 1000.0], ETA)


def test_full_grid_kept_when_every_cell_passes():
    assert scan.propose({"mt": counts([9, 9, 9, 9])}, PT, ETA, PT, [ETA, [0.0, 2.4]], 5) == (PT, ETA)


def test_highest_failing_bin_merges_into_its_lower_neighbour():
    # [100,1000) fails: merged into [50,100) -> [50,1000) holds 6 per cell
    pt, eta = scan.propose({"mt": counts([9, 9, 3, 3])}, PT, ETA, PT, [ETA, [0.0, 2.4]], 5)
    assert pt == [20.0, 30.0, 50.0, 1000.0] and eta == ETA


def test_lowest_bin_merges_upwards_and_channels_are_combined_by_minimum():
    pt, _ = scan.propose({"mt": counts([3, 9, 9, 9]), "tt": counts([9, 9, 9, 9])}, PT, ETA, PT, [ETA], 5)
    assert pt == [20.0, 50.0, 100.0, 1000.0]


def test_single_eta_bin_wins_when_it_keeps_more_cells():
    # forward cells hold 1 jet per region: two |eta| bins fail even as one pT bin (4 < 5), one |eta| bin keeps all
    grid = counts([4, 4, 4, 4])
    grid[:, :, 1] = 1
    assert scan.propose({"mt": grid}, PT, ETA, PT, [ETA, [0.0, 2.4]], 5) == (PT, [0.0, 2.4])


def test_tie_in_cells_goes_to_two_eta_bins():
    # two |eta| bins: [20,50,1000] x 2 = 4 cells; one |eta| bin: all four pT bins x 1 = 4 cells
    pt, eta = scan.propose({"mt": counts([3, 3, 3, 3])}, PT, ETA, PT, [ETA, [0.0, 2.4]], 5)
    assert eta == ETA and pt == [20.0, 50.0, 1000.0]


def test_no_proposal_when_even_one_cell_is_too_thin():
    assert scan.propose({"mt": counts([0.5, 0.5, 0.5, 0.5])}, PT, ETA, PT, [ETA, [0.0, 2.4]], 5) is None


def test_region_names_follow_the_working_points():
    wps = scan.wps_loose_to_tight({"btag_working_points": {"T": 0.54, "M": 0.16}})
    assert [name for name, _ in wps] == ["M", "T"]
    assert scan.region_names(wps) == ["fail M", "M not T", "pass T"]


def test_a_process_without_events_counts_as_empty(tmp_path):
    # preselection.py writes a process without selected events as a tree without branches
    ROOT = pytest.importorskip("ROOT")
    path = str(tmp_path / "vbf_hbb.root")
    handle = ROOT.TFile(path, "RECREATE")
    ROOT.TTree("ntuple", "ntuple").Write()
    handle.Close()
    config = {
        "btag_working_points": {"T": 0.54, "M": 0.16},
        "jet_flavor_categories": {"b": 5, "c": 4, "light": 0},
        "jet_pt_column": "btag_probe_jet_pt", "jet_eta_column": "btag_probe_jet_eta",
        "jet_flavor_column": "btag_probe_jet_hadron_flavour", "jet_btag_column": "btag_probe_jet_upart",
    }
    assert scan.filled_files([path], "ntuple") == []
    result = scan.count_jets([path], config, PT, ETA, "ntuple")
    assert set(result) == {"b", "c", "light"}
    assert all(grid.shape == (3, len(PT) - 1, len(ETA) - 1) and not grid.any() for grid in result.values())
