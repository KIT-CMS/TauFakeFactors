"""Validate a b-tag efficiency map against a strict consumer's requirements.

Some consumers of the efficiency payload refuse an efficiency they cannot turn
into a per-jet weight, and abort the job rather than substitute a value.  The
bbtautau CROWN 2018 UParT consumer
(``cpp_addons/src/btag_sf_strict.cxx::multi_wp_event_weight``) is one: it
throws unless, for every bin it evaluates,

    isfinite(eff) and 0 < eff <= 1                  (value range)
    eff(looser) - eff(tighter) > 0, strictly        (no ties)
    1 - eff(loosest) > 0                            (untagged probability)

``btag_efficiency.py`` can produce all three violations without failing: an
empty bin is written as ``eff = 0.0`` (with a warning), a bin where no jet
passes the working point is written as ``eff = 0.0`` (silently), and a bin
where every jet passes is written as ``eff = 1.0`` (also silently).  A payload
built from thin per-process samples therefore looks successful and aborts the
production instead.

The checks here are deliberately self-contained: they need only the efficiency
grids the generator just produced, no scale-factor payload.  The consumer has
a fourth requirement -- that the resulting jet weight ``num / denom`` stays
finite and positive -- which cannot be decided without the experiment's scale
factors, so it lives in the standalone
``/work/sdaigler/bbtautau/check_btag_payload.py --sf-check`` instead.  Passing
the checks here is necessary, not sufficient.
"""
from typing import Dict, List, Tuple


def working_points_loose_to_tight(config: Dict) -> List[str]:
    """Order the configured working points by increasing discriminator cut.

    Sorted by threshold rather than by the order they happen to appear in the
    YAML, so a reordered config cannot silently invert the monotonicity test.
    """
    return [
        name for name, _ in sorted(
            config["btag_working_points"].items(), key=lambda item: item[1]
        )
    ]


def _edge_label(edges: List[float], index: int, name: str) -> str:
    """Name one bin in physical units, e.g. ``pt=[140,200)``."""
    try:
        return f"{name}=[{edges[index]:g},{edges[index + 1]:g})"
    except (IndexError, TypeError):
        return f"{name}_bin={index}"


def _cell_label(
    pt_bins: List[float],
    eta_bins: List[float],
    i_pt: int,
    i_eta: int,
) -> str:
    return (
        _edge_label(eta_bins, i_eta, "eta")
        + "/"
        + _edge_label(pt_bins, i_pt, "pt")
    )


def validate_efficiencies(
    all_efficiencies: Dict[str, Dict[str, Dict[str, List[List[float]]]]],
    all_bins: Dict[str, Dict[str, Tuple[List[float], List[float]]]],
    config: Dict,
) -> List[str]:
    """Return one human-readable message per violation; empty means valid.

    Bins are named in GeV and |eta| rather than by index, because the remedy is
    a coarser ``jet_pt_bins`` / ``jet_eta_bins`` for that process and flavour,
    and an index does not say which one.
    """
    problems: List[str] = []

    required = config.get("required_sample_types") or []
    for name in sorted(set(required) - set(all_efficiencies)):
        problems.append(
            f"MISSING SAMPLE TYPE  '{name}' is required but no measurement "
            f"produced it (a skipped process leaves the payload short, and the "
            f"consumer fails at run time on the first event that needs it)"
        )

    wps = working_points_loose_to_tight(config)
    flavors = config["jet_flavor_categories"]

    for sample_type, per_wp in sorted(all_efficiencies.items()):
        flavor_bins = all_bins[sample_type]
        for flavor_name in flavors:
            pt_bins, eta_bins = flavor_bins[flavor_name]
            n_pt, n_eta = len(pt_bins) - 1, len(eta_bins) - 1

            for i_pt in range(n_pt):
                for i_eta in range(n_eta):
                    cell = _cell_label(pt_bins, eta_bins, i_pt, i_eta)
                    where = f"{sample_type}/flavor={flavors[flavor_name]}/{cell}"

                    ladder = []
                    for wp in wps:
                        value = per_wp[wp][flavor_name][i_pt][i_eta]
                        ladder.append(value)
                        if not (value == value) or value in (
                            float("inf"), float("-inf")
                        ):
                            problems.append(
                                f"NOT FINITE           {where}/{wp}: {value}"
                            )
                        elif value <= 0.0:
                            problems.append(
                                f"eff <= 0             {where}/{wp}: {value}"
                            )
                        elif value > 1.0:
                            problems.append(
                                f"eff > 1              {where}/{wp}: {value}"
                            )

                    # Strictly decreasing from the loosest to the tightest
                    # working point: the consumer's denominator for a jet that
                    # passes `low` but not `high` is eff(low) - eff(high), so a
                    # tie divides by zero.
                    for looser, tighter, low, high in zip(
                        wps[:-1], wps[1:], ladder[:-1], ladder[1:]
                    ):
                        if low - high <= 0.0:
                            problems.append(
                                f"NOT MONOTONIC        {where}: "
                                f"eff({looser})={low} must exceed "
                                f"eff({tighter})={high}"
                            )

                    # A jet that fails even the loosest working point is
                    # weighted with 1 - eff(loosest) in the denominator.
                    if 1.0 - ladder[0] <= 0.0:
                        problems.append(
                            f"DEGENERATE           {where}: "
                            f"1 - eff({wps[0]})={1.0 - ladder[0]} "
                            f"must be positive"
                        )

    return problems


def resolve_mode(config: Dict, cli_override: str = None) -> str:
    """Return the validation mode: ``off``, ``warn`` or ``error``.

    ``--validate-payload`` on the command line wins over the config key
    ``validate_payload``; the default is ``off`` so that eras whose consumer
    tolerates zero bins keep working unchanged.
    """
    mode = cli_override or config.get("validate_payload", "off")
    mode = str(mode).lower()
    if mode not in ("off", "warn", "error"):
        raise ValueError(
            f"validate_payload must be 'off', 'warn' or 'error', got '{mode}'."
        )
    return mode
