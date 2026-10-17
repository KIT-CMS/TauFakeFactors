"""Pool definitions of the b-tag efficiency configs: who is measured together, under which payload keys.

Standard library only (no ROOT, matplotlib or the fake-factor helpers with their wurlitzer import), so the config
tests run validate_process_config in every environment; btag_efficiency.py imports these from here.
"""
from typing import Dict, List


def get_input_processes(proc_conf: Dict, process: str) -> List[str]:
    """Return the input processes whose jets are pooled into one measurement.

    ``input_processes`` lets several MC processes share a single efficiency
    measurement.  Thin processes (``vbf_hbb`` has single-digit probe-jet counts
    in 2018) cannot fill a (pT, |eta|) x working-point grid on their own, and
    the resulting empty bins are written as ``eff = 0.0``
    (``btag_efficiency.histograms_to_efficiencies``), which the strict 2018 consumer
    rejects outright.  Without the key the process stands alone, as before.
    """
    inputs = proc_conf.get("input_processes", [process])
    if isinstance(inputs, str):
        inputs = [inputs]
    if not inputs:
        raise ValueError(
            f"Process '{process}': 'input_processes' must not be empty."
        )
    return list(inputs)


def get_output_sample_types(proc_conf: Dict, process: str) -> List[str]:
    """Return the correctionlib ``sample_type`` keys this measurement fills.

    A pooled measurement is written out once per member so that the payload
    keeps one entry per MC sample type -- the 2018 consumer looks the
    efficiency up under the sample's own name (``bjet_eff_sample_type`` is the
    identity for the SM profile), so every sample type it can be run on must be
    present as a key.  ``sample_types`` (list) and ``sample_type`` (str) are
    both accepted; the singular form is the pre-existing spelling.
    """
    if "sample_types" in proc_conf:
        keys = proc_conf["sample_types"]
        if isinstance(keys, str):
            keys = [keys]
        keys = list(keys)
    elif "sample_type" in proc_conf:
        keys = [proc_conf["sample_type"]]
    else:
        raise ValueError(
            f"Process '{process}': needs either 'sample_type' (str) or "
            f"'sample_types' (list) naming the payload key(s) it fills."
        )
    if not keys:
        raise ValueError(
            f"Process '{process}': 'sample_types' must not be empty."
        )
    return keys


def validate_process_config(config: Dict) -> None:
    """Fail loudly on process definitions that would silently corrupt a payload.

    ``run_for_channel`` accumulates results into ``all_efficiencies[key]``.  Two
    process entries claiming the same payload key therefore overwrite each
    other -- and because the multiprocessing path collects futures in
    completion order, *which* of them survives is not even deterministic.  The
    same holds for an input process pooled into two measurements: its jets
    would be counted twice, into two different grids.  Both are configuration
    errors with no legitimate use, so they abort here rather than produce a
    plausible-looking payload.
    """
    processes: Dict[str, Dict] = config["processes"]

    key_owner: Dict[str, str] = {}
    input_owner: Dict[str, str] = {}
    for process, proc_conf in processes.items():
        proc_conf = proc_conf or {}
        for key in get_output_sample_types(proc_conf, process):
            if key in key_owner:
                raise ValueError(
                    f"Payload key '{key}' is claimed by both process "
                    f"'{key_owner[key]}' and process '{process}'. Each "
                    f"correctionlib sample_type may be filled exactly once."
                )
            key_owner[key] = process
        for member in get_input_processes(proc_conf, process):
            if member in input_owner:
                raise ValueError(
                    f"Input process '{member}' is pooled into both "
                    f"'{input_owner[member]}' and '{process}'. Pools must be "
                    f"disjoint, otherwise its jets are counted twice."
                )
            input_owner[member] = process

    required = config.get("required_sample_types")
    if required:
        missing = sorted(set(required) - set(key_owner))
        extra = sorted(set(key_owner) - set(required))
        if missing or extra:
            raise ValueError(
                "Configured processes do not cover 'required_sample_types' "
                f"exactly: missing {missing or 'none'}, unexpected "
                f"{extra or 'none'}."
            )
