"""Composable YAML loading for TauFakeFactors configuration files."""

from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

from ruamel.yaml import YAML


def _read_yaml(path: Path, yaml: YAML) -> Dict[str, Any]:
    with path.open("r") as stream:
        document = yaml.load(stream)
    if document is None:
        return {}
    if not isinstance(document, Mapping):
        raise ValueError(f"Configuration '{path}' must contain a YAML mapping.")
    return dict(document)


def _load_specific(path: Path, yaml: YAML, stack: Iterable[Path]) -> Dict[str, Any]:
    resolved = path.resolve()
    ancestors = tuple(stack)
    if resolved in ancestors:
        chain = " -> ".join(str(item) for item in (*ancestors, resolved))
        raise ValueError(f"Cyclic base_config_file chain: {chain}")

    document = _read_yaml(resolved, yaml)
    merged: Dict[str, Any] = {}

    base_file = document.pop("base_config_file", None)
    if base_file:
        merged.update(
            _load_specific(
                resolved.parent / str(base_file), yaml, (*ancestors, resolved)
            )
        )

    processes_file = document.pop("processes_file", None)
    if processes_file:
        processes_document = _read_yaml(resolved.parent / str(processes_file), yaml)
        if "processes" not in processes_document:
            raise ValueError(
                f"Processes file '{processes_file}' referenced by '{resolved}' "
                "does not define a top-level 'processes' mapping."
            )
        merged["processes"] = dict(processes_document["processes"])

    local_processes = document.pop("processes", None)
    if local_processes is not None:
        combined_processes = dict(merged.get("processes", {}))
        combined_processes.update(local_processes)
        # A local ``name: null`` removes an inherited process/pool; the
        # measurement requires disjoint pools, so re-pooling a member needs
        # the ability to drop the pool it came from.
        document["processes"] = {
            name: entry for name, entry in combined_processes.items() if entry is not None
        }

    inherited_channels = "channels" in merged
    merged.update(document)
    if "channel" in document and "channels" not in document and inherited_channels:
        merged.pop("channels", None)
    return merged


def load_config(
    config_file: Any,
    *,
    defaults: Optional[Mapping[str, Any]] = None,
    yaml: Optional[YAML] = None,
) -> Dict[str, Any]:
    """Load common settings, optional inherited config, and one config file.

    ``base_config_file`` reuses all settings from another YAML file in the same
    directory. ``processes_file`` imports only its top-level ``processes`` map;
    a local ``processes`` map can then add or replace individual processes, and
    a local ``<name>: null`` removes an inherited one. Both references are
    resolved relative to the file containing them.
    """

    path = Path(config_file).resolve()
    yaml_loader = yaml or YAML(typ="rt")
    config: Dict[str, Any] = dict(defaults or {})

    common_path = path.parent / "common_settings.yaml"
    if common_path.exists():
        config.update(_read_yaml(common_path, yaml_loader))

    config.update(_load_specific(path, yaml_loader, ()))
    return config
