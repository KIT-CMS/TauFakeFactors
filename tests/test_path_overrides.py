"""Tests for CLI/env/config path-override precedence (Task 16).

Covers two things:

1. ``helper.functions.resolve_path_setting`` — the pure precedence function
   (CLI > env > config > default) used by ``preselection.py`` (for
   ``ntuple_path``/``output_path``) and ``btag_efficiency.py`` (for
   ``file_path``/``output_base``), plus the new CLI flags that feed it.
2. The "no literal user paths" rule for versioned 2018 b-tag efficiency
   configs under ``configs/btag_efficiency/2018/`` — Run-3 era configs
   predate this rule and still ship literal ``/ceph/<user>`` and
   ``/store/user/<user>`` paths in their ``common_settings.yaml``, so they
   are explicitly exempted via ``LEGACY_USER_PATH_ERAS``.
"""

import re
from pathlib import Path
from typing import List, Tuple

import pytest

import helper.functions as func

REPO_ROOT = Path(__file__).resolve().parents[1]
BTAG_EFFICIENCY_DIR = REPO_ROOT / "configs" / "btag_efficiency"

# Run-3 eras predate the CLI/env path-override work (Task 16) and still ship
# literal /ceph/<user> and /store/user/<user> paths in their
# common_settings.yaml (verified: e.g. configs/btag_efficiency/2024/
# common_settings.yaml has output_path: "/ceph/sgiappic/CMS_FF/btag/"). They
# are explicitly exempted from the no-user-paths rule below. Every era added
# from 2018 onward (starting with Task 19) must NOT contain literal user
# paths -- CLI/env overrides (or portable relative defaults) must be used
# instead.
LEGACY_USER_PATH_ERAS = (
    "2022preEE",
    "2022postEE",
    "2023preBPix",
    "2023postBPix",
    "2024",
    "2025",
)

# Matches literal "/ceph/<user>" and "/store/user/<user>" path prefixes.
USER_PATH_RE = re.compile(r"(/ceph|/store/user)/[A-Za-z0-9]+")


def _iter_yaml_files(directory: Path) -> List[Path]:
    return sorted(directory.rglob("*.yaml"))


def _user_path_violations(yaml_paths: List[Path]) -> List[Tuple[Path, str]]:
    """Return (path, matched_string) for every literal user path found."""
    violations: List[Tuple[Path, str]] = []
    for path in yaml_paths:
        text = path.read_text()
        for match in USER_PATH_RE.finditer(text):
            violations.append((path, match.group(0)))
    return violations


# ---------------------------------------------------------------------------
# resolve_path_setting: precedence unit tests
# ---------------------------------------------------------------------------
ENV_VAR = "TFF_TEST_PATH_OVERRIDE"


@pytest.fixture(autouse=True)
def _clean_test_env(monkeypatch):
    """Ensure ENV_VAR never leaks in from the real environment/other tests."""
    monkeypatch.delenv(ENV_VAR, raising=False)


def test_cli_only():
    result = func.resolve_path_setting({}, "some_path", "/cli/value", ENV_VAR)
    assert result == "/cli/value"


def test_env_only(monkeypatch):
    monkeypatch.setenv(ENV_VAR, "/env/value")
    result = func.resolve_path_setting({}, "some_path", None, ENV_VAR)
    assert result == "/env/value"


def test_config_only():
    config = {"some_path": "/config/value"}
    result = func.resolve_path_setting(config, "some_path", None, ENV_VAR)
    assert result == "/config/value"


def test_default_only():
    result = func.resolve_path_setting({}, "some_path", None, ENV_VAR, default="/default/value")
    assert result == "/default/value"


def test_none_set_raises_missing_path_setting_error():
    with pytest.raises(func.MissingPathSettingError):
        func.resolve_path_setting({}, "some_path", None, ENV_VAR)


def test_missing_path_setting_error_is_a_key_error():
    # "KeyError-derived" per spec: catchable as a plain KeyError too.
    assert issubclass(func.MissingPathSettingError, KeyError)
    with pytest.raises(KeyError):
        func.resolve_path_setting({}, "some_path", None, ENV_VAR)


def test_falsy_cli_value_falls_through_to_env(monkeypatch):
    """An empty-string / falsy CLI value must not shadow the env var."""
    monkeypatch.setenv(ENV_VAR, "/env/value")
    result = func.resolve_path_setting({}, "some_path", "", ENV_VAR)
    assert result == "/env/value"


# --- pairwise precedence -----------------------------------------------------
def test_cli_beats_env(monkeypatch):
    monkeypatch.setenv(ENV_VAR, "/env/value")
    result = func.resolve_path_setting({}, "some_path", "/cli/value", ENV_VAR)
    assert result == "/cli/value"


def test_cli_beats_config():
    config = {"some_path": "/config/value"}
    result = func.resolve_path_setting(config, "some_path", "/cli/value", ENV_VAR)
    assert result == "/cli/value"


def test_env_beats_config(monkeypatch):
    monkeypatch.setenv(ENV_VAR, "/env/value")
    config = {"some_path": "/config/value"}
    result = func.resolve_path_setting(config, "some_path", None, ENV_VAR)
    assert result == "/env/value"


def test_cli_beats_env_and_config(monkeypatch):
    monkeypatch.setenv(ENV_VAR, "/env/value")
    config = {"some_path": "/config/value"}
    result = func.resolve_path_setting(config, "some_path", "/cli/value", ENV_VAR)
    assert result == "/cli/value"


def test_config_beats_default():
    config = {"some_path": "/config/value"}
    result = func.resolve_path_setting(
        config, "some_path", None, ENV_VAR, default="/default/value"
    )
    assert result == "/config/value"


def test_env_beats_default(monkeypatch):
    monkeypatch.setenv(ENV_VAR, "/env/value")
    result = func.resolve_path_setting({}, "some_path", None, ENV_VAR, default="/default/value")
    assert result == "/env/value"


def test_pure_function_does_not_mutate_config_or_environment(monkeypatch):
    """resolve_path_setting must have no side effects: config dict and the
    environment are left exactly as they were.
    """
    import os

    config = {"some_path": "/config/value"}
    config_copy = dict(config)
    func.resolve_path_setting(config, "some_path", "/cli/value", ENV_VAR)
    assert config == config_copy
    assert ENV_VAR not in os.environ


# ---------------------------------------------------------------------------
# CLI args exist: preselection.py / btag_efficiency.py
# ---------------------------------------------------------------------------
def test_preselection_cli_has_ntuple_path_and_output_path_flags():
    import preselection

    help_text = preselection.build_arg_parser().format_help()
    assert "--ntuple-path" in help_text
    assert "--output-path" in help_text


def test_btag_efficiency_cli_has_file_path_and_output_path_flags():
    import btag_efficiency

    help_text = btag_efficiency.build_arg_parser().format_help()
    assert "--file-path" in help_text
    assert "--output-path" in help_text


def test_preselection_parses_new_flags():
    import preselection

    args = preselection.build_arg_parser().parse_args(
        ["--config-file", "cfg.yaml", "--ntuple-path", "/n", "--output-path", "/o"]
    )
    assert args.ntuple_path == "/n"
    assert args.output_path == "/o"


def test_preselection_new_flags_default_to_none():
    import preselection

    args = preselection.build_arg_parser().parse_args(["--config-file", "cfg.yaml"])
    assert args.ntuple_path is None
    assert args.output_path is None


def test_btag_efficiency_parses_new_flags():
    import btag_efficiency

    args = btag_efficiency.build_arg_parser().parse_args(
        ["--config-file", "cfg.yaml", "--file-path", "/f", "--output-path", "/o"]
    )
    assert args.file_path == "/f"
    assert args.output_path == "/o"


def test_btag_efficiency_new_flags_default_to_none():
    import btag_efficiency

    args = btag_efficiency.build_arg_parser().parse_args(["--config-file", "cfg.yaml"])
    assert args.file_path is None
    assert args.output_path is None


# ---------------------------------------------------------------------------
# No-user-paths rule for versioned 2018 configs
# ---------------------------------------------------------------------------
def iter_era_dirs():
    if not BTAG_EFFICIENCY_DIR.is_dir():
        return
    for era_dir in sorted(BTAG_EFFICIENCY_DIR.iterdir()):
        if era_dir.is_dir():
            yield era_dir.name, era_dir


ERA_DIRS = list(iter_era_dirs())
ERA_DIR_IDS = [era for era, _ in ERA_DIRS]


@pytest.mark.parametrize("era,era_dir", ERA_DIRS, ids=ERA_DIR_IDS)
def test_no_literal_user_paths_outside_legacy_eras(era: str, era_dir: Path) -> None:
    """Every era directory must contain no literal /ceph or /store/user/<user>
    path, UNLESS it is listed in LEGACY_USER_PATH_ERAS (pre-Task-16 Run-3
    configs).
    """
    if era in LEGACY_USER_PATH_ERAS:
        pytest.skip(
            f"{era} is a legacy Run-3 era, exempt from the no-user-paths rule (Task 16)"
        )
    violations = _user_path_violations(_iter_yaml_files(era_dir))
    assert not violations, f"{era_dir}: literal user paths found: {violations}"


def test_run3_eras_are_the_ones_actually_exempted() -> None:
    """Sanity check the exemption is not vacuous: at least one legacy Run-3
    era currently DOES contain a literal user path (if this ever stops being
    true, LEGACY_USER_PATH_ERAS should be revisited/shrunk).
    """
    legacy_dirs = [
        era_dir for era, era_dir in ERA_DIRS if era in LEGACY_USER_PATH_ERAS
    ]
    assert legacy_dirs, "expected at least one legacy Run-3 era config dir to exist"

    found_user_path = any(
        _user_path_violations(_iter_yaml_files(era_dir)) for era_dir in legacy_dirs
    )
    assert found_user_path, (
        "expected at least one legacy Run-3 config to contain a literal user "
        "path; if this fails, the LEGACY_USER_PATH_ERAS exemption may be obsolete"
    )


def test_2018_configs_have_no_user_paths_or_do_not_exist_yet() -> None:
    """configs/btag_efficiency/2018 (added by Task 19) must contain no
    literal /ceph/<user> or /store/user/<user> path.

    The directory does not exist yet at the time Task 16 runs, so this
    assertion passes vacuously today -- that is asserted explicitly (rather
    than silently skipped) via the check that "2018" is absent from
    LEGACY_USER_PATH_ERAS, so the rule is active and enforced the moment
    Task 19 adds the directory.
    """
    assert "2018" not in LEGACY_USER_PATH_ERAS

    era_dir = BTAG_EFFICIENCY_DIR / "2018"
    if not era_dir.is_dir():
        return  # vacuous pass, explicitly noted above

    violations = _user_path_violations(_iter_yaml_files(era_dir))
    assert not violations, f"literal user paths found in 2018 configs: {violations}"


def test_rule_fails_for_hypothetical_2018_config_with_user_path(tmp_path: Path) -> None:
    """Prove the no-user-paths rule has teeth (does not pass everything): a
    synthetic 2018-like config containing a literal /ceph/<user> path (and a
    /store/user/<user> ntuple path) must be flagged as a violation.
    """
    fake_era_dir = tmp_path / "configs" / "btag_efficiency" / "2018"
    fake_era_dir.mkdir(parents=True)
    bad_config = fake_era_dir / "common_settings.yaml"
    bad_config.write_text(
        'era: "2018"\n'
        'ntuple_path: "root://cmsdcache-kit-disk.gridka.de//store/user/sgiappic/CROWN/ntuples/foo"\n'
        'output_path: "/ceph/sgiappic/CMS_FF/btag/"\n'
        'file_path: "/ceph/sgiappic/CMS_FF/btag/"\n'
    )

    violations = _user_path_violations(_iter_yaml_files(fake_era_dir))
    assert violations, "expected the no-user-paths rule to flag /ceph/sgiappic"

    matched_strings = {matched for _, matched in violations}
    assert any(m.startswith("/ceph") for m in matched_strings)
    assert any(m.startswith("/store/user") for m in matched_strings)


def test_rule_passes_for_hypothetical_2018_config_without_user_paths(tmp_path: Path) -> None:
    """Counterpart to the previous test: a synthetic 2018-like config using
    only relative/portable paths (as Task 19 is expected to ship, resolved
    via the Task-16 CLI/env overrides) must NOT be flagged.
    """
    clean_era_dir = tmp_path / "configs" / "btag_efficiency" / "2018"
    clean_era_dir.mkdir(parents=True)
    clean_config = clean_era_dir / "common_settings.yaml"
    clean_config.write_text(
        'era: "2018"\n'
        'ntuple_path: "ntuples"\n'
        'output_path: "output"\n'
        'file_path: "output"\n'
    )

    violations = _user_path_violations(_iter_yaml_files(clean_era_dir))
    assert not violations, f"unexpected violations for a portable-paths config: {violations}"
