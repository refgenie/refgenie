"""The looper pre-submit hook, ``refgenie.integrations.looper.populate``."""

import json
import subprocess
import sys

import pytest

from refgenie.integrations import looper


def _namespaces(**var_templates):
    return {
        "pipeline": {"pipeline_name": "demo", "var_templates": dict(var_templates)},
        "sample": {},
        "project": {},
    }


def _looper_merge(namespaces, returned):
    """Copy of looper's leaf-by-leaf merge (looper/conductor.py, `_exec_pre_submit`)."""
    for namespace, mapping in returned.items():
        for key, val in mapping.items():
            namespaces[namespace][key] = val


@pytest.fixture
def bound(refgenie_session, monkeypatch):
    """Make the hook use the session fixture instead of opening a database."""
    monkeypatch.setattr(looper, "_open", lambda _: (refgenie_session, refgenie_session.paths()))
    return refgenie_session


def test_refgenie_namespace_matches_seek(bound):
    ns = _namespaces()
    looper.populate(ns)
    expected = bound.asset.seek_components(bound.parse_asset_registry_path("rCRSd/fasta.fasta"))
    assert ns["refgenie"]["rCRSd"]["fasta"]["fasta"] == expected


def test_refgenie_urls_in_the_pipeline_block_are_resolved(bound):
    ns = _namespaces(x="refgenie://rCRSd/fasta", n=3, missing="refgenie://rCRSd/nope")
    _looper_merge(ns, looper.populate(ns))
    var_templates = ns["pipeline"]["var_templates"]
    assert var_templates["x"] == bound.asset.seek_components(
        bound.parse_asset_registry_path("rCRSd/fasta")
    )
    # Scalars pass through; a missing asset only warns and stays as written.
    assert var_templates["n"] == 3
    assert var_templates["missing"] == "refgenie://rCRSd/nope"
    assert ns["pipeline"]["pipeline_name"] == "demo"


def test_refgenie_is_built_once_per_config(refgenie_session, monkeypatch):
    built = []

    def fake_refgenie(**kwargs):
        built.append(kwargs)
        return refgenie_session

    looper._open.cache_clear()
    monkeypatch.setattr(looper, "Refgenie", fake_refgenie)
    try:
        looper.populate(_namespaces(refgenie_db_config="/some/config.yaml"))
        looper.populate(_namespaces(refgenie_db_config="/some/config.yaml"))
    finally:
        looper._open.cache_clear()
    assert len(built) == 1
    assert str(built[0]["database_config_path"]) == "/some/config.yaml"


# ---------------------------------------------------------------------------
# refgenie_db_config that names no file: PEPATAC sets it to
# "${REFGENIE_DB_CONFIG_PATH}", which looper passes as written when unset.
# ---------------------------------------------------------------------------


@pytest.fixture
def db_config_files(refgenie_session, tmp_path, monkeypatch):
    """Run the hook in an empty working folder, with the real config-path handling.

    Returns the working folder and the configured (fallback) config file.
    """
    from refgenie.config import config
    from refgenie.db.engine import load_database_config

    work = tmp_path / "work"
    work.mkdir()
    monkeypatch.chdir(work)
    monkeypatch.delenv("REFGENIE_DB_CONFIG_PATH", raising=False)
    fallback = tmp_path / "home" / "refgenie_db_config.yaml"
    monkeypatch.setattr(config, "database_config_path", fallback)

    def fake_refgenie(database_config_path=None):
        load_database_config(database_config_path)  # creates the file it reads
        return refgenie_session

    monkeypatch.setattr(looper, "Refgenie", fake_refgenie)
    looper._open.cache_clear()
    yield work, fallback
    looper._open.cache_clear()


@pytest.mark.parametrize("value", ["${REFGENIE_DB_CONFIG_PATH}", "$REFGENIE_DB_CONFIG_PATH"])
def test_unexpanded_db_config_falls_back(db_config_files, value):
    work, fallback = db_config_files
    looper.populate(_namespaces(refgenie_db_config=value))
    assert list(work.iterdir()) == []
    assert fallback.is_file()


@pytest.mark.parametrize("value", ["", "  ", None])
def test_empty_db_config_falls_back(db_config_files, value):
    work, fallback = db_config_files
    looper.populate(_namespaces(refgenie_db_config=value))
    assert list(work.iterdir()) == []
    assert fallback.is_file()


def test_db_config_variable_wins_when_set(db_config_files, tmp_path, monkeypatch):
    work, fallback = db_config_files
    from_env = tmp_path / "env" / "refgenie_db_config.yaml"
    monkeypatch.setenv("REFGENIE_DB_CONFIG_PATH", str(from_env))
    looper.populate(_namespaces(refgenie_db_config="${REFGENIE_DB_CONFIG_PATH}"))
    assert list(work.iterdir()) == []
    assert from_env.is_file()
    assert not fallback.exists()


# ---------------------------------------------------------------------------
# Importing the hook must leave the host's logging alone: looper's own output
# went quiet when it reset the root logger.
# ---------------------------------------------------------------------------

_ROOT_LOGGER_PROBE = """
import json, logging, sys

root = logging.getLogger()
logging.basicConfig(level=logging.INFO)  # the host program's own setup


def state():
    return {
        "handlers": [[type(h).__name__, id(h), h.level] for h in root.handlers],
        "level": root.level,
    }


before = state()
import refgenie
import refgenie.integrations.looper
from refgenie import Refgenie
print(json.dumps({"before": before, "after": state()}))
"""


@pytest.mark.component  # a fresh interpreter, so nothing has imported refgenie yet
def test_importing_the_hook_leaves_the_root_logger_alone():
    result = subprocess.run(
        [sys.executable, "-c", _ROOT_LOGGER_PROBE],
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
    states = json.loads(result.stdout)
    assert states["after"] == states["before"]
