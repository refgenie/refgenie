"""The looper pre-submit hook, ``refgenie.integrations.looper.populate``."""

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
