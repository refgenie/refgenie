"""The lazy asset-path view, ``Refgenie.paths()`` (refgenie/core/paths.py).

Component tier: it reads the session fixture's built rCRSd fasta asset.
"""

import jinja2
import pytest

pytestmark = pytest.mark.component


def _seek(r, registry_path):
    """The path a registry path seeks to: the view is keyed by alias, as it is."""
    return r.asset.seek_components(r.parse_asset_registry_path(registry_path))


def test_leaf_matches_seek_and_is_str(refgenie_session):
    leaf = refgenie_session.paths()["rCRSd"]["fasta"]["fasta"]
    assert isinstance(leaf, str)
    assert leaf == _seek(refgenie_session, "rCRSd/fasta.fasta")


def test_missing_levels_raise_key_error(refgenie_session):
    paths = refgenie_session.paths()
    with pytest.raises(KeyError):
        paths["nope"]
    with pytest.raises(KeyError):
        paths["rCRSd"]["nope"]
    with pytest.raises(KeyError):
        paths["rCRSd"]["fasta"]["nope"]
    assert "nope" not in paths
    assert "rCRSd" in paths
    assert "fasta" in paths["rCRSd"]


def test_jinja_defined_guards(refgenie_session):
    template = jinja2.Template(
        "{% if refgenie['rCRSd'].fasta is defined %}{{ refgenie['rCRSd'].fasta.fasta }}"
        "{% endif %}|{% if refgenie['rCRSd'].nope is defined %}x{% endif %}"
    )
    rendered = template.render(refgenie=refgenie_session.paths())
    expected = _seek(refgenie_session, "rCRSd/fasta.fasta")
    assert rendered == f"{expected}|"


def test_each_leaf_is_looked_up_once(refgenie_session, monkeypatch):
    calls = []
    real_seek = refgenie_session.asset.seek_components

    def spy(*args, **kwargs):
        calls.append((args, kwargs))
        return real_seek(*args, **kwargs)

    monkeypatch.setattr(refgenie_session.asset, "seek_components", spy)
    paths = refgenie_session.paths()
    first = paths["rCRSd"]["fasta"]["fasta"]
    second = paths["rCRSd"]["fasta"]["fasta"]
    assert first == second
    assert len(calls) == 1


def test_to_dict_matches_lookups(refgenie_session):
    paths = refgenie_session.paths()
    walked = paths.to_dict()
    assert set(walked) == {"rCRSd"}
    assert set(walked["rCRSd"]) == {"fasta"}
    assert "fasta" in walked["rCRSd"]["fasta"]
    for seek_key, value in walked["rCRSd"]["fasta"].items():
        assert isinstance(value, str)
        assert value == paths["rCRSd"]["fasta"][seek_key]
