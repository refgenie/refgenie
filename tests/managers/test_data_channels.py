"""
Tests for refgenie.managers.sources data channels.

Covers data-channel add/remove/get round-trips, the documented ValueError
guards, and index-file fetching/parsing.

The index-parsing group reproduces a bug where a channel URL pointing at an
HTML page (instead of index.yaml) caused an opaque error:
    IndexFile() argument after ** must be a mapping, not str
"""

from unittest.mock import patch

import pytest

from refgenie.db.tables import DataChannel, DataChannelType


def _add(refgenie, name, type=DataChannelType.https, **kwargs):
    kwargs.setdefault("index_address", f"https://example.com/{name}/index.yaml")
    return refgenie.sources.add_channel(name=name, type=type, **kwargs)


class TestDataChannelCRUD:
    """add_channel / remove_channel / get_channel / list_channels round-trips."""

    def test_add_channel_persists_and_returns_it(self, refgenie_minimal):
        """add_channel returns the created channel and it appears in list_channels."""
        before = len(refgenie_minimal.sources.list_channels())

        channel = _add(
            refgenie_minimal,
            "test_channel",
            index_address="https://example.com/test_channel/index.yaml",
            description="a test channel",
        )
        assert isinstance(channel, DataChannel)
        assert channel.name == "test_channel"
        assert channel.index_address == "https://example.com/test_channel/index.yaml"
        assert channel.description == "a test channel"

        names = [c.name for c in refgenie_minimal.sources.list_channels()]
        assert len(names) == before + 1
        assert "test_channel" in names

    def test_get_channel_returns_none_when_absent(self, refgenie_minimal):
        """get_channel yields None for a name that was never added."""
        assert refgenie_minimal.sources.get_channel("nonexistent_12345") is None

    def test_remove_channel_deletes_and_reports_true(self, refgenie_minimal):
        """remove_channel returns True and the channel is gone afterward."""
        _add(refgenie_minimal, "removable")
        assert "removable" in [c.name for c in refgenie_minimal.sources.list_channels()]

        assert refgenie_minimal.sources.remove_channel("removable") is True
        assert "removable" not in [c.name for c in refgenie_minimal.sources.list_channels()]

    def test_remove_missing_channel_reports_false(self, refgenie_minimal):
        """Removing a channel that does not exist returns False, not an error."""
        assert refgenie_minimal.sources.remove_channel("nonexistent_12345") is False

    @pytest.mark.parametrize(
        "type_",
        [DataChannelType.http, DataChannelType.https, DataChannelType.local],
    )
    def test_add_channel_accepts_each_handler_type(self, refgenie_minimal, type_):
        """Each channel type maps to a handler and persists with that type."""
        channel = _add(refgenie_minimal, f"chan_{type_.value}", type=type_)
        assert channel.type == type_
        assert refgenie_minimal.sources.get_channel(f"chan_{type_.value}").type == type_


class TestDataChannelValidation:
    """Documented guards on add_channel."""

    def test_duplicate_name_raises_value_error(self, refgenie_minimal):
        """A second add with the same name raises ValueError (channel already exists)."""
        _add(refgenie_minimal, "dup")
        with pytest.raises(ValueError, match="already exists"):
            _add(refgenie_minimal, "dup")

    def test_double_underscore_in_name_raises_value_error(self, refgenie_minimal):
        """'__' is reserved as a channel/file separator and is rejected in names."""
        with pytest.raises(ValueError, match="__"):
            _add(refgenie_minimal, "bad__name")


HTML_RESPONSE = """<!DOCTYPE html>
<html>
<head><title>Recipes</title></head>
<body><h1>Refgenie Recipes</h1></body>
</html>
"""

VALID_INDEX_YAML = """
asset_class:
  dir: asset_classes
  files:
    - fasta.yaml
recipe:
  dir: recipes
  files:
    - fasta.yaml
"""


def _add_http_channel(rgc, name, url):
    """Helper to add an HTTP data channel."""
    return rgc.sources.add_channel(
        name=name,
        type=DataChannelType.http,
        index_address=url,
    )


class TestDataChannelIndexParsing:
    """Tests for index file parsing edge cases."""

    def test_html_response_gives_clear_error(self, refgenie_minimal):
        """
        When the channel URL returns HTML instead of YAML, get_index_file
        should return None (or raise a clear error), not crash with
        'argument after ** must be a mapping, not str'.
        """
        _add_http_channel(
            refgenie_minimal,
            "html-channel",
            "https://refgenie.github.io/refgenie-registry/",
        )

        with patch(
            "refgenie.managers.sources.handlers.HTTPChannelHandler.fetch_index_content",
            return_value=HTML_RESPONSE,
        ):
            # This should NOT raise TypeError about "argument after ** must be a mapping"
            result = refgenie_minimal.sources.get_index_file("html-channel")
            assert result is None

    def test_valid_yaml_index_parses_correctly(self, refgenie_minimal):
        """A well-formed index.yaml should parse into an IndexFile."""
        _add_http_channel(
            refgenie_minimal,
            "good-channel",
            "https://example.com/index.yaml",
        )

        with patch(
            "refgenie.managers.sources.handlers.HTTPChannelHandler.fetch_index_content",
            return_value=VALID_INDEX_YAML,
        ):
            result = refgenie_minimal.sources.get_index_file("good-channel")
            assert result is not None
            assert result.asset_class.dir == "asset_classes"
            assert [str(f) for f in result.recipe.files] == ["fasta.yaml"]


class _FakeRegistry:
    """An in-memory stand-in for the asset-class or recipe manager.

    ``add`` records the URL, raising ``exists_error`` for names already seen
    (unless overwriting) and ``fail_error`` for URLs listed in ``failing``.
    """

    def __init__(self, exists_error, existing=(), failing=()):
        self.exists_error = exists_error
        self.added: list[str] = list(existing)
        self.failing = set(failing)

    def add(self, url, exists_overwrite=False):
        from refgenie.exceptions import RefgenieError

        if url in self.failing:
            raise RefgenieError(f"cannot read {url}")
        if url in self.added and not exists_overwrite:
            raise self.exists_error(url, "1.0.0")
        self.added.append(url)


class TestSyncChannel:
    """``SourceManager.sync_channel`` with in-memory asset-class and recipe managers."""

    pytestmark = pytest.mark.unit

    @staticmethod
    def _sources(asset_classes, recipes, *, reachable=True, **fakes):
        from refgenie.exceptions import AssetClassExistsError, RecipeExistsError
        from refgenie.managers.sources import SourceManager
        from tests.helpers import make_engine

        asset_class = _FakeRegistry(
            AssetClassExistsError,
            existing=fakes.get("existing_asset_classes", ()),
            failing=fakes.get("failing_asset_classes", ()),
        )
        recipe = _FakeRegistry(RecipeExistsError, existing=fakes.get("existing_recipes", ()))
        sources = SourceManager(make_engine(), asset_class, recipe)
        sources.test_channel = lambda name: reachable
        sources.iter_asset_classes = lambda name: iter(asset_classes)
        sources.iter_recipes = lambda name: iter(recipes)
        return sources, asset_class, recipe

    def test_unreachable_channel_raises(self):
        from refgenie.exceptions import RefgenieError

        sources, asset_class, recipe = self._sources(["ac1"], ["r1"], reachable=False)
        with pytest.raises(RefgenieError, match="not accessible"):
            sources.sync_channel("chan")
        assert asset_class.added == [] and recipe.added == []

    def test_asset_class_failure_stops_the_recipes(self):
        sources, asset_class, recipe = self._sources(
            ["ac1", "ac2"], ["r1"], failing_asset_classes=["ac2"]
        )
        report = sources.sync_channel("chan")
        assert report.asset_classes_added == 1
        assert report.asset_classes_failed == 1
        assert recipe.added == []
        assert report.recipes_added == 0
        assert not report.ok

    def test_exists_ok_counts_skips(self):
        sources, _, recipe = self._sources(
            ["ac1", "ac2"],
            ["r1", "r2"],
            existing_asset_classes=["ac1"],
            existing_recipes=["r2"],
        )
        report = sources.sync_channel("chan", exists_ok=True)
        assert report.asset_classes_added == 1
        assert report.asset_classes_skipped == 1
        assert report.recipes_added == 1
        assert report.recipes_skipped == 1
        assert report.errors == []
        assert recipe.added == ["r2", "r1"]
