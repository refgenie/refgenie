"""PluginSettingsManager (``rg.plugins``): per-plugin settings in the database."""

import pytest

from refgenie import Refgenie
from refgenie.plugins import HookEvent
from refgenie.plugins.hooks import POST_UPDATE
from tests.plugins import fake_plugin


def test_fresh_database_has_no_settings(refgenie_minimal):
    assert refgenie_minimal.plugins.settings("nfcore") == {}
    assert refgenie_minimal.plugins.all_settings() == {}


def test_set_round_trips_and_merges(refgenie_minimal):
    plugins = refgenie_minimal.plugins
    assert plugins.set("nfcore", config_path="/x/nf.config") == {"config_path": "/x/nf.config"}
    assert plugins.set("nfcore", other=3) == {"config_path": "/x/nf.config", "other": "3"}
    assert plugins.settings("nfcore") == {"config_path": "/x/nf.config", "other": "3"}


def test_settings_returns_a_copy(refgenie_minimal):
    plugins = refgenie_minimal.plugins
    plugins.set("nfcore", config_path="/x")
    plugins.settings("nfcore")["config_path"] = "/changed"
    assert plugins.settings("nfcore") == {"config_path": "/x"}


def test_unset_one_key_keeps_the_rest(refgenie_minimal):
    plugins = refgenie_minimal.plugins
    plugins.set("nfcore", config_path="/x", other="y")
    plugins.set("myplugin", greeting="hi")
    assert plugins.unset("nfcore", "other") == {"config_path": "/x"}
    assert plugins.all_settings() == {
        "nfcore": {"config_path": "/x"},
        "myplugin": {"greeting": "hi"},
    }


def test_unset_with_no_keys_removes_the_entry(refgenie_minimal):
    plugins = refgenie_minimal.plugins
    plugins.set("nfcore", config_path="/x")
    plugins.set("myplugin", greeting="hi")
    assert plugins.unset("nfcore") == {}
    assert plugins.all_settings() == {"myplugin": {"greeting": "hi"}}


def test_changes_are_committed(refgenie_minimal):
    refgenie_minimal.plugins.set("nfcore", config_path="/x")
    other = Refgenie(database_engine=refgenie_minimal.database_engine, suppress_migrations=True)
    assert other.plugins.settings("nfcore") == {"config_path": "/x"}


@pytest.mark.parametrize(
    "plugin, key", [("bad name", "k"), ("nfcore", "bad key"), ("", "k"), ("a/b", "k")]
)
def test_invalid_names_are_refused(refgenie_minimal, plugin, key):
    with pytest.raises(ValueError):
        refgenie_minimal.plugins.set(plugin, **{key: "v"})


def test_rerunning_init_keeps_the_settings(refgenie_minimal, tmp_path):
    refgenie_minimal.plugins.set("nfcore", config_path="/x")
    refgenie_minimal.database.init(genome_folder=tmp_path / "genomes")
    assert refgenie_minimal.plugins.settings("nfcore") == {"config_path": "/x"}


def test_a_plugin_reads_its_settings_inside_its_hook(refgenie_minimal, install_plugins):
    install_plugins({POST_UPDATE: [("fake", "tests.plugins.fake_plugin:read_settings")]})
    refgenie_minimal.plugins.set("fake", greeting="hi")
    # Dispatch through the instance's own sink, as a manager's record would.
    refgenie_minimal._events.emit(HookEvent(hook=POST_UPDATE))
    assert fake_plugin.SEEN_SETTINGS == [{"greeting": "hi"}]


def test_enabled_follows_the_instance(refgenie_minimal):
    assert refgenie_minimal.plugins.enabled is True
    off = Refgenie(
        database_engine=refgenie_minimal.database_engine, suppress_migrations=True, plugins=False
    )
    assert off.plugins.enabled is False
