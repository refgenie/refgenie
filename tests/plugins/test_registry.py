"""Plugin discovery, loading, the off switch and describe()."""

import logging

from refgenie.plugins import registry
from tests.plugins.fake_plugin import RECORD


def test_load_returns_plugins_in_name_order(install_plugins):
    install_plugins({"post_update": [("zeta", RECORD), ("alpha", RECORD)]})
    assert [name for name, _ in registry.load("post_update")] == ["alpha", "zeta"]


def test_unknown_hook_is_described_but_never_loaded(install_plugins):
    install_plugins({"pre_tag": [("legacy", RECORD)]})
    assert registry.load("pre_tag") == []
    (info,) = registry.describe()
    assert (info.hook, info.name, info.status) == (
        "pre_tag",
        "legacy",
        registry.STATUS_UNKNOWN_HOOK,
    )


def test_plugin_that_does_not_import_is_skipped_and_reported(install_plugins, caplog):
    install_plugins({"post_update": [("broken", "no_such_module_xyz:func"), ("ok", RECORD)]})
    with caplog.at_level(logging.WARNING, logger="refgenie"):
        loaded = registry.load("post_update")
    assert [name for name, _ in loaded] == ["ok"]
    assert "refgenie plugin 'broken'" in caplog.text
    statuses = {i.name: i.status for i in registry.describe()}
    assert statuses["broken"].startswith("load error:")
    assert statuses["ok"] == registry.STATUS_OK


def test_disabled_by_name(install_plugins, disable_plugins):
    install_plugins({"post_update": [("one", RECORD), ("two", RECORD)]})
    disable_plugins("one")
    assert [name for name, _ in registry.load("post_update")] == ["two"]
    statuses = {i.name: i.status for i in registry.describe()}
    assert statuses == {"one": registry.STATUS_DISABLED, "two": registry.STATUS_OK}


def test_disabled_all(install_plugins, disable_plugins):
    install_plugins({"post_update": [("one", RECORD)]})
    for value in ("1", "true", "all", "ALL"):
        disable_plugins(value)
        assert registry.disabled() == "all"
        assert registry.load("post_update") == []


def test_zero_and_false_disable_nothing(install_plugins, disable_plugins):
    install_plugins({"post_update": [("one", RECORD)]})
    for value in ("", "0", "false"):
        disable_plugins(value)
        assert registry.disabled() == set()
        assert [name for name, _ in registry.load("post_update")] == ["one"]


def test_describe_reports_no_dist_for_a_bare_entry_point(install_plugins):
    install_plugins({"post_build": [("one", RECORD)]})
    (info,) = registry.describe()
    assert (info.target, info.dist, info.version) == (RECORD, None, None)
