"""PluginHost: errors are logged, not raised; plugins cannot re-trigger hooks."""

import logging

import pytest

from refgenie.plugins.hooks import PRE_PULL, HookEvent
from refgenie.plugins.host import PluginHost
from tests.plugins import fake_plugin
from tests.plugins.fake_plugin import RECORD


class StubRefgenie:
    """What a plugin receives as `rg`; `dispatch` lets a plugin re-enter."""

    dispatch = None


def _host():
    rg = StubRefgenie()
    host = PluginHost(rg)
    rg.dispatch = host
    return host


def test_a_failing_plugin_is_logged_and_the_next_still_runs(install_plugins, caplog):
    install_plugins(
        {PRE_PULL: [("a_boom", "tests.plugins.fake_plugin:boom"), ("b_record", RECORD)]}
    )
    with caplog.at_level(logging.WARNING, logger="refgenie"):
        _host()(HookEvent(hook=PRE_PULL))
    assert "refgenie plugin 'a_boom' failed in pre_pull" in caplog.text
    assert "REFGENIE_DISABLE_PLUGINS=a_boom" in caplog.text
    assert [hook for hook, _ in fake_plugin.CALLS] == [PRE_PULL]


def test_a_plugin_cannot_re_dispatch(install_plugins):
    install_plugins({PRE_PULL: [("reenter", "tests.plugins.fake_plugin:reenter")]})
    _host()(HookEvent(hook=PRE_PULL))
    assert [hook for hook, _ in fake_plugin.CALLS] == ["reenter"]


def test_keyboard_interrupt_propagates(install_plugins):
    install_plugins({PRE_PULL: [("stop", "tests.plugins.fake_plugin:interrupt")]})
    host = _host()
    with pytest.raises(KeyboardInterrupt):
        host(HookEvent(hook=PRE_PULL))
    # The guard is released, so the next dispatch runs again.
    install_plugins({PRE_PULL: [("record", RECORD)]})
    host(HookEvent(hook=PRE_PULL))
    assert [hook for hook, _ in fake_plugin.CALLS] == [PRE_PULL]


def test_no_plugins_means_nothing_runs(install_plugins):
    install_plugins({})
    _host()(HookEvent(hook=PRE_PULL))
    assert fake_plugin.CALLS == []
