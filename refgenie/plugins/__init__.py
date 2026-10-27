"""refgenie's plugin system: outside packages that react to refgenie events.

A plugin is a function ``hook(rg, event)`` registered in an entry-point group
named ``refgenie.hooks.<hook>``, where ``<hook>`` is one of ``HOOKS``. ``rg`` is
the ``Refgenie`` instance and ``event`` is a ``HookEvent``.

How it fits together:

- ``hooks.py``: the hook names and the ``HookEvent`` / ``Change`` types.
- ``events.py``: ``EventSink`` and ``@update_scope``. Managers record into the
  sink; they never load or call plugins.
- ``registry.py``: entry-point discovery, loading, the off switch, and
  ``describe()`` for ``refgenie plugins``.
- ``host.py``: ``PluginHost``, the one place plugin code is called.

This package sits at layer 2 (see ``tests/test_layering.py``): it imports only
``config``, ``logger`` and ``const``, so managers and the root can both use it.
Discovery is lazy; importing this package scans nothing.

Only the plugin-author types are re-exported here, so ``from refgenie.plugins
import HookEvent`` stays cheap.
"""

from refgenie.plugins.hooks import HOOKS, Change, HookEvent

__all__ = ["HOOKS", "Change", "HookEvent"]
