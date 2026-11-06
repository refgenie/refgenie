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

It imports only ``config``, ``logger`` and ``const`` (see
``tests/test_layering.py``), and discovery is lazy.
"""

from refgenie.plugins.hooks import HOOKS, Change, HookEvent

__all__ = ["HOOKS", "Change", "HookEvent"]
