"""The local-mode command surface: actions API, security, error envelope.

Mounted only by ``create_app(mode="local")`` (``refgenie dash``), never in
server mode: some request models accept local filesystem paths, which is safe
only on a loopback-bound, header-guarded, single-user app. Routes live under
``/v1/actions`` and are listed in the OpenAPI docs at ``/docs``.

``recipe.add`` and ``asset_class.add`` are not exposed: they take a path or
URL from the request body (local file read and SSRF).

Import from the submodules; this package imports nothing, so
``refgenie.server.jobs`` can use ``local.security`` without a cycle.
"""
