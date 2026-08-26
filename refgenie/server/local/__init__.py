"""The local-mode command surface: actions API, security, error envelope.

This package exists only in ``create_app(mode="local")`` (``refgenie dash``).
It is never included in server mode -- several request models accept
server-local filesystem paths, which is safe only on a loopback-bound,
header-guarded, single-user app.

Endpoint contract (v1)
----------------------

Router prefix ``/actions``, included once at ``/v1``. No GET routes: reads are
``/v4``, job polling is ``/v1/jobs``.

=======  ===================================  ======================  ==========================
Method   Path                                 Body model              Result
=======  ===================================  ======================  ==========================
POST     /v1/actions/pull                     PullRequest             202 JobRef
POST     /v1/actions/build                    BuildRequest            202 JobRef
POST     /v1/actions/build/preflight          BuildRequest            200 PreflightResult
POST     /v1/actions/genomes                  GenomeInitRequest       202 JobRef (genome init)
DELETE   /v1/actions/assets/{asset_digest}    --                      200 ActionResult
DELETE   /v1/actions/genomes/{genome_ref}     --                      200 ActionResult
POST     /v1/actions/aliases                  AliasSetRequest         200 ActionResult
DELETE   /v1/actions/aliases/{alias_name}     --                      200 ActionResult
POST     /v1/actions/subscriptions            SubscribeRequest        200 ActionResult
DELETE   /v1/actions/subscriptions            UnsubscribeRequest      200 ActionResult
POST     /v1/actions/assets/default           SetDefaultAssetRequest  200 ActionResult
=======  ===================================  ======================  ==========================

``recipe.add`` / ``asset_class.add`` are deliberately NOT here: they take a
server-side path or URL from an HTTP body (arbitrary local-file read + SSRF).
When a UI needs them, the endpoint must accept an uploaded YAML document.

The ``X-Refgenie-Action`` requirement
-------------------------------------

Every route above 403s (code ``missing_action_header``) unless the request
carries the ``X-Refgenie-Action`` header. Presence is the signal; the value is
ignored; clients send ``1``. Because it is not a CORS-safelisted header, any
cross-origin request carrying it is forced into a preflight, where the origin
allowlist in :mod:`refgenie.server.local.security` decides. The host guard
(421 ``forbidden_host`` for non-loopback ``Host`` headers) closes the
DNS-rebinding hole that CORS cannot.

Envelopes
---------

* 200: ``ActionResult`` -- ``{"ok": true, "message": str, "data": {...}|null}``.
* 202: the jobs package's ``JobRef`` --
  ``{job_id, kind, status, created_at, duplicate, links:{self, events, cancel}}``.
  ``duplicate: true`` means an identical in-flight job was coalesced (not an
  error); poll ``links.self`` or subscribe to ``/v1/jobs/events``.
* non-2xx: ``{"ok": false, "error": {"code", "message", "detail"?, "field"?}}``.
  Codes are :class:`refgenie.server.errors.ErrorCode`; the same code appears on
  a failed job's ``error.code``, so a failure means the same thing on both
  transports. Highlights: ``asset_exists`` 409, ``no_archive`` 404,
  ``pull_failed`` 502, ``no_subscriptions`` 409, ``build_failed`` 500,
  ``*_not_found`` 404, ``conflict`` 409, ``missing_action_header`` 403,
  ``forbidden_host`` 421, ``validation_error`` 422.
"""

from refgenie.server.errors import install_error_handlers
from refgenie.server.local.actions import router
from refgenie.server.local.security import (
    HostHeaderGuardMiddleware,
    LocalSecuritySettings,
    install_local_security,
    require_action_header,
)

#: The name ``create_app`` imports.
actions_router = router

__all__ = [
    "HostHeaderGuardMiddleware",
    "LocalSecuritySettings",
    "actions_router",
    "install_error_handlers",
    "install_local_security",
    "require_action_header",
    "router",
]
