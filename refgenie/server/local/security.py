"""Local-mode security: origin allowlist, host guard, action header.

The local app manages one user's genomes over loopback HTTP, which makes it a
CSRF and DNS-rebinding target: any web page the user visits can ``fetch()``
``http://127.0.0.1:<port>`` unless something here says no. Three mechanisms,
installed together by :func:`install_local_security` because none is
sufficient alone:

* **CORS allowlist** -- only the deployed SPA and the Vite dev server may make
  cross-origin requests, without credentials, and never with ``DELETE``.
* **Custom header** (``X-Refgenie-Action``) -- required on every state-changing
  route. It is not a CORS-safelisted header, so a cross-origin request that
  carries it is *forced* into a preflight, where the allowlist decides. A
  malicious page cannot mount a simple-request CSRF write.
* **Host guard** -- the reason the header is not enough. A DNS-rebinding
  attacker resolves ``evil.example`` to ``127.0.0.1``; the browser then treats
  the local server as same-origin, and CORS plus the preflight-forcing header
  no longer apply. Rejecting any request whose ``Host`` is not a loopback name
  (421 Misdirected Request) closes that hole. This is the same class of guard
  ``mcp``'s ``TransportSecuritySettings`` provides and which the *public*
  server deliberately disables -- local mode is the case where it must be on.

A fourth concern, the **localhost bridge**, extends the same settings object
rather than adding another middleware path: the deployed public SPA
(``https://refgenie.org``) may probe and read a local instance when
``bridge_mode`` is ``read`` (the default), and may additionally submit a pull
when it is ``full``. The bridge contributes origins to the one CORS allowlist,
one extra middleware (the Local Network Access preflight header), and one
extra dependency on the actions router (:func:`require_action_origin`). See
``REFGENIE_BRIDGE_*`` in ``refgenie/config``.

This module is the single owner of local-mode security. ``create_app`` calls
:func:`install_local_security` and adds no local CORS, host guard or header
dependency of its own.
"""

import re
from typing import Annotated, Literal

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from pydantic import AliasChoices, Field, field_validator
from pydantic_settings import BaseSettings, NoDecode, SettingsConfigDict

from refgenie.config import (
    REFGENIE_BRIDGE_MODE,
    REFGENIE_BRIDGE_ORIGIN_REGEX,
    REFGENIE_BRIDGE_ORIGINS,
)
from refgenie.server.errors import ErrorCode, error_response

__all__ = [
    "BRIDGE_CROSS_ORIGIN_ACTIONS",
    "HostHeaderGuardMiddleware",
    "LocalNetworkAccessMiddleware",
    "LocalSecuritySettings",
    "install_local_security",
    "require_action_header",
    "require_action_origin",
]

#: The only actions reachable from an allowlisted *bridge* origin, and only
#: under ``bridge_mode == "full"``. Everything else (delete, alias, subscribe,
#: build, ...) is same-origin only -- enforced here by an explicit allowed-path
#: set, not by hoping the remote UI never calls them. Pull is the one action
#: whose natural trigger lives on the remote page; keeping destructive verbs
#: off the cross-origin surface caps the blast radius of a misconfigured
#: bridge at disclosure plus an unwanted download.
BRIDGE_CROSS_ORIGIN_ACTIONS: frozenset[tuple[str, str]] = frozenset({("POST", "/v1/actions/pull")})

#: Regexes that would match (nearly) everything -- refused at construction,
#: same as a literal "*" origin. Best-effort: a user can still write a
#: wildcard-equivalent regex we cannot recognize; that is documented as a
#: deliberately-not-defended-against misconfiguration.
_WILDCARD_REGEXES = frozenset({".*", "^.*$", ".+", "^.+$", ".*$", "^.*"})


class LocalSecuritySettings(BaseSettings):
    """Local-mode security knobs, overridable via ``REFGENIE_LOCAL_*`` env vars.

    The ``bridge_*`` fields are the localhost bridge's configuration and read
    the unprefixed ``REFGENIE_BRIDGE_*`` env vars (their defaults come from
    ``refgenie.config``, which reads the same variables at import). Bridge
    origins are *public* origins: they get read access under ``read`` mode and
    exactly one action (``POST /v1/actions/pull``) under ``full`` mode --
    unlike ``allowed_origins``, whose members (the developer's own Vite dev
    server) are trusted like the same-origin SPA.

    ``require_action_header`` exists for tests and is not a supported
    production configuration.
    """

    model_config = SettingsConfigDict(env_prefix="REFGENIE_LOCAL_", populate_by_name=True)

    #: Fully trusted cross-origin callers: the Vite dev server serving the
    #: *local* SPA during development. ``https://refgenie.org`` is NOT here
    #: any more -- the deployed public SPA is a *bridge* origin (below) with a
    #: deliberately smaller surface. ``https://docs.refgenie.org`` is the docs
    #: site and is allowlisted nowhere.
    allowed_origins: list[str] = [
        "http://localhost:5173",
        "http://127.0.0.1:5173",
    ]
    #: Host-header names accepted by the guard, port-insensitive.
    allowed_hosts: list[str] = ["127.0.0.1", "localhost", "::1"]
    #: The header whose *presence* marks a deliberate state-changing request.
    action_header: str = "x-refgenie-action"
    require_action_header: bool = True

    # --- Localhost bridge (REFGENIE_BRIDGE_*) -------------------------------

    #: ``off`` = bridge origins get no CORS at all; ``read`` (default) = they
    #: may read (/ping, /v4, /v1/jobs); ``full`` = additionally
    #: ``POST /v1/actions/pull``.
    bridge_mode: Literal["off", "read", "full"] = Field(
        default=REFGENIE_BRIDGE_MODE,
        validation_alias=AliasChoices("bridge_mode", "REFGENIE_BRIDGE_MODE"),
    )
    #: Exact public origins the bridge admits (comma-separated in the env var).
    bridge_origins: Annotated[list[str], NoDecode] = Field(
        default=[o.strip() for o in REFGENIE_BRIDGE_ORIGINS.split(",") if o.strip()],
        validation_alias=AliasChoices("bridge_origins", "REFGENIE_BRIDGE_ORIGINS"),
    )
    #: Optional regex over origins, for SPA dev servers and Cloudflare preview
    #: deployments only. A careless regex here undoes every bridge protection.
    bridge_origin_regex: str | None = Field(
        default=REFGENIE_BRIDGE_ORIGIN_REGEX or None,
        validation_alias=AliasChoices("bridge_origin_regex", "REFGENIE_BRIDGE_ORIGIN_REGEX"),
    )
    #: Whether /ping reveals filesystem paths (default: it must not -- a home
    #: directory path leaks the OS username to every allowlisted origin).
    bridge_expose_paths: bool = Field(
        default=False,
        validation_alias=AliasChoices("bridge_expose_paths", "REFGENIE_BRIDGE_EXPOSE_PATHS"),
    )

    @field_validator("bridge_origins", mode="before")
    @classmethod
    def _split_comma_separated(cls, value):
        if isinstance(value, str):
            return [origin.strip() for origin in value.split(",") if origin.strip()]
        return value

    @field_validator("bridge_origin_regex", mode="before")
    @classmethod
    def _empty_regex_is_none(cls, value):
        if isinstance(value, str) and not value.strip():
            return None
        return value

    def cross_origin_allowlist(self) -> list[str]:
        """Every exact origin CORS admits: dev origins plus, when the bridge
        is on, the bridge origins."""
        origins = list(self.allowed_origins)
        if self.bridge_mode != "off":
            origins += [o for o in self.bridge_origins if o not in origins]
        return origins

    def is_bridge_origin(self, origin: str) -> bool:
        """Whether ``origin`` is a *public* bridge origin (exact or regex).

        Deliberately independent of ``bridge_mode``: classification does not
        change with the mode, only what a bridge origin is allowed to do.
        """
        if origin in self.bridge_origins:
            return True
        if self.bridge_origin_regex is not None:
            return re.fullmatch(self.bridge_origin_regex, origin) is not None
        return False

    def origin_is_allowlisted(self, origin: str) -> bool:
        """Whether ``origin`` may make cross-origin requests at all right now
        (the same set CORS grants headers to)."""
        if origin in self.allowed_origins:
            return True
        if self.bridge_mode == "off":
            return False
        return self.is_bridge_origin(origin)


def require_action_header(request: Request) -> None:
    """403 unless the request carries the custom action header.

    A dependency rather than middleware on purpose: attached once via
    ``APIRouter(dependencies=[Depends(require_action_header)])`` it covers
    every route on the actions router by construction (a newly added route
    cannot forget it), it never intercepts CORS preflights (``CORSMiddleware``
    answers ``OPTIONS`` before routing), and it is unit-testable without an
    app. The header's presence is what matters, never its value; clients send
    ``X-Refgenie-Action: 1``.
    """
    settings = getattr(request.app.state, "local_security", None)
    if settings is None:
        settings = LocalSecuritySettings()
    if not settings.require_action_header:
        return
    if settings.action_header not in request.headers:
        raise HTTPException(
            status_code=403,
            detail={
                "code": str(ErrorCode.MISSING_ACTION_HEADER),
                "message": (
                    "State-changing requests must carry the "
                    f"'{settings.action_header}' header (any value)."
                ),
            },
        )


def require_action_origin(request: Request) -> None:
    """403 any cross-origin action the bridge policy does not permit.

    Runs on the actions router alongside :func:`require_action_header` (which
    is the load-bearing anti-CSRF control: the custom header forces a
    preflight, so no simple cross-origin request reaches an action at all).
    This dependency is the *policy* layer on top:

    * no ``Origin``, or ``Origin`` == the request's own origin -- the
      same-origin local SPA; allow.
    * a **bridge** origin (the deployed public SPA) -- allow only under
      ``bridge_mode == "full"``, and only for the explicit
      :data:`BRIDGE_CROSS_ORIGIN_ACTIONS` set; the 403 message names the exact
      remedy so the SPA can render it verbatim.
    * a plain allowlisted origin (the Vite dev server) -- allow; it is the
      developer's own machine serving the local SPA.
    * anything else -- 403. CORS already stops the *response* from being read,
      but the request would still execute server-side without this check.
    """
    origin = request.headers.get("origin")
    if origin is None:
        return
    own_origin = f"{request.url.scheme}://{request.headers.get('host', '')}"
    if origin == own_origin:
        return
    settings = getattr(request.app.state, "local_security", None)
    if settings is None:
        settings = LocalSecuritySettings()
    if settings.is_bridge_origin(origin):
        if settings.bridge_mode != "full":
            raise HTTPException(
                status_code=403,
                detail={
                    "code": str(ErrorCode.FORBIDDEN_ORIGIN),
                    "message": (
                        "Cross-origin actions are disabled on this refgenie "
                        f"(bridge mode is '{settings.bridge_mode}'). "
                        "Restart with `refgenie dash --bridge full` to enable them."
                    ),
                },
            )
        path = request.scope["path"]
        root_path = request.scope.get("root_path", "")
        if root_path and path.startswith(root_path):
            path = path[len(root_path) :]
        if (request.method, path) not in BRIDGE_CROSS_ORIGIN_ACTIONS:
            raise HTTPException(
                status_code=403,
                detail={
                    "code": str(ErrorCode.FORBIDDEN_ORIGIN),
                    "message": (
                        "This action is not available cross-origin. Open your "
                        "local refgenie directly (http://127.0.0.1:<port>) to "
                        "perform it."
                    ),
                },
            )
        return
    if origin in settings.allowed_origins:
        return
    raise HTTPException(
        status_code=403,
        detail={
            "code": str(ErrorCode.FORBIDDEN_ORIGIN),
            "message": f"Origin {origin!r} is not allowed to perform actions.",
        },
    )


def _hostname(host_header: str) -> str:
    """The hostname part of a ``Host`` header value, port stripped.

    ``127.0.0.1:8080`` -> ``127.0.0.1``; ``[::1]:8080`` -> ``::1``.
    """
    host_header = host_header.strip()
    if host_header.startswith("["):
        end = host_header.find("]")
        return host_header[1:end] if end != -1 else host_header
    return host_header.split(":", 1)[0]


class HostHeaderGuardMiddleware:
    """Reject requests whose ``Host`` is not a loopback name -- 421, envelope code
    ``forbidden_host``.

    Pure ASGI, and deliberately not Starlette's ``TrustedHostMiddleware``:
    that one strips the port before matching and supports patterns this design
    does not want. This is the single host-guard implementation in the
    codebase.
    """

    def __init__(self, app, allowed_hosts: list[str]):
        self.app = app
        self.allowed = {host.lower() for host in allowed_hosts}

    async def __call__(self, scope, receive, send):
        if scope["type"] not in ("http", "websocket"):
            await self.app(scope, receive, send)
            return
        host = ""
        for name, value in scope.get("headers", []):
            if name == b"host":
                host = value.decode("latin-1")
                break
        if _hostname(host).lower() not in self.allowed:
            response = error_response(
                421,
                ErrorCode.FORBIDDEN_HOST,
                f"Host {host!r} is not a loopback name; refusing (DNS-rebinding guard).",
            )
            await response(scope, receive, send)
            return
        await self.app(scope, receive, send)


#: The request-header spellings a browser may use for the private/local
#: network preflight, mapped to the response header each expects. Chrome's
#: original Private Network Access spec used "-Private-Network"; the LNA-era
#: successor renames it to "-Local-Network". Feature-detected per request
#: rather than hardcoding one spelling.
_LNA_HEADER_PAIRS: tuple[tuple[bytes, bytes], ...] = (
    (b"access-control-request-private-network", b"access-control-allow-private-network"),
    (b"access-control-request-local-network", b"access-control-allow-local-network"),
)


class LocalNetworkAccessMiddleware:
    """Answer the private/local-network CORS preflight for allowlisted origins.

    Starlette's ``CORSMiddleware`` handles the PNA-era
    ``Access-Control-Request-Private-Network`` spelling itself (via
    ``allow_private_network=True``, set in :func:`install_local_security`; the
    flag is required -- without it Starlette *400s* the preflight before this
    middleware could help). What Starlette does not know is the LNA-era
    ``-Local-Network`` successor spelling; this middleware appends whichever
    grant header is still missing -- **only** for origins the bridge allowlist
    admits; emitting it unconditionally would hand a blanket local-network
    grant to any caller.

    Compatibility note: this is the belt-and-braces path for Chrome's older
    PNA behavior. In Chrome 142+ the gate is a one-time *user permission
    prompt* that no server header can satisfy; the header is harmless there.
    Revisit when the LNA spec settles.

    Pure ASGI (matching :class:`HostHeaderGuardMiddleware`), installed by
    :func:`install_local_security` when ``bridge_mode != "off"``, outermost so
    it sees the preflight response ``CORSMiddleware`` produces.
    """

    def __init__(self, app, settings: LocalSecuritySettings):
        self.app = app
        self.settings = settings

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http" or scope.get("method") != "OPTIONS":
            await self.app(scope, receive, send)
            return
        origin = ""
        lna_response_headers: list[bytes] = []
        for name, value in scope.get("headers", []):
            if name == b"origin":
                origin = value.decode("latin-1")
            for request_header, response_header in _LNA_HEADER_PAIRS:
                if name == request_header and value.strip().lower() == b"true":
                    lna_response_headers.append(response_header)
        if not lna_response_headers or not (origin and self.settings.origin_is_allowlisted(origin)):
            await self.app(scope, receive, send)
            return

        async def send_with_lna_grant(message):
            if message["type"] == "http.response.start":
                headers = list(message.get("headers", []))
                # CORSMiddleware (allow_private_network=True) already answers
                # the "-Private-Network" spelling; only add what is missing so
                # the response never carries a duplicated grant header.
                present = {name.lower() for name, _ in headers}
                headers += [
                    (header, b"true")
                    for header in lna_response_headers
                    if header not in present
                ]
                message = {**message, "headers": headers}
            await send(message)

        await self.app(scope, receive, send_with_lna_grant)


def install_local_security(app: FastAPI, settings: LocalSecuritySettings | None = None) -> None:
    """Install the host guard, the local CORS policy and the bridge layers.

    Order matters: Starlette's ``add_middleware`` prepends, so the **last**
    added runs outermost. CORS is added after the host guard so it answers
    preflights before the guard sees them; the LNA middleware is added after
    CORS (outermost) so it can decorate the preflight *response* CORS builds.

    ``allow_credentials=False`` is deliberate -- local mode has no cookies or
    auth, so an allowlisted origin can never ride ambient authority. ``DELETE``
    is deliberately absent from ``allow_methods``: destructive verbs are not on
    the cross-origin surface at all (same-origin callers are unaffected --
    CORS does not apply to them).

    The bridge contributes origins to this one CORS policy (under ``read`` or
    ``full`` mode); under ``off`` it contributes nothing, so only the
    same-origin/dev defaults remain and a public origin gets no CORS grant.

    Raises:
        ValueError: If any origin allowlist contains ``"*"`` (or the regex
            escape hatch is a match-everything regex) -- a wildcard would
            defeat the entire design, so it fails at construction, not at
            request time.
    """
    settings = settings or LocalSecuritySettings()
    if "*" in settings.allowed_origins or "*" in settings.bridge_origins:
        raise ValueError(
            "allow_origins=['*'] is not permitted in local mode: the origin "
            "allowlist is the CSRF boundary. List explicit origins."
        )
    if settings.bridge_origin_regex is not None:
        re.compile(settings.bridge_origin_regex)  # fail fast on a bad regex
        if settings.bridge_origin_regex.strip() in _WILDCARD_REGEXES:
            raise ValueError(
                "REFGENIE_BRIDGE_ORIGIN_REGEX matches every origin, which "
                "defeats the origin allowlist entirely. List explicit origins "
                "or narrow the regex."
            )
    app.state.local_security = settings
    bridge_on = settings.bridge_mode != "off"
    app.add_middleware(HostHeaderGuardMiddleware, allowed_hosts=settings.allowed_hosts)
    app.add_middleware(
        CORSMiddleware,
        allow_origins=settings.cross_origin_allowlist(),
        allow_origin_regex=settings.bridge_origin_regex if bridge_on else None,
        allow_credentials=False,
        allow_methods=["GET", "POST", "OPTIONS"],  # no DELETE cross-origin
        allow_headers=["content-type", settings.action_header],
        # Chrome's Private Network Access preflight carries
        # Access-Control-Request-Private-Network, and Starlette 400s that
        # preflight outright unless this flag is set -- our LNA middleware
        # decorating the response afterwards cannot un-400 it. Only when the
        # bridge is on: the origin allowlist above still refuses everyone
        # else, so the broader echo is harmless.
        allow_private_network=bridge_on,
        max_age=600,
    )
    if bridge_on:
        app.add_middleware(LocalNetworkAccessMiddleware, settings=settings)
