"""Local-mode security: origin allowlist, host guard, action header.

The local app listens on loopback, so any web page could ``fetch()`` it.
:func:`install_local_security` installs the defenses together: a CORS
allowlist, a required ``X-Refgenie-Action`` header on state-changing routes
(forces a preflight), a loopback-only ``Host`` guard against DNS rebinding,
and the localhost bridge (``REFGENIE_BRIDGE_*``) that lets the public SPA
read, and optionally pull. :func:`require_local_origin` keeps ``/v1/files``
off the bridge.
"""

import re
from typing import Annotated, Literal

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from pydantic import AliasChoices, Field, field_validator
from pydantic_settings import BaseSettings, NoDecode, SettingsConfigDict

from refgenie.server.errors import ErrorCode, error_response

__all__ = [
    "BRIDGE_CROSS_ORIGIN_ACTIONS",
    "HostHeaderGuardMiddleware",
    "LocalNetworkAccessMiddleware",
    "LocalSecuritySettings",
    "install_local_security",
    "is_bridge_request",
    "require_action_header",
    "require_action_origin",
    "require_local_origin",
]

#: The only actions a bridge origin may call, and only under ``bridge_mode ==
#: "full"``. Destructive actions stay same-origin.
BRIDGE_CROSS_ORIGIN_ACTIONS: frozenset[tuple[str, str]] = frozenset({("POST", "/v1/actions/pull")})

#: Match-everything regexes, refused like a literal "*" origin (best effort).
_WILDCARD_REGEXES = frozenset({".*", "^.*$", ".+", "^.+$", ".*$", "^.*"})


class LocalSecuritySettings(BaseSettings):
    """Local-mode security settings, from ``REFGENIE_LOCAL_*`` env vars.

    The ``bridge_*`` fields read ``REFGENIE_BRIDGE_*``. Bridge origins are
    public: read access under ``read``, plus pull under ``full``.
    ``require_action_header=False`` is for tests only.
    """

    model_config = SettingsConfigDict(env_prefix="REFGENIE_LOCAL_", populate_by_name=True)

    #: Fully trusted cross-origin callers: the Vite dev server for the local SPA.
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

    #: ``off``: no CORS for bridge origins; ``read``: /ping, /v4, /v1/jobs;
    #: ``full``: also ``POST /v1/actions/pull``.
    bridge_mode: Literal["off", "read", "full"] = Field(
        default="read",
        validation_alias=AliasChoices("bridge_mode", "REFGENIE_BRIDGE_MODE"),
    )
    #: Exact public origins the bridge admits (comma-separated in the env var).
    bridge_origins: Annotated[list[str], NoDecode] = Field(
        default=["https://refgenie.org", "https://api.refgenie.org"],
        validation_alias=AliasChoices("bridge_origins", "REFGENIE_BRIDGE_ORIGINS"),
    )
    #: Optional origin regex, for SPA dev servers and preview deployments only.
    bridge_origin_regex: str | None = Field(
        default=None,
        validation_alias=AliasChoices("bridge_origin_regex", "REFGENIE_BRIDGE_ORIGIN_REGEX"),
    )
    #: Whether /ping reveals filesystem paths (off: a home path leaks the username).
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
        """Whether ``origin`` is a public bridge origin (exact or regex), whatever the mode."""
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

    A router dependency, so it covers every actions route and never sees
    CORS preflights. Only the header's presence matters; clients send
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


def _request_settings(request: Request) -> LocalSecuritySettings:
    """The app's security settings, or fresh ones when the app installed none."""
    settings = getattr(request.app.state, "local_security", None)
    return settings if settings is not None else LocalSecuritySettings()


def _is_same_origin(request: Request, origin: str | None) -> bool:
    """No ``Origin`` header, or the request's own origin: the local SPA."""
    if origin is None:
        return True
    return origin == f"{request.url.scheme}://{request.headers.get('host', '')}"


def is_bridge_request(request: Request) -> bool:
    """Whether the request comes from a public bridge origin (not same-origin, not a trusted dev origin)."""
    origin = request.headers.get("origin")
    if _is_same_origin(request, origin):
        return False
    return _request_settings(request).is_bridge_origin(origin)


def require_action_origin(request: Request) -> None:
    """403 any cross-origin action the bridge policy does not permit.

    Same-origin and ``allowed_origins`` callers pass. Bridge origins pass only
    under ``bridge_mode == "full"`` and only for
    :data:`BRIDGE_CROSS_ORIGIN_ACTIONS`. Everything else gets 403, since CORS
    blocks reading the response but not running the request.
    """
    origin = request.headers.get("origin")
    if _is_same_origin(request, origin):
        return
    settings = _request_settings(request)
    if is_bridge_request(request):
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


#: ``Sec-Fetch-Site`` values that mean "another site started this request".
_CROSS_SITE_FETCH = frozenset({"cross-site", "same-site"})


def require_local_origin(request: Request) -> None:
    """403 unless the request comes from the local SPA or a trusted dev origin.

    Guards routes that serve bytes from the user's disk. Bridge origins are
    refused under every ``bridge_mode``. With no ``Origin`` header,
    ``Sec-Fetch-Site: cross-site``/``same-site`` (an embed on another page)
    is refused; anything else, including curl, is allowed.
    """
    origin = request.headers.get("origin")
    settings = _request_settings(request)
    if origin is None:
        fetch_site = request.headers.get("sec-fetch-site", "").strip().lower()
        if fetch_site not in _CROSS_SITE_FETCH:
            return
    elif _is_same_origin(request, origin) or origin in settings.allowed_origins:
        return
    raise HTTPException(
        status_code=403,
        detail={
            "code": str(ErrorCode.FORBIDDEN_ORIGIN),
            "message": (
                "Local files are readable only from this refgenie's own web UI. "
                "Open your local refgenie directly (http://127.0.0.1:<port>)."
            ),
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
    """Reject requests whose ``Host`` is not a loopback name (421 ``forbidden_host``).

    Not Starlette's ``TrustedHostMiddleware``, which supports patterns this
    guard does not want.
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


#: Request header for the private/local-network preflight, mapped to the
#: response header that grants it. Browsers use either spelling.
_LNA_HEADER_PAIRS: tuple[tuple[bytes, bytes], ...] = (
    (b"access-control-request-private-network", b"access-control-allow-private-network"),
    (b"access-control-request-local-network", b"access-control-allow-local-network"),
)


class LocalNetworkAccessMiddleware:
    """Answer the private/local-network CORS preflight for allowlisted origins.

    Starlette handles the older ``-Private-Network`` spelling; this adds the
    newer ``-Local-Network`` grant, only for origins the allowlist admits.
    Installed outermost when ``bridge_mode != "off"``.
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
                # Add only the grant headers CORSMiddleware did not already set.
                present = {name.lower() for name, _ in headers}
                headers += [
                    (header, b"true") for header in lna_response_headers if header not in present
                ]
                message = {**message, "headers": headers}
            await send(message)

        await self.app(scope, receive, send_with_lna_grant)


def install_local_security(app: FastAPI, settings: LocalSecuritySettings | None = None) -> None:
    """Install the host guard, the local CORS policy and the bridge layers.

    ``add_middleware`` prepends, so the last one added runs outermost: CORS
    answers preflights before the host guard, and the LNA middleware wraps
    the preflight response CORS builds. No credentials (local mode has no
    auth) and no cross-origin ``DELETE``.

    Raises:
        ValueError: If an origin allowlist contains ``"*"`` or the origin
            regex matches everything.
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
        # Starlette 400s a Private Network Access preflight unless this is set.
        allow_private_network=bridge_on,
        max_age=600,
    )
    if bridge_on:
        app.add_middleware(LocalNetworkAccessMiddleware, settings=settings)
