"""Serving the built web UI (the React SPA) from the refgenie app.

The UI source lives in ``frontend/`` at the repo root and vite builds it
**directly into the package** (``build.outDir = ../refgenie/server/webui``,
``build.assetsDir = _app``, ``base = /``). There is therefore one bundle
location, not a packaged one and a dev one: ``refgenie/server/webui/``.

There is deliberately no ``refgenie/server/web.py`` and no
``refgenie/server/web/`` directory -- a module and a package of the same name
shadow each other on import. Everything that serves the bundle lives here.

**root_path.** ``frontend/index.html`` ships a literal ``<base href="/">``.
Mounting behind a reverse-proxy prefix (``create_app(root_path="/refgenie")``)
rewrites that one token at app-construction time and caches the rendered
string, so a sub-path deployment needs no rebuild. The SPA takes its router
basename and API base from ``/service-info``, not from a build-time constant.

**A missing bundle is never fatal.** A source checkout that has not run
``npm --prefix frontend run build`` still starts, still serves ``/v4``,
``/docs`` and the rest of the API, and answers page requests with a 503 that
says how to fix it.
"""

import json
import logging
import os
import re
from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import FileResponse, HTMLResponse, Response

from refgenie.server.const import API_PATH_PREFIXES
from refgenie.server.errors import ErrorCode, error_response

logger = logging.getLogger(__name__)

#: Directory name of the built bundle inside the package.
WEBUI_DIRNAME = "webui"
#: Vite's ``build.assetsDir``: content-hashed files, safe to cache forever.
ASSETS_DIRNAME = "_app"
#: Written by the frontend build (``frontend/scripts/stamp-build.mjs``).
BUILD_INFO_FILENAME = "build-info.json"

_IMMUTABLE = "public, max-age=31536000, immutable"
_NO_CACHE = "no-cache"

#: ``frontend/index.html`` ships a literal ``<base href="/" />``; this matches it
#: however it is spelled (quoting and self-closing slash vary with formatting).
_BASE_HREF_RE = re.compile(r"""<base\s+href=["'][^"']*["']\s*/?>""", re.IGNORECASE)

#: Vite's ``base`` is ``/``, so the bundle's own script/style tags reference
#: ``/_app/...`` absolutely -- and ``<base href>`` does not rewrite absolute
#: URLs. A sub-path deployment therefore needs these rewritten too.
_ABS_ASSET_RE = re.compile(rf"""(src|href)=(["'])/{ASSETS_DIRNAME}/""")

MISSING_BUNDLE_HTML = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>refgenie &mdash; web UI not built</title>
<style>
  body { font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
         max-width: 42rem; margin: 4rem auto; padding: 0 1.5rem;
         line-height: 1.6; color: #1c1c1c; background: #fff; }
  h1 { font-size: 1.6rem; margin-bottom: 0.25rem; }
  code { background: #f2f2f2; padding: 0.15rem 0.35rem; border-radius: 3px; }
  pre { background: #f2f2f2; padding: 0.75rem 1rem; border-radius: 4px;
        overflow-x: auto; }
  .muted { color: #666; }
  @media (prefers-color-scheme: dark) {
    body { color: #e8e8e8; background: #1a1a1a; }
    code, pre { background: #2a2a2a; }
    .muted { color: #a0a0a0; }
    a { color: #7fb2ff; }
  }
</style>
</head>
<body>
<h1>refgenie</h1>
<p class="muted">The web UI bundle is not built, so there is no page to show.
The API is running normally.</p>
<p>Running from a source checkout? Build the bundle:</p>
<pre>npm --prefix frontend run build</pre>
<p>Running from an installed wheel? The bundle ships inside it, so a missing
one means the wheel build is broken &mdash; please report a bug.</p>
<p>In the meantime:
<a href="docs">API documentation</a> &middot;
<a href="v4/genomes">/v4/genomes</a></p>
</body>
</html>
"""


def resolve_web_dist(explicit: "Path | None" = None) -> "Path | None":
    """Locate the built web UI directory. First hit wins; never raises.

    1. ``explicit`` (tests pass a tmp dir),
    2. ``$REFGENIE_WEB_DIST``,
    3. the packaged/dev layout, ``refgenie/server/webui/``.

    There is no ``frontend/dist`` fallback: vite builds into the package.

    Returns ``None`` when no ``index.html`` is found.
    """
    candidates = []
    if explicit is not None:
        candidates.append(Path(explicit))
    env_dist = os.environ.get("REFGENIE_WEB_DIST")
    if env_dist:
        candidates.append(Path(env_dist))
    candidates.append(Path(__file__).parent / WEBUI_DIRNAME)

    for candidate in candidates:
        try:
            if (candidate / "index.html").is_file():
                return candidate
        except OSError:  # unreadable path; treat as absent
            continue
    return None


def read_build_info(dist: "Path | None") -> dict:
    """The frontend build stamp (``build-info.json``), or ``{}``.

    Reported under ``/service-info`` as ``refgenie.web_ui`` so a deployed UI can
    be identified without guessing from asset hashes. Never raises.
    """
    if dist is None:
        return {}
    try:
        return json.loads((dist / BUILD_INFO_FILENAME).read_text())
    except (OSError, ValueError):
        return {}


def _render_index(dist: Path, root_path: str) -> str:
    """``index.html`` with its ``<base href>`` pointed at ``root_path``.

    Rendered once, at app construction, and cached: a sub-path deployment is a
    one-token difference, not a reason to rebuild the frontend.
    """
    html = (dist / "index.html").read_text()
    if root_path:
        prefix = root_path.rstrip("/")
        html = _BASE_HREF_RE.sub(f'<base href="{prefix}/" />', html, count=1)
        # The bundle's own <script>/<link> tags point at /_app/... absolutely,
        # which <base href> does not touch. Rewrite them in the same pass, or a
        # sub-path deployment loads its HTML and none of its code.
        html = _ABS_ASSET_RE.sub(rf"\1=\g<2>{prefix}/{ASSETS_DIRNAME}/", html)
    return html


def _is_api_path(path: str) -> bool:
    """True if ``path`` belongs to the API rather than to the SPA."""
    return any(path.startswith(prefix) for prefix in API_PATH_PREFIXES)


def _json_404(full_path: str):
    """A JSON 404 in the standard ``{"ok": false, "error": {...}}`` envelope.

    The catch-all builds its responses directly, so the local-mode exception
    handlers never see these -- the envelope has to be emitted here or an
    unmatched API path is the one non-2xx that breaks the contract.
    """
    return error_response(404, ErrorCode.NOT_FOUND, f"Not Found: /{full_path}")


def mount_spa(app: FastAPI, dist: "Path | None") -> None:
    """Register the SPA catch-all on ``app``. Must be registered LAST.

    The catch-all matches every path, so any route added after it is dead. See
    the ordering rule in ``refgenie.server.main.create_app``.
    """
    root_path = app.root_path or ""

    index_html: str | None = None
    if dist is not None:
        try:
            index_html = _render_index(dist, root_path)
        except OSError:
            logger.warning("Web UI index.html at %s is unreadable; serving the 503 page", dist)
            dist = None
    if dist is None:
        logger.warning(
            "No web UI bundle found; / will answer 503. "
            "On a source checkout, build it with `npm --prefix frontend run build`; "
            "on an installed wheel the bundle ships inside the package, so a "
            "missing one is a broken wheel build -- please report a bug."
        )

    dist_root = os.path.realpath(dist) if dist is not None else None

    def _index_response() -> Response:
        return HTMLResponse(index_html, headers={"Cache-Control": _NO_CACHE})

    @app.api_route("/{full_path:path}", methods=["GET", "HEAD"], include_in_schema=False)
    async def spa_catch_all(full_path: str) -> Response:
        # A real file in the bundle wins first -- this is what serves the
        # hashed assets under /_app/, which are content-addressed and so are
        # cacheable forever.
        if dist_root is not None and full_path:
            candidate = os.path.realpath(os.path.join(dist_root, full_path))
            # Containment check: `..` segments and symlinks out of the bundle
            # must not become a file read.
            if not (candidate == dist_root or candidate.startswith(dist_root + os.sep)):
                return _json_404(full_path)
            if os.path.isfile(candidate):
                immutable = full_path.startswith(f"{ASSETS_DIRNAME}/")
                return FileResponse(
                    candidate,
                    headers={"Cache-Control": _IMMUTABLE if immutable else _NO_CACHE},
                )

        # No such file. A path under an API prefix is then a 404, in JSON:
        # handing back index.html with a 200 would make a client's typo look
        # like a successful request that returned unparseable data.
        if _is_api_path("/" + full_path):
            return _json_404(full_path)

        if dist_root is None:
            return HTMLResponse(MISSING_BUNDLE_HTML, status_code=503)

        # Anything else is a client-side route: serve the shell.
        return _index_response()
