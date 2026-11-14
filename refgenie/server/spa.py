"""Serve the built web UI (the React SPA) from the refgenie app.

Vite builds ``frontend/`` straight into ``refgenie/server/webui/``. Behind a
proxy prefix, ``create_app(root_path=...)`` rewrites ``<base href>`` once at
startup, so no rebuild is needed. A missing bundle is not fatal: the API
still runs and pages answer 503 with build instructions.
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

# Python before 3.13 has no built-in type for .md, and macOS has no
# /etc/mime.types to supply one, so SKILL.md would go out as octet-stream.
_MEDIA_TYPES = {".md": "text/markdown"}

#: Directory name of the built bundle inside the package.
WEBUI_DIRNAME = "webui"
#: Vite's ``build.assetsDir``: content-hashed files, safe to cache forever.
ASSETS_DIRNAME = "_app"
#: Written by the frontend build (``frontend/scripts/stamp-build.mjs``).
BUILD_INFO_FILENAME = "build-info.json"

_IMMUTABLE = "public, max-age=31536000, immutable"
_NO_CACHE = "no-cache"

#: Matches the ``<base href="/" />`` in ``frontend/index.html``, however it is spelled.
_BASE_HREF_RE = re.compile(r"""<base\s+href=["'][^"']*["']\s*/?>""", re.IGNORECASE)

#: Absolute ``/_app/...`` asset URLs, which ``<base href>`` does not rewrite.
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


#: What a public server (``refgenie serve``) shows at ``/``. The web UI is
#: hosted separately, at refgenie.org.
SERVER_LANDING_HTML = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>refgenie API</title>
<style>
  body { font-family: system-ui, sans-serif; max-width: 40rem; margin: 4rem auto;
         padding: 0 1.5rem; line-height: 1.6; }
</style>
</head>
<body>
<h1>refgenie API</h1>
<p>This is a refgenie server. To browse genomes, go to
<a href="https://refgenie.org">refgenie.org</a>.</p>
<ul>
  <li><a href="https://docs.refgenie.org">Documentation</a></li>
  <li><a href="docs">API reference</a> (<a href="openapi.json">OpenAPI</a>)</li>
  <li>For AI agents: <a href="llms.txt">llms.txt</a>, <a href="SKILL.md">SKILL.md</a></li>
</ul>
</body>
</html>
"""


def resolve_web_dist(explicit: "Path | None" = None) -> "Path | None":
    """Locate the built web UI directory, or None. Never raises.

    Tries ``explicit``, then ``$REFGENIE_WEB_DIST``, then ``refgenie/server/webui/``.
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
    """The frontend build stamp (``build-info.json``), or ``{}``. Never raises."""
    if dist is None:
        return {}
    try:
        return json.loads((dist / BUILD_INFO_FILENAME).read_text())
    except (OSError, ValueError):
        return {}


def _render_index(dist: Path, root_path: str) -> str:
    """``index.html`` with its ``<base href>`` and asset URLs pointed at ``root_path``."""
    html = (dist / "index.html").read_text()
    if root_path:
        prefix = root_path.rstrip("/")
        html = _BASE_HREF_RE.sub(f'<base href="{prefix}/" />', html, count=1)
        html = _ABS_ASSET_RE.sub(rf"\1=\g<2>{prefix}/{ASSETS_DIRNAME}/", html)
    return html


def _is_api_path(path: str) -> bool:
    """True if ``path`` belongs to the API rather than to the SPA."""
    return any(path.startswith(prefix) for prefix in API_PATH_PREFIXES)


def _json_404(full_path: str):
    """A JSON 404 in the standard error envelope (exception handlers do not see these)."""
    return error_response(404, ErrorCode.NOT_FOUND, f"Not Found: /{full_path}")


def mount_spa(app: FastAPI, dist: "Path | None", shell: bool = True) -> None:
    """Register the SPA catch-all on ``app``. Must be last: it matches every path.

    With ``shell=False`` (a public server) ``/`` is a short landing page, other
    page paths 404, and only static bundle files such as ``SKILL.md`` are served.
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
        # A real bundle file wins; hashed /_app/ assets are cached forever.
        if dist_root is not None and full_path:
            candidate = os.path.realpath(os.path.join(dist_root, full_path))
            # Refuse `..` and symlinks that leave the bundle.
            if not (candidate == dist_root or candidate.startswith(dist_root + os.sep)):
                return _json_404(full_path)
            if os.path.isfile(candidate):
                immutable = full_path.startswith(f"{ASSETS_DIRNAME}/")
                return FileResponse(
                    candidate,
                    media_type=_MEDIA_TYPES.get(os.path.splitext(candidate)[1]),
                    headers={"Cache-Control": _IMMUTABLE if immutable else _NO_CACHE},
                )

        if not shell:
            return HTMLResponse(SERVER_LANDING_HTML) if not full_path else _json_404(full_path)

        # An unmatched API path is a JSON 404, not the HTML shell.
        if _is_api_path("/" + full_path):
            return _json_404(full_path)

        if dist_root is None:
            return HTMLResponse(MISSING_BUNDLE_HTML, status_code=503)

        # Anything else is a client-side route: serve the shell.
        return _index_response()
