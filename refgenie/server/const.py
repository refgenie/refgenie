"""Web-layer constants for the refgenie app factory (both modes)."""

import os
from typing import Literal

from importlib.metadata import version

ALL_VERSIONS = {"version": version("refgenie")}

# --- App modes -------------------------------------------------------------
#
# One factory, two modes. `server` is the public REST API (DRS, data channels,
# seqcol, MCP, download counting); `local` is the single-user management app
# behind `refgenie dash`. These name a *web-layer* concern and are unrelated to
# refgenie/core/mode.py's LocalMode/ServerMode strategy objects -- see create_app.
APP_MODE_SERVER = "server"
APP_MODE_LOCAL = "local"
AppMode = Literal["server", "local"]

# --- Local mode -------------------------------------------------------------
#
# There are no local-security constants here on purpose:
# refgenie/server/local/security.py is the single owner of local CORS, the host
# guard, the action header and the bridge policy, and it carries its defaults
# on LocalSecuritySettings -- change them there.

# --- Web UI (SPA) -----------------------------------------------------------
#
# The client-side routes the SPA owns. The root namespace is reserved for them:
# the JSON API is mounted at /v4 and /v1 only. The non-collision guard in
# tests/test_web.py keeps it that way.
#
# Localhost-bridge deep-link contract (treat any change as a coordinated
# cross-repo change with the public SPA -- it hands off to these URLs):
#
#   http://localhost:{port}/genomes/{digest}
#   http://localhost:{port}/genomes/{digest}/{assetDigest}
#   http://localhost:{port}/pull?server={remoteApiBase}&genome={digest}&asset_group={name}
#
# The /pull entry lands on a prefilled, UNSUBMITTED pull confirmation. It must
# never auto-execute from a URL: that would recreate the CSRF hole the
# X-Refgenie-Action header closes.
SPA_CLIENT_ROUTES: tuple[str, ...] = (
    "genomes",
    "assets",
    "asset-groups",
    "asset-classes",
    "recipes",
    "aliases",
    "remote",
    "manage",
    "jobs",
    "build",
    "pull",
    "about",
)

#: Path prefixes owned by the API. The SPA catch-all answers an unmatched path
#: under one of these with a JSON 404 instead of handing back index.html -- an
#: API typo must never look like a successful page load.
API_PATH_PREFIXES: tuple[str, ...] = (
    "/v4/",
    "/v1/",
    "/seqcol",
    "/ga4gh/",
    "/data_channel",
    "/mcp",
    "/openapi.json",
    "/docs",
    "/redoc",
    "/service-info",
    "/ping",
    "/_app/",
)

# Shared service-info identity, used by the root /service-info document and by
# the seqcol sub-service so the two agree.
SERVICE_ORGANIZATION = {"name": "Refgenie", "url": "https://refgenie.databio.org"}
SERVICE_CONTACT_URL = "https://github.com/refgenie/refgenie/issues"
SERVICE_DOCUMENTATION_URL = "https://refgenie.databio.org"
SERVICE_GITHUB_URL = "https://github.com/refgenie/refgenie"

# Where the GA4GH sequence collections service is mounted. Kept off the root so
# refget's /collection, /list/*, /aliases/{kind}/{digest} and /service-info do
# not collide with refgenie's own resource API.
SEQCOL_MOUNT_PATH = "/seqcol"
SEQCOL_SERVICE_ID = "org.refgenie.seqcol"
SEQCOL_SERVICE_NAME = "Refgenie Sequence Collections"
DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS = int(
    os.environ.get("DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS", 60)
)

# Publish-catalog artifact to import at startup and daily (see
# refgenie/catalog_transfer.py). Unset => the server serves whatever its SQL
# catalog already holds.
REFGENIE_CATALOG_URL = os.environ.get("REFGENIE_CATALOG_URL", None)
CATALOG_IMPORT_INTERVAL_SECONDS = int(os.environ.get("CATALOG_IMPORT_INTERVAL_SECONDS", 24 * 3600))
