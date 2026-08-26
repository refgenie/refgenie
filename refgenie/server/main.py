try:
    # the following imports are only included in extras 'dash'
    import uvicorn
    from fastapi import FastAPI
    from fastapi.middleware.cors import CORSMiddleware
except ImportError as e:
    raise ImportError(
        "The 'dash' extras are not installed. Please install refgenie with the 'dash' extras to use the web interface."
    ) from e
import logging
from contextlib import asynccontextmanager
from pathlib import Path

from refgenie.server.dependencies import get_refgenie as _get_refgenie
from refgenie.server.routers.shared import router as refgenie_core_router
from refgenie.core import Refgenie
from refget.router import create_refget_router, setup_backend

from refgenie.const import API_VERSION
from refgenie.server.const import (
    ALL_VERSIONS,
    APP_MODE_LOCAL,
    APP_MODE_SERVER,
    AppMode,
    CATALOG_IMPORT_INTERVAL_SECONDS,
    DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS,
    REFGENIE_CATALOG_URL,
    SEQCOL_MOUNT_PATH,
    SEQCOL_SERVICE_ID,
    SEQCOL_SERVICE_NAME,
    SERVICE_CONTACT_URL,
    SERVICE_DOCUMENTATION_URL,
    SERVICE_GITHUB_URL,
    SERVICE_ORGANIZATION,
)
from refgenie.server.spa import mount_spa, read_build_info, resolve_web_dist
from refgenie.server.stats import (
    EndpointCollector,
    EndpointHitCounterMiddleware,
    handle_endpoint_collector_hits,
)
from refgenie.db.tables import StoreType

logger = logging.getLogger(__name__)


SERVER_TAGS_METADATA = [
    {
        "name": "Version 4",
        "description": "Version 4 of the Refgenieserver API, which includes support for pangenomes and sequences.",
    },
    {
        "name": "Default",
        "description": "Default endpoints for the Refgenieserver API, which includes basic genome asset management.",
    },
    {
        "name": "Data Channel",
        "description": "Endpoints for the Data Channel, which provides access to genome assets in a data channel format.",
    },
    {
        "name": "GA4GH DRS",
        "description": "Endpoints for the GA4GH DRS (Data Repository Service), which provides access to genomic data in a standardized format. More details: [GA4GH DRS standard](https://ga4gh.github.io/data-repository-service-schemas/)",
    },
]

LOCAL_TAGS_METADATA = [
    {
        "name": "Version 4",
        "description": "The refgenie JSON API: genomes, assets, recipes and asset classes.",
    },
    {
        "name": "Actions",
        "description": "State-changing operations on this local refgenie instance (pull, build, delete, aliases, subscriptions).",
    },
    {
        "name": "Remote",
        "description": "Browsing the catalogs of the servers this instance subscribes to.",
    },
]

#: Every capability key the web UI knows about. The UI gates each affordance on
#: a flag, never on the mode, so "server mode" is just "local mode minus the
#: command surface". A key missing from a /service-info document reads as false;
#: adding a key here is a simultaneous change in the frontend's Capabilities type.
CAPABILITY_KEYS = (
    "pull",
    "build",
    "delete",
    "aliases_write",
    "subscriptions",
    "recipes_write",
    "asset_classes_write",
    "remote_browse",
    "genome_init",
    "jobs",
    "jobs_cancel",
    "downloads",
    "archives",
    "seqcol",
    "drs",
)


def _capabilities(mode: AppMode) -> dict[str, bool]:
    """The capability flags for ``mode``.

    Server mode serves data and nothing else: no writes, no remote browsing, no
    jobs. Local mode is the mirror image -- it manages one user's assets, but
    has no archives to download and no seqcol or DRS service.

    ``recipes_write`` / ``asset_classes_write`` are false in both modes: adding
    a recipe or asset class takes a server-side path or URL, which over HTTP is
    an arbitrary local-file read and an SSRF. Those stay CLI-only.
    """
    is_local = mode == APP_MODE_LOCAL
    return {
        "pull": is_local,
        "build": is_local,
        "delete": is_local,
        "aliases_write": is_local,
        "subscriptions": is_local,
        "recipes_write": False,
        "asset_classes_write": False,
        "remote_browse": is_local,
        "genome_init": is_local,
        "jobs": is_local,
        "jobs_cancel": is_local,
        "downloads": not is_local,
        "archives": not is_local,
        "seqcol": not is_local,
        "drs": not is_local,
    }


def _configure_cors(app: FastAPI, mode: AppMode) -> None:
    """Install CORS for ``server`` mode. Local mode is handled elsewhere.

    The public API is read-only and unauthenticated, so it is open to every
    origin -- with ``allow_credentials=False``, because a wildcard origin
    combined with credentials is an invalid combination that browsers reject
    outright (the previous configuration did exactly that).

    Local-mode CORS is *not* set here. It belongs to
    ``refgenie.server.local.security.install_local_security``, which is the
    single owner of the local origin allowlist, the host guard and the
    action-header dependency -- three things that only make sense together.
    """
    if mode != APP_MODE_SERVER:
        return
    app.add_middleware(
        CORSMiddleware,
        allow_origins=["*"],
        allow_credentials=False,
        allow_methods=["GET", "HEAD", "OPTIONS"],
        allow_headers=["*"],
    )


def create_app(
    mode: AppMode = APP_MODE_SERVER,
    root_path: str = "",
    refgenie_instance: "Refgenie | None" = None,
    web_dist: "Path | None" = None,
) -> FastAPI:
    """Create and configure a Refgenie FastAPI application.

    One factory, two modes:

    * ``server`` -- the public REST API: the v4 JSON API, GA4GH DRS, data
      channels, a mounted seqcol service, ``/mcp``, download counting and the
      background scheduler.
    * ``local`` -- ``refgenie dash``: the same v4 JSON API over the user's own
      refgenie, plus the local-only command surface under ``/v1``.

    Both modes serve the same web UI bundle; the UI discovers which mode it is
    talking to from ``/service-info``.

    Note that this ``mode`` is a *web-layer* concept. It is unrelated to
    ``refgenie/core/mode.py``'s ``LocalMode`` / ``ServerMode`` strategy objects on
    the ``Refgenie`` instance, and must never be derived from or validated
    against ``rg.mode``: tests legitimately build a server-mode app over a
    local-mode ``Refgenie``.

    Args:
        mode: ``"server"`` (default) or ``"local"``.
        root_path: URL prefix for sub-path mounting (e.g. "/refgenie").
            Passed to FastAPI(root_path=...) and used to rewrite the SPA's
            ``<base href>`` so a sub-path deployment needs no frontend rebuild.
        refgenie_instance: Optional pre-built Refgenie instance. If provided,
            it is used instead of creating a new one (useful for testing).
        web_dist: Optional path to a built web UI bundle, overriding the
            packaged ``refgenie/server/webui/``. Tests pass a tmp dir here.

    Returns:
        A fully configured FastAPI app ready to serve or mount as a sub-app.
    """
    if mode not in (APP_MODE_SERVER, APP_MODE_LOCAL):
        raise ValueError(f"Unknown app mode: {mode!r}")
    is_local = mode == APP_MODE_LOCAL

    if refgenie_instance is not None:
        rg = refgenie_instance
    elif is_local:
        # The user's real config singleton -- this app manages the caller's own
        # refgenie installation.
        rg = _get_refgenie()
    else:
        rg = Refgenie(server_mode=True)
        # Register declared stores (REFGENIE_STORES) before enumerating them.
        from refgenie.server.store_bootstrap import bootstrap_stores_from_env

        bootstrap_stores_from_env(rg)

    # The set of stores this server federates over comes from the `store`
    # registry table -- the single source of truth. Empty (or local mode) means
    # DB-backed seqcol with refget_store.enabled=false, exactly as a store-less
    # server behaved before. Read here at construction because the seqcol mounts
    # must be registered before the SPA catch-all.
    store_rows = []
    if not is_local:
        try:
            store_rows = rg.store.enabled_stores()
        except Exception:  # noqa: BLE001 - a missing/empty registry is DB mode
            logger.exception("Could not read the store registry; serving DB-mode seqcol")
            store_rows = []

    # Each federated store is a mounted seqcol sub-application (see the /seqcol
    # wiring below). Mounted sub-apps do not receive lifespan events, so their
    # backends are bound from this app's lifespan; the list is filled in during
    # app construction and closed over there. Each entry is
    # ``(store_row, seqcol_app, mount_path)``.
    seqcol_mounts: list = []

    def _import_catalog():
        """Refresh the SQL catalog from the published artifact; never fatal.

        A failed import leaves the previous catalog serving (or, on first
        boot, an empty one) -- the store-backed seqcol endpoints work either
        way, so a bad artifact must not take the server down with it.
        """
        from refgenie.catalog_transfer import import_publish_catalog

        try:
            import_publish_catalog(rg.database_engine, REFGENIE_CATALOG_URL)
        except Exception:
            logger.exception(f"Publish-catalog import from {REFGENIE_CATALOG_URL} failed")

    @asynccontextmanager
    async def local_lifespan(app: FastAPI):
        """Local mode: make sure the user's database exists, and nothing else.

        No scheduler, no catalog import, no seqcol backend, no MCP session
        manager -- none of those exist in this mode.

        The in-process JobManager is created here rather than at app
        construction so its worker threads and its handler on the `refgenie`
        logger exist only while the app is actually serving -- a test that
        builds an app without entering its lifespan leaves no threads behind.
        """
        from refgenie.server.jobs import JobManager

        rg._create_db_and_tables()
        app.state.job_manager = JobManager(rg)
        try:
            yield
        finally:
            # wait=False: at process exit nobody is watching, and blocking
            # uvicorn's shutdown on a two-hour build helps no one. Abandoned
            # jobs are named in a warning.
            app.state.job_manager.shutdown(wait=False)
            app.state.job_manager = None

    @asynccontextmanager
    async def server_lifespan(app: FastAPI):
        # apscheduler is a `server`-extra dependency and is imported here, not
        # at module scope: `refgenie dash` imports this module on an install
        # that has only the `dash` extra.
        from apscheduler.schedulers.background import BackgroundScheduler

        print("Creating database and tables")
        rg._create_db_and_tables()
        if REFGENIE_CATALOG_URL:
            _import_catalog()
        if seqcol_mounts:
            # Give each seqcol mount its own readonly snapshot of its store.
            #
            # Each must be a *second* store instance, not one from rg's router:
            # gtars' into_readonly() consumes its receiver, so converting a
            # shared instance would leave the v4 router and managers holding an
            # empty store, and rg's stores also get mutated elsewhere
            # (_import_remote_collection, _fetch_remote_sequence), which a
            # readonly snapshot cannot support.
            import tempfile

            from refget.seqcolapi import prepare_store

            for store_row, sc_app, _mount_path in seqcol_mounts:
                cache_dir = tempfile.mkdtemp(prefix="refgenie_seqcol_cache_")
                setup_backend(
                    sc_app,
                    store=prepare_store(
                        store_row.url,
                        remote=store_row.type == StoreType.remote,
                        cache_dir=cache_dir,
                    ),
                )
        else:
            setup_backend(app, engine=rg.database_engine)
        app.state.endpoint_hits_collector = EndpointCollector()
        scheduler = BackgroundScheduler()
        scheduler.add_job(
            handle_endpoint_collector_hits,
            "interval",
            seconds=DOWNLOAD_COUNT_DUMP_JOB_INTERVAL_SECONDS,
            args=[app, rg],
        )
        if REFGENIE_CATALOG_URL:
            # Idempotent upsert, so the per-worker scheduler is harmless.
            scheduler.add_job(
                _import_catalog,
                "interval",
                seconds=CATALOG_IMPORT_INTERVAL_SECONDS,
            )
        scheduler.start()
        # Start the MCP StreamableHTTP session manager's task group. The /mcp
        # streamable app is mounted as a sub-application (below), and mounted
        # sub-apps do not receive lifespan events -- same reason the seqcol
        # backend is bound from this lifespan. Use the per-app instance captured
        # at mount time: the module-level MCP server is a singleton, and every
        # streamable_http_app() call overwrites its .session_manager, so reading
        # it here would race across app instances; .run() is also once-per-manager.
        # Without this, every /mcp request fails with "Task group is not initialized".
        async with app.state.mcp_session_manager.run():
            yield
        scheduler.shutdown()

    app = FastAPI(
        title="Refgenie" if is_local else "Refgenieserver REST API",
        description=(
            "the local management API for one refgenie installation"
            if is_local
            else "a web interface and RESTful API for reference genome assets"
        ),
        version=ALL_VERSIONS["version"],
        lifespan=local_lifespan if is_local else server_lifespan,
        openapi_tags=LOCAL_TAGS_METADATA if is_local else SERVER_TAGS_METADATA,
        swagger_ui_parameters={"docExpansion": "none", "defaultModelsExpandDepth": 0},
        root_path=root_path,
    )
    app.state.mode = mode

    # --- Registration order matters ----------------------------------------
    #
    #   middleware -> dependency overrides -> API routers -> sub-app mounts ->
    #   SPA fallback LAST
    #
    # The SPA catch-all matches every path, so anything registered after it is
    # dead. Keep new routes above mount_spa().

    _configure_cors(app, mode)

    # Override the shared dependency so the routers use this app's Refgenie
    # instance. get_db_session takes get_refgenie through Depends, so this also
    # guarantees the handler's session is on the same database.
    app.dependency_overrides[_get_refgenie] = lambda: rg

    if not is_local:
        app.add_middleware(EndpointHitCounterMiddleware)

    if is_local:
        # Local CORS, the Host-header guard and the X-Refgenie-Action
        # dependency travel together; local/security.py is their single owner.
        # The error handlers make every non-2xx from this app use the
        # {"ok": false, "error": {...}} envelope the SPA codes against.
        from refgenie.server.local import (
            actions_router,
            install_error_handlers,
            install_local_security,
        )

        install_local_security(app)
        install_error_handlers(app)

        # The state-changing command surface. Local mode only -- several of its
        # request models accept server-local filesystem paths, which is safe
        # only on a loopback-bound, header-guarded, single-user app.
        app.include_router(actions_router, prefix="/v1")

        from refgenie.server.routers.remote import router as remote_router

        app.include_router(remote_router, prefix="/v1")

        # Background jobs: pull and build run on dedicated worker threads and
        # report over /v1/jobs/events. The manager itself is created and torn
        # down in local_lifespan, and the router reads it off app.state, so
        # importing this router costs nothing at construction time.
        from refgenie.server.jobs import jobs_router

        app.include_router(jobs_router, prefix="/v1")

    app.include_router(refgenie_core_router, prefix="/v4", tags=["Version 4"])

    if not is_local:
        # Import routers after the app is defined to avoid circular imports
        from refgenie.server.routers.version4 import router as version4

        app.include_router(version4, prefix="/v4", tags=["Version 4"])

        from refgenie.server.routers.ga4gh_drs import router as ga4gh_drs

        app.include_router(ga4gh_drs, prefix="/v4/ga4gh/drs", tags=["GA4GH DRS"])
        app.include_router(ga4gh_drs, prefix="/ga4gh/drs", tags=["GA4GH DRS"])

        from refgenie.server.routers.data_channel import router as data_channel

        app.include_router(data_channel, prefix="/v4/data_channel", tags=["Data Channel"])
        app.include_router(data_channel, prefix="/data_channel", tags=["Data Channel"])

    # --- Sequence collections (seqcol) --- (server mode only)
    #
    # Store mode mounts refget's shared store-backed seqcol app, which brings
    # the readonly-store contract, freshness reloading, the seqcol JSON schema
    # and a complete GA4GH service-info. It is *mounted*, not wired into this
    # app: refget's setup_backend binds to app.state.backend, so a router
    # included twice at two prefixes would still serve one store. A mounted
    # sub-app owns its own backend, which is what would make serving several
    # stores possible later. Only one store is served today.
    #
    # DB mode has no store, so it keeps the plain router plus a minimal
    # service-info that reports refget_store.enabled: false.
    if not is_local:
        if store_rows:
            from refget.seqcolapi import create_seqcol_app

            # Highest-priority store (index 0) is the default mount at
            # SEQCOL_MOUNT_PATH; every other store gets a name-qualified mount.
            for idx, store_row in enumerate(store_rows):
                mount_path = (
                    SEQCOL_MOUNT_PATH if idx == 0 else f"{SEQCOL_MOUNT_PATH}/{store_row.name}"
                )
                service_id = (
                    SEQCOL_SERVICE_ID if idx == 0 else f"{SEQCOL_SERVICE_ID}.{store_row.name}"
                )
                service_name = (
                    SEQCOL_SERVICE_NAME
                    if idx == 0
                    else f"{SEQCOL_SERVICE_NAME} ({store_row.name})"
                )
                sc_app = create_seqcol_app(
                    store_path=store_row.url,
                    remote=store_row.type == StoreType.remote,
                    store_url=store_row.url,
                    service_info_id=service_id,
                    service_info_name=service_name,
                    organization=SERVICE_ORGANIZATION,
                    contact_url=SERVICE_CONTACT_URL,
                    documentation_url=SERVICE_DOCUMENTATION_URL,
                    freshness=True,
                    cors=False,  # the host app's CORS middleware does not reach a mount
                    # Routes now, store on startup -- see the lifespan above.
                    defer_backend=True,
                )
                app.mount(mount_path, sc_app)
                seqcol_mounts.append((store_row, sc_app, mount_path))
        else:
            app.include_router(
                create_refget_router(
                    sequences=False,
                    pangenomes=False,
                    refget_store_url=None,
                    mount_prefix=SEQCOL_MOUNT_PATH,
                ),
                prefix=SEQCOL_MOUNT_PATH,
            )

            # This handler only exists in DB mode, so the only true answer is
            # "no store". Store mode's service-info comes from the mounted
            # sub-app above, which reports the store URL itself.
            @app.get(f"{SEQCOL_MOUNT_PATH}/service-info", tags=["Sequence Collections"])
            async def seqcol_service_info():
                return {
                    "id": SEQCOL_SERVICE_ID,
                    "name": SEQCOL_SERVICE_NAME,
                    "type": {"group": "org.ga4gh", "artifact": "refget-seqcol", "version": "1.0.0"},
                    "seqcol": {"refget_store": {"enabled": False}},
                }

    # --- The web UI bundle --------------------------------------------------
    web_dist = resolve_web_dist(web_dist)
    build_info = read_build_info(web_dist)
    capabilities = _capabilities(mode)
    # /ping (the localhost-bridge presence probe, both modes) reads these off
    # app.state rather than importing this module -- no circular import.
    app.state.capabilities = capabilities

    from refgenie.server.routers.ping import router as ping_router

    app.include_router(ping_router)

    @app.get("/service-info", tags=["Default"])
    async def service_info():
        """GA4GH-style discovery document, and the web UI's bootstrap.

        Two audiences, cleanly separated. The GA4GH fields stay at the top
        level; everything refgenie-specific nests under ``refgenie`` so a
        standards client sees the shape it expects.

        For a refgenie client this is the one well-known fetch that turns a
        server URL into a RefgetStore URL; everything after that goes straight
        to the store. ``seqcol.refget_store.url`` is the bootstrap field.

        ``seqcol.services`` lists the seqcol mounts. It is a list today with at
        most one entry, so that serving more than one store later adds entries
        rather than changing the document's shape. The flat ``seqcol.url`` and
        ``seqcol.refget_store`` keys describe the first (currently only) mount.

        For the web UI this is how one bundle serves both modes: it reads
        ``refgenie.capabilities`` and gates every affordance on a flag. Local
        mode has no seqcol mount, so the ``seqcol`` block is absent rather than
        describing a service that is not there.
        """
        document = {
            "id": "org.refgenie.api",
            "name": "Refgenie",
            "type": {
                "group": "org.refgenie",
                "artifact": "refgenie",
                "version": API_VERSION,
            },
            "description": (
                "A web interface and RESTful API for reference genome assets, "
                "with a GA4GH sequence collections service"
            ),
            "organization": SERVICE_ORGANIZATION,
            "contactUrl": SERVICE_CONTACT_URL,
            "documentationUrl": SERVICE_DOCUMENTATION_URL,
            "environment": "production",
            "version": ALL_VERSIONS["version"],
            "refgenie": {
                "mode": mode,
                "api_base": "/v4",
                "root_path": root_path,
                "refgenie_version": ALL_VERSIONS["version"],
                "service_name": ("refgenie local dashboard" if is_local else "refgenie server"),
                "web_ui": build_info,
                "capabilities": capabilities,
                "links": {
                    "docs": SERVICE_DOCUMENTATION_URL,
                    "github": SERVICE_GITHUB_URL,
                    "openapi": "/openapi.json",
                },
            },
        }
        if not is_local:
            if seqcol_mounts:
                mounts = []
                for idx, (store_row, _sc_app, mount_path) in enumerate(seqcol_mounts):
                    mounts.append(
                        {
                            "id": (
                                SEQCOL_SERVICE_ID
                                if idx == 0
                                else f"{SEQCOL_SERVICE_ID}.{store_row.name}"
                            ),
                            "name": (
                                SEQCOL_SERVICE_NAME
                                if idx == 0
                                else f"{SEQCOL_SERVICE_NAME} ({store_row.name})"
                            ),
                            "store": store_row.name,
                            "url": mount_path,
                            "service_info": f"{mount_path}/service-info",
                            "refget_store": {"enabled": True, "url": store_row.url},
                        }
                    )
                # The flat keys describe the highest-priority (default) mount.
                default = mounts[0]
                document["seqcol"] = {**default, "services": mounts}
            else:
                mount = {
                    "id": SEQCOL_SERVICE_ID,
                    "name": SEQCOL_SERVICE_NAME,
                    "url": SEQCOL_MOUNT_PATH,
                    "service_info": f"{SEQCOL_MOUNT_PATH}/service-info",
                    "refget_store": {"enabled": False},
                }
                document["seqcol"] = {**mount, "services": [mount]}
        return document

    if not is_local:
        # MCP endpoint (Streamable HTTP)
        from mcp.server.transport_security import TransportSecuritySettings
        from refgenie.mcp.tools import mcp as mcp_server, set_refgenie as set_mcp_refgenie

        set_mcp_refgenie(rg)
        # mcp 2.x defaults streamable_http_app(host="127.0.0.1"), which enables
        # DNS-rebinding Host-header validation and 421s any real host. This API is
        # public, proxied, and already allow_origins=["*"], so that browser-oriented
        # guard does not apply -- disable it or the mount rejects api.refgenie.org.
        #
        # Local mode is the opposite case: it is browser-facing and loopback-only,
        # so it gets the same class of guard *enabled*, as
        # HostHeaderGuardMiddleware (refgenie/server/local/security.py).
        mcp_app = mcp_server.streamable_http_app(
            transport_security=TransportSecuritySettings(enable_dns_rebinding_protection=False)
        )
        # Bind THIS app's session manager (the lifespan starts its task group).
        # streamable_http_app() just created a fresh manager and stored it on the
        # module-singleton server; capture it now, before another create_app() call
        # can overwrite the singleton's reference.
        app.state.mcp_session_manager = mcp_server.session_manager
        app.mount("/mcp", mcp_app)

    # LAST: the catch-all that serves the SPA and deep-links its client routes.
    mount_spa(app, web_dist)

    return app


def create_local_app() -> FastAPI:
    """Factory entry point: uvicorn refgenie.server.main:create_local_app --factory"""
    return create_app(mode=APP_MODE_LOCAL)


def __getattr__(name: str):
    """Build the module-level server app on first access.

    ``uvicorn refgenie.server.main:app`` still works, but the app is no longer
    constructed as an import side effect: ``create_app()`` in server mode mounts
    the MCP sub-app and so needs the ``server`` extras, while ``refgenie dash``
    imports this same module on an install that has only the ``dash`` extra.
    """
    if name == "app":
        application = create_app()
        globals()["app"] = application
        return application
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def run_server(port: int, reload: bool = False):
    """Serve the public API on all interfaces."""
    uvicorn.run("refgenie.server.main:app", host="0.0.0.0", port=port, reload=reload)


def run_local(port: int):
    """Serve the local management app on the loopback interface only.

    The app object is passed rather than an import string on purpose: that keeps
    the caller's real ``Refgenie`` in this process instead of re-resolving the
    user's config in a reloader subprocess.

    Passing an app *object* also pins uvicorn to a single worker, which local
    mode requires: the JobManager holds every pull's and build's state in this
    process's memory, so a second worker would answer half the browser's
    requests with a 404 for a job it had just submitted.
    """
    uvicorn.run(create_app(mode=APP_MODE_LOCAL), host="127.0.0.1", port=port)
