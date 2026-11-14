"""Seqcol serving that follows the store registry at runtime.

refget binds one backend per app (``app.state.backend``), so each store is
served by its own seqcol sub-app. The set of stores is not fixed: every catalog
import may add, change or disable one. So the host app mounts a single
:class:`SeqcolServing` at ``SEQCOL_MOUNT_PATH``. It is an ASGI app that sends
each request to the sub-app for the store the path names, looked up in a
:class:`SeqcolTable`. :meth:`SeqcolServing.refresh` builds a new table when the
enabled store set changes and swaps it in; a request reads the table once, so
it is served wholly by the old set or wholly by the new one.

URL layout: the highest-priority store is served at ``SEQCOL_MOUNT_PATH``,
every other store at ``SEQCOL_MOUNT_PATH/<name>``. With no stores, the
DB-backed refget router is served at ``SEQCOL_MOUNT_PATH``.
"""

import threading
from dataclasses import dataclass

from fastapi import FastAPI
from refget.router import create_refget_router, setup_backend
from starlette.types import Receive, Scope, Send

from refgenie.db.tables import Store, StoreType
from refgenie.logger import logger
from refgenie.server.const import (
    SEQCOL_MOUNT_PATH,
    SEQCOL_SERVICE_ID,
    SEQCOL_SERVICE_NAME,
    SERVICE_CONTACT_URL,
    SERVICE_DOCUMENTATION_URL,
    SERVICE_ORGANIZATION,
)


@dataclass(frozen=True)
class SeqcolService:
    """Where one store's seqcol service is served, as ``/service-info`` lists it."""

    store: str
    url: str
    id: str
    name: str
    store_url: str

    @classmethod
    def for_row(cls, row: Store, *, default: bool) -> "SeqcolService":
        """The service for ``row``; the default store is unqualified."""
        if default:
            return cls(row.name, SEQCOL_MOUNT_PATH, SEQCOL_SERVICE_ID, SEQCOL_SERVICE_NAME, row.url)
        return cls(
            row.name,
            f"{SEQCOL_MOUNT_PATH}/{row.name}",
            f"{SEQCOL_SERVICE_ID}.{row.name}",
            f"{SEQCOL_SERVICE_NAME} ({row.name})",
            row.url,
        )


@dataclass(frozen=True)
class SeqcolMount(SeqcolService):
    """A :class:`SeqcolService` with the sub-app that serves it."""

    app: FastAPI


@dataclass(frozen=True)
class SeqcolTable:
    """One store set's seqcol sub-apps. Never mutated; replaced whole."""

    #: (name, url, type, priority) of each enabled store, in priority order;
    #: a change to any of them rebuilds the table.
    stores: tuple[tuple[str, str, str, int], ...]
    #: One per store, default first. Empty in DB mode.
    mounts: tuple[SeqcolMount, ...]


def _db_seqcol_app(engine) -> FastAPI:
    """The DB-backed seqcol app, used when no store is enabled."""
    app = FastAPI(title=SEQCOL_SERVICE_NAME)
    app.include_router(create_refget_router(sequences=False, pangenomes=False))

    @app.get("/service-info", tags=["Sequence Collections"])
    async def seqcol_service_info():
        return {
            "id": SEQCOL_SERVICE_ID,
            "name": SEQCOL_SERVICE_NAME,
            "type": {"group": "org.ga4gh", "artifact": "refget-seqcol", "version": "1.0.0"},
            "seqcol": {"refget_store": {"enabled": False}},
        }

    setup_backend(app, engine=engine)
    return app


def _store_mount(row: Store, store, *, default: bool) -> SeqcolMount:
    """A seqcol sub-app serving ``store``, the opened store for ``row``."""
    from refget.seqcolapi import create_seqcol_app

    service = SeqcolService.for_row(row, default=default)
    app = create_seqcol_app(
        store_path=row.url,
        remote=row.type == StoreType.remote,
        store_url=row.url,
        service_info_id=service.id,
        service_info_name=service.name,
        organization=SERVICE_ORGANIZATION,
        contact_url=SERVICE_CONTACT_URL,
        documentation_url=SERVICE_DOCUMENTATION_URL,
        # The router's stores are frozen; a store-set change rebuilds them.
        freshness=False,
        cors=False,  # the host app's CORS middleware wraps the dispatcher
        # Bound below to the router's own store object.
        defer_backend=True,
    )
    setup_backend(app, store=store)
    return SeqcolMount(**vars(service), app=app)


def _qualifier(scope: Scope) -> str:
    """The first segment of the path below this mount, or ``""``."""
    path = scope["path"]
    root_path = scope.get("root_path", "")
    if root_path and path.startswith(root_path):
        path = path[len(root_path) :]
    parts = path.split("/", 2)
    return parts[1] if len(parts) > 1 and parts[0] == "" else ""


class SeqcolServing:
    """The ASGI app at ``SEQCOL_MOUNT_PATH``, kept in step with the store registry.

    Mounted sub-apps get no lifespan events, so the host app's lifespan calls
    :meth:`refresh` once to load the stores, and the catalog import calls it
    again after each import.
    """

    def __init__(self, rg):
        self._rg = rg
        self._table: SeqcolTable | None = None
        # Serializes refreshes; requests never take it.
        self._lock = threading.Lock()
        self._db_app: FastAPI | None = None

    def services(self) -> tuple[SeqcolService, ...]:
        """The store services to advertise; empty in DB mode."""
        return self._table.mounts

    def refresh(self) -> bool:
        """Serve the registry's current enabled stores. Returns whether anything changed.

        A changed store set reloads the store router (the v4 API reads it too)
        and builds a new table bound to its stores, then swaps both in. Until
        the swap, requests keep using the old ones. A store that fails to load,
        or a registry that cannot be read, raises, and the old router and
        table stay in place.
        """
        with self._lock:
            current = self._table
            rows = list(self._rg.store.enabled_stores())
            keys = tuple((r.name, r.url, StoreType(r.type).value, r.priority) for r in rows)
            if current is not None and current.stores == keys:
                return False
            if rows:
                # The first build reuses a router that already exists; later
                # ones reload it, since the rows it was built from changed.
                router = (
                    self._rg.store_router if current is None else self._rg.reload_store_router()
                )
                mounts = tuple(
                    _store_mount(row, router.get_store(row.name), default=idx == 0)
                    for idx, row in enumerate(rows)
                )
                self._table = SeqcolTable(stores=keys, mounts=mounts)
            else:
                if current is not None:
                    # The router must stop holding the dropped stores, too.
                    self._rg.reload_store_router()
                if self._db_app is None:
                    self._db_app = _db_seqcol_app(self._rg.database_engine)
                self._table = SeqcolTable(stores=keys, mounts=())
        names = ", ".join(k[0] for k in keys) or "none; DB-backed seqcol"
        logger.info(f"Serving seqcol for stores: {names}")
        return True

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        mounts = self._table.mounts  # read once: one request, one store set
        if scope["type"] in ("http", "websocket"):
            name = _qualifier(scope)
            for mount in mounts[1:]:
                if mount.store == name:
                    # As a Mount would: the store's segment joins root_path. The
                    # scope is updated in place so the hit counter sees the route.
                    scope["root_path"] = scope.get("root_path", "") + "/" + name
                    await mount.app(scope, receive, send)
                    return
        # The default store's sub-app, or the DB-backed one with no stores.
        await (mounts[0].app if mounts else self._db_app)(scope, receive, send)
