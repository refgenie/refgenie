"""Endpoint hit tracking: collector, ASGI middleware, and persistence handlers."""

from collections import defaultdict
from datetime import datetime
from collections.abc import Callable
from typing import Any

from fastapi import FastAPI
from pydantic import BaseModel, Field
from starlette.types import ASGIApp, Receive, Scope, Send

from refgenie.core import Refgenie
from refgenie.exceptions import MissingStagedAssetError

ParamsKey = frozenset[tuple[str, Any]]


class EndpointStats(BaseModel):
    """Statistics for a specific endpoint + path params combination"""

    total_hits: int = Field(default=0)
    updated_at: datetime = Field(default_factory=datetime.now)

    def increment(self) -> None:
        """Increment the total hits and update the last hit timestamp"""
        self.total_hits += 1
        self.updated_at = datetime.now()

    def reset(self) -> None:
        """Reset the total hits and update the last hit timestamp"""
        self.total_hits = 0
        self.updated_at = datetime.now()


class EndpointCollector:
    def __init__(self):
        # Structure: {route_name: {frozen_path_params: EndpointStats}}
        self._stats: defaultdict[str, defaultdict[ParamsKey, EndpointStats]] = defaultdict(
            lambda: defaultdict(EndpointStats)
        )

    def __repr__(self):
        return f"EndpointCollector({self._stats})"

    def __str__(self):
        return str(self._stats)

    def record_hit(self, route_name: str, path_params: dict[str, Any]) -> None:
        """Record a hit for a specific route and path parameters"""
        # Convert path params to frozenset for dict key (for hashability)
        params_key = frozenset(path_params.items())

        self._stats[route_name][params_key].increment()

    def get_stats(self, route_name: str) -> dict[ParamsKey, EndpointStats]:
        """Get stats for a specific route"""
        if route_name not in self._stats:
            return {}
        return dict(self._stats[route_name])

    def get_total_hits(self, route_name: str) -> int:
        """Get total hits for a specific route"""
        return sum(stats.total_hits for stats in self._stats[route_name].values())

    def reset_stats(self, route_name: str, path_params: frozenset[tuple[str, Any]]) -> None:
        """Reset stats for a specific route"""
        print(f"Resetting stats for {route_name} with {path_params}")
        del self._stats[route_name][path_params]


class EndpointHitCounterMiddleware:
    """
    Pure ASGI middleware to record hits for specific endpoints and path params.

    Implemented as a bare ASGI callable (rather than subclassing Starlette's
    ``BaseHTTPMiddleware``) so it never receives a ``_CachedRequest`` and remains
    compatible across Starlette/FastAPI versions.
    """

    def __init__(self, app: ASGIApp):
        self.app = app

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        # Only act on HTTP requests; pass through lifespan/websocket unchanged.
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return

        # Dispatch first: Starlette's router resolves the route and writes
        # `route`/`path_params` into this same scope dict. Reading them after the
        # fact avoids re-implementing route matching here, and works whether
        # include_router flattens its children (fastapi <0.137) or wraps them in
        # a lazy _IncludedRouter (>=0.137).
        await self.app(scope, receive, send)

        app = scope.get("app")
        route = scope.get("route")
        path_params = scope.get("path_params")
        if app is not None and route is not None and path_params is not None:
            app.state.endpoint_hits_collector.record_hit(
                path_params=path_params, route_name=route.name
            )


HandlerFun = Callable[[EndpointCollector, Refgenie], None]


def persist_archive_download_count(endpoint_hits_collector: EndpointCollector, refgenie: Refgenie):
    """
    Persist the download count for each staged archive to the database.
    Calls `refgenie.stage.increment_download_count` for each archive.
    """
    downaload_stats = endpoint_hits_collector.get_stats("download_archive")
    for path_params, stats in downaload_stats.items():
        if (asset_digest := dict(path_params).get("asset_digest")) is None or stats.total_hits == 0:
            continue
        try:
            refgenie.stage.increment_download_count(
                asset_digest, mode="archive", count=stats.total_hits
            )
        except MissingStagedAssetError as e:
            print(f"Error incrementing download count for {asset_digest}: {e}")


HANDLERS_CATALOG: dict[str, HandlerFun] = {
    "increment_archive_download_count": persist_archive_download_count,
}


def handle_endpoint_collector_hits(app: FastAPI, refgenie: Refgenie):
    """
    Execute all handlers from the catalog to handle the endpoint hits collected
    After the handlers are run, reset the endpoint hits collector.
    """
    endpoint_hits_collector = app.state.endpoint_hits_collector
    for handler_name, handler_fun in HANDLERS_CATALOG.items():
        print(f"Running endpoint hits collector {handler_name=}")
        handler_fun(endpoint_hits_collector, refgenie)
    # reset the handled & unhandled stats
    app.state.endpoint_hits_collector = EndpointCollector()
