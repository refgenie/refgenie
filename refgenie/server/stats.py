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
        params_key = frozenset(path_params.items())

        self._stats[route_name][params_key].increment()

    def get_stats(self, route_name: str) -> dict[ParamsKey, EndpointStats]:
        """Get stats for a specific route"""
        if route_name not in self._stats:
            return {}
        return dict(self._stats[route_name])


class EndpointHitCounterMiddleware:
    """Pure ASGI middleware that records hits per route and path params.

    Not ``BaseHTTPMiddleware``, to stay compatible across Starlette versions.
    """

    def __init__(self, app: ASGIApp):
        self.app = app

    async def __call__(self, scope: Scope, receive: Receive, send: Send) -> None:
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return

        # Mounted sub-apps overwrite scope["app"], so keep the host app (which
        # owns the collector) from before dispatch.
        host_app = scope.get("app")

        # The router writes `route` and `path_params` into scope during dispatch.
        await self.app(scope, receive, send)

        app = host_app or scope.get("app")
        collector = getattr(app.state, "endpoint_hits_collector", None) if app else None
        route = scope.get("route")
        path_params = scope.get("path_params")
        if collector is not None and route is not None and path_params is not None:
            collector.record_hit(path_params=path_params, route_name=route.name)


HandlerFun = Callable[[EndpointCollector, Refgenie], None]


def persist_archive_download_count(endpoint_hits_collector: EndpointCollector, refgenie: Refgenie):
    """Add the collected download counts to each staged archive in the database."""
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
    """Run every handler on the collected hits, then reset the collector."""
    endpoint_hits_collector = app.state.endpoint_hits_collector
    for handler_name, handler_fun in HANDLERS_CATALOG.items():
        print(f"Running endpoint hits collector {handler_name=}")
        handler_fun(endpoint_hits_collector, refgenie)
    app.state.endpoint_hits_collector = EndpointCollector()
