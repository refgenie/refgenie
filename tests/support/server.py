"""Server app and TestClient construction, extras guards and route assertions."""

import importlib.util
import re
from collections import Counter, defaultdict
from contextlib import contextmanager
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from tests.support.factories import (
    GENOME,
    GROUP,
    TESTS_DATA_DIR,
    fake_digest,
    make_engine,
    register_fasta,
    seed_staged_asset,
)


def requires_server() -> None:
    """Skip the caller unless the 'server' extras are importable."""
    pytest.importorskip("fastapi", reason="server extras not installed (fastapi)")
    pytest.importorskip("apscheduler", reason="server extras not installed (apscheduler)")


def requires_dash() -> None:
    """Skip the caller unless the 'dash' extras are importable."""
    pytest.importorskip("uvicorn", reason="dash extras not installed (uvicorn)")


#: Mark form of ``requires_server`` for a class or a single test. Calling
#: ``requires_server()`` in a class body runs at import and skips the module.
needs_server = pytest.mark.skipif(
    not all(importlib.util.find_spec(m) for m in ("fastapi", "apscheduler")),
    reason="server extras not installed",
)


def web_assets_present() -> bool:
    """Whether the web UI bundle is built in this checkout.

    The single predicate consulted by both the ``requires_web_assets`` marker
    below and the ``REFGENIE_REQUIRE_WEB_ASSETS`` anti-rot guard in
    ``tests/conftest.py``.
    """
    try:
        from refgenie.server.spa import resolve_web_dist
    except ImportError:  # no server extras, so no web UI to serve
        return False
    return resolve_web_dist() is not None


#: Skip a web-UI test when this checkout has no built bundle -- a fresh
#: checkout has none, CI builds one first. See tests/conftest.py for the
#: REFGENIE_REQUIRE_WEB_ASSETS guard that turns a missing bundle into a hard
#: CI failure instead of a silent skip.
requires_web_assets = pytest.mark.skipif(
    not web_assets_present(),
    reason="web UI assets not built (cd frontend && npm run build)",
)


def make_server_app(rg, *, web_dist=None):
    """The REAL refgenie server app over ``rg``. Requires the 'server' extras.

    Plain function (not a fixture) so fastapi stays a lazy import and the
    mock-only test classes keep running without the server extras installed.
    """
    from refgenie.server.main import create_app

    return create_app(refgenie_instance=rg, web_dist=web_dist)


def make_local_app(rg, *, web_dist=None):
    """The REAL local-mode app (``refgenie dash``) over ``rg``.

    Same factory as the server app, one argument different -- that is the point
    of the mode argument. Nothing here hand-builds a ``FastAPI()``: a test must
    exercise what ships.
    """
    from refgenie.server.const import APP_MODE_LOCAL
    from refgenie.server.main import create_app

    return create_app(mode=APP_MODE_LOCAL, refgenie_instance=rg, web_dist=web_dist)


def make_server_client(rg, *, raise_server_exceptions=True):
    """A (not yet entered) TestClient over the REAL server app for ``rg``."""
    from fastapi.testclient import TestClient

    app = make_server_app(rg)
    return TestClient(app, raise_server_exceptions=raise_server_exceptions)


def make_local_client(rg, *, raise_server_exceptions=True):
    """A (not yet entered) TestClient over the REAL local-mode app for ``rg``."""
    from fastapi.testclient import TestClient

    app = make_local_app(rg)
    # base_url: the local app's Host-header guard (local/security.py) admits
    # loopback names only, so the default "testserver" host would 421 everything.
    return TestClient(
        app, base_url="http://localhost", raise_server_exceptions=raise_server_exceptions
    )


def stub_rgc():
    """A MagicMock Refgenie whose alias manager answers deterministically.

    ``list_genomes``/``list_aliases`` source aliases from ``rgc.alias``, not the
    SQL table, so an empty test database needs a stub to stay deterministic.
    """
    from refgenie.exceptions import MissingAliasError

    rgc = MagicMock()
    rgc.alias.list_all.return_value = []
    rgc.alias.get_for_genome.return_value = []
    rgc.alias.get_for_genomes.side_effect = lambda digests: {d: [] for d in digests}
    rgc.alias.resolve.side_effect = MissingAliasError("stub")
    return rgc


@contextmanager
def serve_refgenie(client_rg, server_rg, *urls):
    """Serve the real app for ``server_rg`` over an in-process ASGI TestClient
    and register it on ``client_rg.servers`` under each URL.

    No sockets, no uvicorn: FastAPI's TestClient IS an httpx2.Client over an
    ASGI transport, and RefgenieserverClient accepts an injected http_client.
    Yields inside the live TestClient context.
    """
    from fastapi.testclient import TestClient

    from refgenie.managers.sources.client import RefgenieserverClient

    with TestClient(make_server_app(server_rg)) as tc:
        for url in urls:
            client_rg.servers.clients[url] = RefgenieserverClient(url, http_client=tc)
        yield


# The route-registration guards walk a FastAPI app's *effective* routes: since
# 0.137, include_router no longer flattens child routes onto app.routes but
# hides them behind an opaque _IncludedRouter, so the prefixed
# path/methods/operation_id must be recovered with FastAPI's own walker.
# `iter_route_contexts` is the supported (non-underscore) walker get_openapi
# itself uses. It is a version-specific internal, so it is imported
# defensively: a FastAPI that does not expose it makes the guards skip (see
# `effective_api_routes`) instead of crashing collection of every module that
# imports from here.


def _normalize_route_path(path):
    """Drop Starlette path converters ({p:path} -> {p}) as OpenAPI does."""
    return re.sub(r"\{([^}:]+):[^}]+\}", r"{\1}", path)


def effective_api_routes(app):
    """Yield an effective view of every APIRoute reachable from `app`.

    As of fastapi 0.137 `include_router` no longer copies child routes onto
    `app.routes`; it appends one opaque `_IncludedRouter` and computes the
    effective (prefixed) routes lazily. `iter_route_contexts` is the walker
    `get_openapi` itself uses; each yielded `RouteContext` exposes the effective
    path, methods and operation id, delegating to the underlying route for
    routes declared directly on the app (which need no prefixing).

    The duplicate-registration bug these guards exist to catch is invisible in
    the public OpenAPI schema -- including a router twice yields a byte-identical
    spec -- so there is no public alternative to this walker. If the running
    FastAPI does not expose it, skip rather than assert on a partial view.
    """
    from fastapi.routing import APIRoute

    try:
        from fastapi.routing import iter_route_contexts
    except ImportError:  # pragma: no cover - depends on FastAPI internals
        pytest.skip(
            "FastAPI does not expose routing.iter_route_contexts; cannot walk "
            "effective routes to guard against duplicate registration"
        )
    for ctx in iter_route_contexts(app.router.routes):
        if isinstance(ctx.original_route, APIRoute):
            yield ctx


def route_keys(app):
    """Yield (method, path) for every APIRoute on `app`."""
    for route in effective_api_routes(app):
        for method in sorted(route.methods or []):
            yield method, _normalize_route_path(route.path)


def assert_no_duplicate_routes(app, label):
    counts = Counter(route_keys(app))
    duplicates = {key: n for key, n in counts.items() if n > 1}
    assert not duplicates, f"{label}: route(s) registered more than once: " + ", ".join(
        f"{m} {p} (x{n})" for (m, p), n in sorted(duplicates.items())
    )


def assert_unique_operation_ids(app, label):
    """Duplicate operationIds break generated clients even when paths differ."""
    seen = defaultdict(list)
    for route in effective_api_routes(app):
        if route.include_in_schema:
            seen[route.operation_id or route.unique_id].append(route.path)
    duplicates = {op: paths for op, paths in seen.items() if len(paths) > 1}
    assert not duplicates, f"{label}: duplicate operationIds: {duplicates}"


def server_with_asset(
    tmp_path,
    *,
    genome_digest=None,
    alias="test_alias",
    assets=(),
    fixtures_path=None,
):
    """A server-side Refgenie with one genome and zero or more seeded assets.

    ``assets`` is a sequence of dicts consumed by :func:`seed_staged_asset`,
    each optionally carrying a ``key`` used in the returned digest mapping
    (defaulting to the asset group name). Returns ``(rgc, {key: asset_digest})``.
    The genome digest defaults to ``fake_digest("test_digest_12345")``.
    """
    from refgenie import Refgenie

    genome_digest = genome_digest or fake_digest("test_digest_12345")

    tmp_path = Path(tmp_path)
    rgc = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rgc.database.init(
        genome_folder=tmp_path / "genomes",
        genome_stage_folder=tmp_path / "archives",
    )
    register_fasta(rgc, fixtures_path)
    rgc.genome.add(
        digest=genome_digest,
        description="Test genome",
        alias_names=[alias] if alias else [],
    )

    digests = {}
    for spec in assets:
        spec = dict(spec)
        key = spec.pop("key", spec["asset_group_name"])
        digests[key] = seed_staged_asset(rgc, genome_digest=genome_digest, **spec)
    return rgc, digests


def make_server_client_world(
    engine,
    tmp_path,
    fixtures_path=None,
    *,
    stage=True,
    register_client_fasta=True,
):
    """A server with a built fasta asset plus a fresh client; ``(server, client)``.

    Note the asymmetry: the *server* engine is created internally,
    and only the *client* uses the passed-in ``engine`` -- callers assert
    against the client's database.

    Serve the pair with ``serve_refgenie(client_rg, server_rg, *urls)``.
    """
    from refgenie import Refgenie

    tmp_path = Path(tmp_path)
    fixtures_path = Path(fixtures_path) if fixtures_path is not None else TESTS_DATA_DIR

    server_genome_folder = tmp_path / "server_genomes"
    server_archive_folder = tmp_path / "server_archives"
    server_rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    server_rg.database.init(
        genome_folder=server_genome_folder,
        genome_stage_folder=server_archive_folder,
    )
    register_fasta(server_rg, fixtures_path)
    server_rg.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=[GENOME],
        description="rCRSd mitochondrial reference",
    )
    server_rg.build.run(
        recipe_name="fasta",
        genome_alias=GENOME,
        asset_group_name=GROUP,
        asset_name="default",
    )
    if stage:
        assets = list(server_rg.asset.list_assets())
        assert len(assets) >= 1
        server_rg.stage.create(assets[0], server_genome_folder, server_archive_folder)

    client_rg = Refgenie(database_engine=engine, suppress_migrations=True)
    client_rg.database.init(genome_folder=tmp_path / "client_genomes")
    if register_client_fasta:
        register_fasta(client_rg, fixtures_path)

    return server_rg, client_rg
