"""Shared test helpers: plain functions, classes and constants.

Two homes, one rule:

* ``tests/helpers.py`` (this module) -- every plain importable function, class
  and constant. **No fixtures.**
* ``tests/conftest.py`` -- fixtures and pytest hooks only. It imports from this
  module whatever its fixtures need.

Top-level imports here are limited to the stdlib, ``pytest`` and
``unittest.mock``. Imports of ``refgenie``, ``sqlmodel``, ``fastapi`` and
friends go *inside* the functions that need them, for two reasons:

* ``tests/e2e/cli_compat/conftest.py`` imports the CLI launcher from here and
  that suite is contractually library-free (it also runs against the Rust
  binary via ``REFGENIE_MODE=rust``);
* the lazy-import pattern is what keeps the mock-only test classes running
  without the ``server`` extras installed.
"""

import hashlib
import json
import os
import shutil
import socket
import subprocess
import sys
import threading
import time
from contextlib import ExitStack, contextmanager
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pytest

# ---------------------------------------------------------------------------
# Paths and well-known names
# ---------------------------------------------------------------------------

TESTS_DATA_DIR = Path(__file__).parent / "data"
RCRSD_FASTA = TESTS_DATA_DIR / "rCRSd.fa"
DEMO_FASTA = TESTS_DATA_DIR / "demo.fa"

# The standard built world: one rCRSd genome with one fasta asset named "test".
GENOME = "rCRSd"
GROUP = "fasta"
ASSET = "test"

#: The in-process ASGI server URLs used by ``serve_refgenie``.
MOCK_SERVER_URL = "http://mock-refgenie-server"
MOCK_MIRROR_URL = "http://mock-refgenie-mirror"

#: Sentinel: leave the key out of a payload / the attribute unset on a mock.
#: An *unset* MagicMock attribute returns a truthy MagicMock, which is
#: behaviourally different from ``[]`` -- some tests depend on that difference.
OMIT = object()


# ---------------------------------------------------------------------------
# Optional-extras guards
# ---------------------------------------------------------------------------


def requires_server() -> None:
    """Skip the caller unless the 'server' extras are importable."""
    pytest.importorskip("fastapi", reason="server extras not installed (fastapi)")
    pytest.importorskip("apscheduler", reason="server extras not installed (apscheduler)")


def requires_dash() -> None:
    """Skip the caller unless the 'dash' extras are importable."""
    pytest.importorskip("uvicorn", reason="dash extras not installed (uvicorn)")


def web_assets_present() -> bool:
    """Whether the web UI bundle is built in this checkout.

    The single predicate consulted by both the ``requires_web_assets`` marker
    below and the ``REFGENIE_REQUIRE_WEB_ASSETS`` anti-rot guard in
    ``tests/conftest.py``.
    """
    from refgenie.server.spa import resolve_web_dist

    return resolve_web_dist() is not None


#: Skip a web-UI test when this checkout has no built bundle -- a fresh
#: checkout has none, CI builds one first. See tests/conftest.py for the
#: REFGENIE_REQUIRE_WEB_ASSETS guard that turns a missing bundle into a hard
#: CI failure instead of a silent skip.
requires_web_assets = pytest.mark.skipif(
    not web_assets_present(),
    reason="web UI assets not built (cd frontend && npm run build)",
)


# ---------------------------------------------------------------------------
# Engines and Refgenie worlds
# ---------------------------------------------------------------------------


def make_engine():
    """A fresh in-memory SQLite engine (StaticPool, thread-check off).

    Plain function so module-level constants and multi-engine tests
    (server + client pairs, dual-world tests) can build engines directly;
    single-engine tests should use the ``engine`` fixture instead.
    """
    from sqlmodel import create_engine
    from sqlmodel.pool import StaticPool

    return create_engine(
        "sqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
        echo=False,
    )


def register_fasta(r, fixtures_path=None):
    """Register fasta asset class and recipe from test data files.

    Idempotent -- safe to call multiple times on the same database.
    """
    from refgenie.exceptions import AssetClassExistsError, RecipeExistsError

    fixtures_path = Path(fixtures_path) if fixtures_path is not None else TESTS_DATA_DIR
    try:
        r.asset_class.add(fixtures_path / "fasta_asset_class.yaml")
    except AssetClassExistsError:
        pass
    try:
        r.recipe.add(fixtures_path / "fasta_asset_recipe.yaml")
    except RecipeExistsError:
        pass


def samtools_available() -> bool:
    """Check if samtools is available on PATH.

    If you need samtools from bulker, activate the crate before running tests:
        bulker activate refgenie/refgenie
        pytest
    """
    return shutil.which("samtools") is not None


def build_rcrsd(r, asset_name=ASSET):
    """Build the standard fasta asset for the rCRSd genome."""
    return r.build_asset(
        recipe_name="fasta",
        genome_name=GENOME,
        asset_group_name=GROUP,
        asset_name=asset_name,
    )


def stage_rcrsd(r, asset_name=ASSET):
    """Stage the built rCRSd fasta asset (default serving_modes); return it."""
    genome_digest = r.alias.resolve(GENOME)
    built_asset = r.asset.get(
        genome_digest=genome_digest,
        asset_group_name=GROUP,
        asset_name=asset_name,
    )
    r.stage.create(
        asset=built_asset,
        genome_folder=r.genome_folder,
        genome_stage_folder=r.genome_stage_folder,
    )
    return built_asset


def make_built_refgenie(root, fixtures_path=None, *, asset_name=ASSET, build=True, stage=False):
    """A Refgenie on a real filesystem under ``root`` with the rCRSd genome.

    Factory form of the refgenie_fs/refgenie_built/staged_refgenie fixture
    chain, for tests that need two independent worlds or a named root.
    """
    from refgenie import Refgenie

    fixtures_path = Path(fixtures_path) if fixtures_path is not None else TESTS_DATA_DIR
    root = Path(root)
    r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    r.init(
        genome_folder=root / "genomes",
        genome_stage_folder=root / "archives",
    )
    register_fasta(r, fixtures_path)
    r.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=[GENOME],
        description="rCRSd genome",
    )
    if build:
        build_rcrsd(r, asset_name=asset_name)
    if stage:
        stage_rcrsd(r, asset_name=asset_name)
    return r


def make_server_rgc(tmp_path, fixtures_path=None, *, genomes=()):
    """A Refgenie over a fresh in-memory catalog with the fasta definitions.

    ``genomes`` is an iterable of ``(digest, alias_names)`` to register.
    """
    from refgenie import Refgenie

    rgc = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rgc.init(genome_folder=Path(tmp_path) / "genomes")
    register_fasta(rgc, fixtures_path)
    for digest, alias_names in genomes:
        rgc.genome.add(digest=digest, description="Test genome", alias_names=list(alias_names))
    return rgc


# ---------------------------------------------------------------------------
# Server app / TestClient plumbing
# ---------------------------------------------------------------------------


def _override_session(app, engine):
    """Point ``get_db_session`` at ``engine`` on an already-built app."""
    from sqlmodel import Session

    from refgenie.server.dependencies import get_db_session

    def _session():
        with Session(engine) as session:
            yield session

    app.dependency_overrides[get_db_session] = _session


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


def make_server_client(rg, *, engine=None, raise_server_exceptions=True):
    """A (not yet entered) TestClient over the REAL server app for ``rg``."""
    from fastapi.testclient import TestClient

    app = make_server_app(rg)
    if engine is not None:
        _override_session(app, engine)
    return TestClient(app, raise_server_exceptions=raise_server_exceptions)


def make_local_client(rg, *, engine=None, raise_server_exceptions=True):
    """A (not yet entered) TestClient over the REAL local-mode app for ``rg``.

    ``engine`` optionally overrides ``get_db_session`` with sessions on that
    engine, for tests that want an empty in-memory database under a stubbed
    Refgenie.
    """
    from fastapi.testclient import TestClient

    app = make_local_app(rg)
    if engine is not None:
        _override_session(app, engine)
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
    rgc.alias.resolve.side_effect = MissingAliasError("stub")
    return rgc


@contextmanager
def serve_refgenie(client_rg, server_rg, *urls):
    """Serve the real app for ``server_rg`` over an in-process ASGI TestClient
    and register it on ``client_rg``'s source manager under each URL.

    No sockets, no uvicorn: FastAPI's TestClient IS an httpx.Client over an
    ASGI transport, and RefgenieserverClient accepts an injected http_client.
    Yields inside the live TestClient context.
    """
    from fastapi.testclient import TestClient

    from refgenie.managers.sources.client import RefgenieserverClient

    with TestClient(make_server_app(server_rg)) as tc:
        for url in urls:
            client_rg.sources._server_clients[url] = RefgenieserverClient(url, http_client=tc)
        yield


# ---------------------------------------------------------------------------
# CLI subprocess launcher
# ---------------------------------------------------------------------------

#: The one way to launch the Python CLI in a subprocess.
PY_CLI_ARGV = [
    sys.executable,
    "-c",
    "from refgenie.cli.main import main; import sys; main(sys.argv[1:])",
]


def cli_argv(*args: str) -> list[str]:
    """``PY_CLI_ARGV`` plus ``args``, as a fresh list."""
    return PY_CLI_ARGV + [str(a) for a in args]


def _merged_env(env: dict | None) -> dict:
    """``os.environ`` updated with ``env``; a ``None`` value *unsets* the name."""
    full_env = os.environ.copy()
    for key, value in (env or {}).items():
        if value is None:
            full_env.pop(key, None)
        else:
            full_env[key] = value
    return full_env


def run_cli(
    *args: str,
    env: dict | None = None,
    check: bool = True,
    timeout: int = 120,
    text: bool = True,
) -> subprocess.CompletedProcess:
    """Run the refgenie CLI in a subprocess, capturing output.

    ``env`` is merged into ``os.environ``. ``check`` raises
    ``CalledProcessError`` *after* the output is captured (rather than using
    subprocess's own ``check=``), because that is what ``run_refgenie``'s
    callers depend on.
    """
    cmd = cli_argv(*args)
    result = subprocess.run(
        cmd, capture_output=True, text=text, env=_merged_env(env), timeout=timeout
    )
    if check and result.returncode != 0:
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stdout, result.stderr)
    return result


def popen_cli(*args: str, env: dict | None = None, text: bool = False) -> subprocess.Popen:
    """Launch the refgenie CLI as a background subprocess with piped output.

    ``env`` is merged into ``os.environ``.
    """
    return subprocess.Popen(
        cli_argv(*args),
        env=_merged_env(env),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=text,
    )


def find_free_port() -> int:
    """Find a free TCP port on localhost."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def wait_for_server(host: str, port: int, timeout: float = 30.0) -> None:
    """Block until ``host:port`` accepts a connection, or raise RuntimeError."""
    start = time.time()
    while time.time() - start < timeout:
        try:
            with socket.create_connection((host, port), timeout=0.5):
                return
        except (ConnectionRefusedError, OSError):
            time.sleep(0.2)
    raise RuntimeError(f"Server failed to start at {host}:{port}")


# ---------------------------------------------------------------------------
# Direct DB seeding for server tests
# ---------------------------------------------------------------------------


def create_asset_class(session, name, serving_modes, version="0.1.0"):
    """Create an AssetClass with the given serving_modes."""
    from refgenie.db.tables import AssetClass

    ac = AssetClass(
        name=name,
        version=version,
        description=f"Test {name} asset class",
        serving_modes=list(serving_modes),
    )
    session.add(ac)
    session.commit()
    session.refresh(ac)
    return ac


def create_asset_on_disk(genome_folder, genome_digest, group_name, asset_name, files) -> Path:
    """Create asset files on disk and return the asset directory path.

    ``files`` is a ``{name: str | bytes}`` mapping.
    """
    asset_path = Path(genome_folder) / "data" / genome_digest / group_name / asset_name
    asset_path.mkdir(parents=True, exist_ok=True)
    for fname, content in files.items():
        if isinstance(content, bytes):
            (asset_path / fname).write_bytes(content)
        else:
            (asset_path / fname).write_text(content)
    return asset_path


def seed_asset(
    session,
    asset_group_id,
    asset_path,
    *,
    name="default",
    description="Test asset",
    size=0,
    is_default=True,
    seed_asset_name=True,
) -> str:
    """Insert a content-addressed Asset (plus its AssetName row); return the digest.

    The content digest is not computed by an insert-time handler, so callers
    that build a directory and seed the DB directly must compute it here.
    ``seed_asset_name=False`` deliberately creates an Asset with *no* AssetName.
    """
    from refgenie.db.tables import Asset, AssetName
    from refgenie.utils.build import get_dir_digest

    digest = get_dir_digest(asset_path)
    session.add(
        Asset(
            digest=digest,
            name=name,
            description=description,
            asset_group_id=asset_group_id,
            size=size,
            path=str(asset_path),
        )
    )
    if seed_asset_name:
        session.add(
            AssetName(
                name=name,
                asset_group_id=asset_group_id,
                asset_digest=digest,
                is_default=is_default,
            )
        )
    session.commit()
    return digest


def seed_staged_asset(
    rgc,
    *,
    genome_digest,
    asset_group_name,
    asset_class,
    files,
    asset_name="default",
    stage_modes=(),
    tarball_bytes=None,
    description="Test asset",
    group_description=None,
    size=0,
    is_default=True,
    seed_asset_name=True,
) -> str:
    """Create an asset group + on-disk asset + DB rows; return the asset digest.

    ``asset_class`` is either an existing class name (looked up) or a
    ``(name, serving_modes[, version])`` tuple to synthesize one.
    ``stage_modes`` is a subset of ``{"file", "archive"}``. ``"file"`` creates
    the stage symlink plus a row with ``directory_contents``; ``"archive"``
    writes ``tarball_bytes`` if given (and records their real sha-256), or
    records a placeholder ``tarball_digest``/``tarball_size``.
    """
    from sqlmodel import Session

    from refgenie.db.tables import AssetGroup, StagedAsset

    engine = rgc.database_engine
    genome_folder = Path(str(rgc.genome_folder))
    stage_folder = Path(str(rgc.genome_stage_folder))
    contents = sorted(files)

    asset_path = create_asset_on_disk(
        genome_folder, genome_digest, asset_group_name, asset_name, files
    )

    with Session(engine) as session:
        if isinstance(asset_class, str):
            asset_class_id = rgc.asset_class.get(asset_class).id
        else:
            ac_name, serving_modes, *rest = asset_class
            asset_class_id = create_asset_class(
                session, ac_name, serving_modes, version=rest[0] if rest else "0.1.0"
            ).id

        group = AssetGroup(
            name=asset_group_name,
            description=group_description or f"Test {asset_group_name} asset group",
            genome_digest=genome_digest,
            asset_class_id=asset_class_id,
        )
        session.add(group)
        session.commit()
        session.refresh(group)

        digest = seed_asset(
            session,
            group.id,
            asset_path,
            name=asset_name,
            description=description,
            size=size,
            is_default=is_default,
            seed_asset_name=seed_asset_name,
        )

        if "file" in stage_modes:
            link_dir = stage_folder / genome_digest / asset_group_name / asset_name
            link_dir.parent.mkdir(parents=True, exist_ok=True)
            link_dir.symlink_to(asset_path)
            session.add(
                StagedAsset(
                    asset_digest=digest,
                    mode="file",
                    directory_contents=contents,
                    tarball_digest=None,
                    tarball_size=None,
                )
            )
            session.commit()

        if "archive" in stage_modes:
            if tarball_bytes is not None:
                archive_dir = stage_folder / genome_digest / asset_group_name
                archive_dir.mkdir(parents=True, exist_ok=True)
                (archive_dir / f"{digest}.tgz").write_bytes(tarball_bytes)
                tarball_digest = hashlib.sha256(tarball_bytes).hexdigest()
                tarball_size = len(tarball_bytes)
            else:
                tarball_digest = "ab" * 32
                tarball_size = 2048
            session.add(
                StagedAsset(
                    asset_digest=digest,
                    mode="archive",
                    directory_contents=contents,
                    tarball_digest=tarball_digest,
                    tarball_size=tarball_size,
                )
            )
            session.commit()

    return digest


def server_with_asset(
    tmp_path,
    *,
    genome_digest="test_digest_12345",
    alias="test_alias",
    assets=(),
    fixtures_path=None,
):
    """A server-side Refgenie with one genome and zero or more seeded assets.

    ``assets`` is a sequence of dicts consumed by :func:`seed_staged_asset`,
    each optionally carrying a ``key`` used in the returned digest mapping
    (defaulting to the asset group name). Returns ``(rgc, {key: asset_digest})``.
    """
    from refgenie import Refgenie

    tmp_path = Path(tmp_path)
    rgc = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rgc.init(
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

    Note the deliberate asymmetry: the *server* engine is created internally,
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
    server_rg.init(
        genome_folder=server_genome_folder,
        genome_stage_folder=server_archive_folder,
    )
    register_fasta(server_rg, fixtures_path)
    server_rg.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=[GENOME],
        description="rCRSd mitochondrial reference",
    )
    server_rg.build_asset(
        recipe_name="fasta",
        genome_name=GENOME,
        asset_group_name=GROUP,
        asset_name="default",
    )
    if stage:
        assets = list(server_rg.asset.list_assets())
        assert len(assets) >= 1
        server_rg.stage.create(assets[0], server_genome_folder, server_archive_folder)

    client_rg = Refgenie(database_engine=engine, suppress_migrations=True)
    client_rg.init(genome_folder=tmp_path / "client_genomes")
    if register_client_fasta:
        register_fasta(client_rg, fixtures_path)

    return server_rg, client_rg


# ---------------------------------------------------------------------------
# Mock server clients and remote sources
# ---------------------------------------------------------------------------


def mock_server_client(
    *,
    server_url="http://test.example.com",
    asset_group_name="bowtie2_index",
    genome_digest="genome_digest_001",
    asset_group_id=1,
    asset_name="default",
    asset_digest="abc123digest",
    serving_modes=("archive",),
    seek_keys=OMIT,
    is_default=OMIT,
    alias_digest=OMIT,
    staged_items=(),
    file_list=(),
    relations=None,
    aliases=OMIT,
    download_side_effect=None,
) -> MagicMock:
    """A MagicMock ServerClient for the pull/seekr dispatch tests.

    Any argument left at :data:`OMIT` omits that key from the payload (or, for
    ``aliases``, leaves the attribute *unset* on the mock -- an unset attribute
    returns a truthy MagicMock, which is not the same as ``[]``).

    ``seek_keys`` populates the asset payload's ``seek_keys`` (seekr reads the
    default/named seek key straight off it). ``alias_digest`` makes
    ``client.get`` answer the read-only alias->digest endpoint with that digest,
    for remote resolution of a non-local genome.
    """
    client = MagicMock()
    client.server_url = server_url

    group = {"id": asset_group_id, "name": asset_group_name}
    if genome_digest is not OMIT:
        group["genome_digest"] = genome_digest
    client.get_asset_groups.return_value = [group]

    asset = {"id": 1, "name": asset_name, "digest": asset_digest}
    if asset_group_id is not OMIT:
        asset["asset_group_id"] = asset_group_id
    if serving_modes is not OMIT:
        asset["serving_modes"] = list(serving_modes)
    if seek_keys is not OMIT:
        asset["seek_keys"] = [dict(sk) for sk in seek_keys]
    if is_default is not OMIT:
        asset["is_default"] = is_default
    client.get_assets.return_value = [asset]

    client.get_staged_assets.return_value = list(staged_items)
    client.get_asset_file_list.return_value = list(file_list)
    _relations_payload = dict(relations) if relations is not None else {"parents": [], "children": []}
    if alias_digest is not OMIT:
        from refgenie.managers.sources.api_ids import API_ID_ALIAS_DIGEST

        def _get(operation_id, params=None, url_format_params=None):
            if operation_id == API_ID_ALIAS_DIGEST:
                return {"digest": alias_digest}
            return _relations_payload

        client.get.side_effect = _get
    else:
        client.get.return_value = _relations_payload
    if aliases is not OMIT:
        client.get_all_aliases.return_value = list(aliases)
    if download_side_effect is not None:
        client.download_with_progress.side_effect = download_side_effect
    return client


def real_client_no_init(url: str = "http://test-server"):
    """A REAL RefgenieserverClient built via ``__new__`` (no network in __init__).

    Deliberately separate from :func:`mock_server_client`: the point is that
    ``get_paginated`` and friends run their real code.
    """
    from refgenie.managers.sources.client import RefgenieserverClient

    client = RefgenieserverClient.__new__(RefgenieserverClient)
    client.server_url = url
    client.openapi_endpoint = "/openapi.json"
    client.api_version = "v4"
    client._http_client = None
    return client


class MockRemoteSource:
    """A faithful in-memory RemoteGenomeSource for testing."""

    def __init__(
        self,
        collections=None,
        store_url="https://store.example.com",
        aliases=None,
        collection_aliases=None,
        fhr=None,
    ):
        self._collections = collections or {}
        self._store_url = store_url
        self._aliases = aliases or {}
        # digest -> [(namespace, alias), ...]
        self._collection_aliases = collection_aliases or {}
        # digest -> camelCase FHR dict
        self._fhr = fhr or {}

    def verify_collection(self, digest):
        return self._collections.get(digest)

    def list_collections(self, page=0, page_size=100):
        return {"results": [{"digest": d} for d in self._collections]}

    def resolve_alias(self, alias):
        return self._aliases.get(alias)

    def get_collection_aliases(self, digest):
        return list(self._collection_aliases.get(digest, []))

    def get_collection_fhr(self, digest):
        return self._fhr.get(digest)

    @property
    def store_url(self):
        return self._store_url

    @property
    def url(self):
        return "https://mock.example.com"


def assert_mock_remote_source_conforms() -> None:
    """The test double must actually satisfy the protocol, or every test using
    it proves nothing about real sources."""
    from refgenie.managers.sources.genomes import RemoteGenomeSource

    assert isinstance(MockRemoteSource(), RemoteGenomeSource)


@contextmanager
def mocked_puller(
    rg,
    mock_client,
    *,
    alias=GENOME,
    subscriptions=("http://test.example.com",),
    remote_source=None,
    mock_genome=True,
    mock_download_modes=True,
    mock_asset_writes=True,
):
    """Point ``rg``'s AssetPuller at ``mock_client`` with the usual patches.

    Yields a namespace with ``puller``, ``client`` and (when patched) the
    ``archive``, ``file`` and ``ensure`` mocks so assertions still work.

    * ``mock_genome=False`` leaves ``_ensure_genome_exists`` (and the asset
      ``exists`` check) real -- required by the rollback tests, which are
      *about* genome creation.
    * ``mock_download_modes=False`` lets the real ``_pull_*_mode`` run (and
      raise).
    * ``mock_asset_writes=False`` leaves the cross-manager writes real.
    """
    puller = rg.asset._asset_puller
    ns = SimpleNamespace(puller=puller, client=mock_client, archive=None, file=None, ensure=None)

    with ExitStack() as stack:
        stack.enter_context(
            patch.object(puller._sources, "get_subscriptions", return_value=list(subscriptions))
        )
        stack.enter_context(
            patch.object(puller._sources, "get_server_client", return_value=mock_client)
        )
        if remote_source is not None:
            stack.enter_context(
                patch(
                    "refgenie.managers.asset.genome_bootstrap.make_source",
                    return_value=remote_source,
                )
            )
        if mock_genome:
            ensure = stack.enter_context(patch.object(puller, "_ensure_genome_exists"))
            ensure.return_value = MagicMock(
                success=True,
                genome_digest=rg.alias.resolve(alias),
                created_genome=False,
                created_alias=False,
                alias_name=alias,
            )
            ns.ensure = ensure
            stack.enter_context(
                patch.object(puller._asset_manager, "exists", return_value=False)
            )
        if mock_download_modes:
            ns.archive = stack.enter_context(patch.object(puller, "_pull_archive_mode"))
            ns.file = stack.enter_context(patch.object(puller, "_pull_file_mode"))
        if mock_asset_writes:
            stack.enter_context(
                patch.object(puller._asset_manager, "add_from_path", return_value=MagicMock())
            )
            stack.enter_context(
                patch.object(puller._asset_manager, "get_default", return_value=None)
            )
            stack.enter_context(patch.object(puller._asset_manager, "set_default"))
            stack.enter_context(patch.object(puller, "_create_symlinks_for_alias"))
        yield ns


# ---------------------------------------------------------------------------
# The standard built world: accessors and snapshots
# ---------------------------------------------------------------------------


def fasta_asset(r, name=ASSET, group=GROUP, genome=GENOME):
    """The built fasta asset ``name`` for ``genome``."""
    return r.asset.get(
        genome_digest=r.alias.resolve(genome),
        asset_group_name=group,
        asset_name=name,
    )


def genome_digest(r, genome=GENOME) -> str:
    return r.alias.resolve(genome)


def only_asset(r):
    """The single asset in a one-asset world (index-based, deliberately)."""
    return list(r.asset.list_assets())[0]


def content_dir(r) -> Path:
    """The digest-addressed content directory of the one built asset."""
    return r.genome_folder / only_asset(r).path


def build_flag(r, group=GROUP, asset=ASSET, genome=GENOME) -> Path:
    """The completion flag at the path advertised to snakemake."""
    template = r.get_asset_build_target_template(group, asset)
    return Path(str(template).replace("{genome_name}", genome))


def build_dir(r, asset_name=ASSET, group=GROUP, genome=GENOME) -> Path:
    from refgenie.utils.build import get_build_dir

    return get_build_dir(
        genome_folder=r.genome_folder,
        genome_name=genome,
        asset_group_name=group,
        asset_name=asset_name,
    )


def boom(*_args, **_kwargs):
    """Monkeypatch target that kills a write path mid-flight."""
    raise RuntimeError("killed before the commit")


def assets_rows(r):
    """Every Asset ORM row."""
    from sqlmodel import Session, select

    from refgenie.db.tables import Asset

    with Session(r.database_engine) as session:
        return session.exec(select(Asset)).unique().all()


def asset_name_rows(r):
    """Every AssetName ORM row."""
    from sqlmodel import Session, select

    from refgenie.db.tables import AssetName

    with Session(r.database_engine) as session:
        return session.exec(select(AssetName)).all()


def db_snapshot(r) -> dict:
    """Every table removal is supposed to touch, as plain comparable values."""
    from sqlmodel import Session, select

    from refgenie.db.tables import AssetGroup, Genome

    with Session(r.database_engine) as session:
        groups = sorted(
            (g.genome_digest, g.name) for g in session.exec(select(AssetGroup)).unique().all()
        )
        genomes = sorted(g.digest for g in session.exec(select(Genome)).unique().all())
    return {
        "genomes": genomes,
        "groups": groups,
        "assets": sorted(a.digest for a in assets_rows(r)),
        "asset_names": sorted(n.name for n in asset_name_rows(r)),
    }


def disk_snapshot(r) -> list[str]:
    """Every path under the genome folder, relative and sorted."""
    root = r.genome_folder
    return sorted(p.relative_to(root).as_posix() for p in root.rglob("*"))


# ---------------------------------------------------------------------------
# Asset construction from plain files
# ---------------------------------------------------------------------------


def add_asset_from_files(
    r,
    *,
    rel_dir,
    files,
    asset_class_name,
    asset_group_name,
    asset_name,
    custom_seek_keys=None,
    genome=GENOME,
    contents="payload\n",
):
    """Materialize ``files`` under ``genome_folder/rel_dir`` and register them.

    ``{genome}`` in a filename is expanded to the genome digest, matching how
    assets are named on disk. Nested paths are created with ``parents=True``.
    """
    digest = r.alias.resolve(genome)
    asset_dir = r.genome_folder / rel_dir
    asset_dir.mkdir(parents=True, exist_ok=True)
    for name in files:
        path = asset_dir / name.replace("{genome}", digest)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents)
    return r.add(
        genome_name=genome,
        asset_group_name=asset_group_name,
        asset_name=asset_name,
        path=Path(rel_dir),
        asset_class_name=asset_class_name,
        custom_seek_keys=custom_seek_keys,
    )


def stage_copy(r, asset, dirname, mutate=b"") -> tuple[Path, Path]:
    """Copy an asset's content dir to a staging dir under the genome folder.

    Appending ``mutate`` bytes to the .fa makes the content digest differ.
    Returns ``(absolute_staging_path, path_relative_to_genome_folder)``; the
    relative form is what ``add_from_path`` expects.
    """
    staging = r.genome_folder / dirname
    shutil.copytree(r.genome_folder / asset.path, staging, symlinks=True)
    if mutate:
        fa = next(staging.glob("*.fa"))
        fa.write_bytes(fa.read_bytes() + mutate)
    return staging, Path(dirname)


def make_command_values(**overrides):
    """A minimal BuildCommandValues, with every field overridable."""
    from refgenie.models import BuildCommandValues

    defaults = dict(
        genome_digest="a" * 32,
        asset_group_name="test_group",
        genome_folder=Path("/tmp/genomes"),
        params={},
        files={},
        assets=None,
        custom_seek_keys={},
    )
    defaults.update(overrides)
    return BuildCommandValues(**defaults)


# ---------------------------------------------------------------------------
# Web UI test helpers (actions / jobs / bridge)
# ---------------------------------------------------------------------------

#: Every state-changing web route carries this anti-CSRF header.
ACTION_HEADERS = {"x-refgenie-action": "1"}


def act(client, method, path, json=None, headers=None):
    """Issue a state-changing request with the action header attached."""
    return client.request(method, path, json=json, headers={**ACTION_HEADERS, **(headers or {})})


class FakeRefgenie:
    """Just enough of a Refgenie for the jobs router's preflight/stage check."""

    def __init__(self, genome_stage_folder=None):
        self.genome_stage_folder = genome_stage_folder


def job_result_ok():
    """The typed ``JobResult`` every fake runner returns."""
    from refgenie.server.jobs.schemas import JobResult

    return JobResult(asset_digest="d" * 64, registry_path="rCRSd/fasta:default")


def instant_runner(ctx):
    """A job runner that succeeds immediately."""
    return job_result_ok()


def web_stub_rgc():
    """A MagicMock Refgenie safe to drive every actions route end to end.

    Return values are plain strings/lists so job results and ActionResult
    payloads serialize; the mocked asset satisfies the runners' ``_asset_result``.
    """
    rgc = stub_rgc()
    asset = MagicMock()
    asset.digest = "assetdigest123"
    asset.name = "test"
    asset.registry_path = "rCRSd/fasta:test"
    asset.asset_group.name = "fasta"
    asset.asset_group.genome.digest = "genomedigest123"
    rgc.pull.return_value = asset
    rgc.build_asset.return_value = asset
    rgc.initialize_and_build.return_value = ("genomedigest123", True)
    rgc.asset.remove_by_digest.return_value = "rCRSd/fasta:test"
    rgc.configuration.get_server_subscriptions.return_value = []
    rgc.preflight_build.return_value = {"ok": True, "errors": [], "resolved": {}}
    return rgc


# ---------------------------------------------------------------------------
# Jobs: params, fake runners, the bare jobs app, and an SSE reader
# ---------------------------------------------------------------------------

#: Minimal valid job bodies, as the HTTP API takes them.
PULL_PARAMS = {"asset_group_name": "fasta", "genome_name": "rCRSd"}
BUILD_PARAMS = {"recipe_name": "fasta", "genome_name": "rCRSd", "asset_group_name": "fasta"}


def pull_params(asset_group_name="fasta", genome_name="rCRSd", **kw):
    from refgenie.server.jobs.schemas import PullJobParams

    return PullJobParams(asset_group_name=asset_group_name, genome_name=genome_name, **kw)


def build_params(asset_group_name="fasta", genome_name="rCRSd", **kw):
    from refgenie.server.jobs.schemas import BuildJobParams

    return BuildJobParams(
        recipe_name="fasta", genome_name=genome_name, asset_group_name=asset_group_name, **kw
    )


def make_slow_runner(gate, steps=3, fail=None):
    """A runner the test paces: it blocks on ``gate`` between steps.

    Progress goes through ``refgenie.progress.emit`` on purpose -- that is the
    path a real pull takes and the path cancellation interrupts.
    """
    from refgenie import progress

    def runner(ctx):
        ctx.phase("resolve", "starting")
        for i in range(steps):
            ctx.check_cancel()
            progress.emit("progress", current=i + 1, total=steps, unit="steps", phase="download")
            gate.wait(timeout=5)
        if fail is not None:
            raise fail
        return job_result_ok()

    return runner


def make_gated_runner(gate):
    """A runner that reaches RUNNING and then blocks until ``gate`` is set."""

    def runner(ctx):
        ctx.phase("download", "working")
        gate.wait(timeout=5)
        return job_result_ok()

    return runner


def open_gate():
    """A gate the runner never has to wait on."""
    gate = threading.Event()
    gate.set()
    return gate


def wait_for(predicate, timeout=5.0):
    """Spin until ``predicate`` holds. Fails the test rather than hanging."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.01)
    raise AssertionError("condition not reached within timeout")


def jobs_app(manager):
    """A bare app with only the jobs router included at ``/v1``.

    Bare on purpose: the contract is that the jobs router carries no dependency
    on the app factory.
    """
    from fastapi import FastAPI

    from refgenie.server.jobs import jobs_router

    app = FastAPI()
    app.state.job_manager = manager
    app.include_router(jobs_router, prefix="/v1")
    return app


def jobs_client(manager, **kwargs):
    """A (not yet entered) TestClient over ``jobs_app`` with the action header."""
    from fastapi.testclient import TestClient

    kwargs.setdefault("headers", ACTION_HEADERS)
    return TestClient(jobs_app(manager), **kwargs)


def submit(client, kind="pull", params=None, **kwargs):
    """POST one job and return the decoded ``JobRef``."""
    body = params if params is not None else (PULL_PARAMS if kind == "pull" else BUILD_PARAMS)
    return client.post("/v1/jobs", json={"kind": kind, "params": body}, **kwargs).json()


def submit_and_wait(client, manager, kind="pull", params=None):
    """POST one job, wait for it to finish, and return its id."""
    job_id = submit(client, kind, params)["job_id"]
    manager.wait(job_id, timeout=5)
    return job_id


def preflight(client, path, origin, method="GET", **extra_headers):
    """Issue a CORS preflight (``OPTIONS``) for ``origin``."""
    headers = {"Origin": origin, "Access-Control-Request-Method": method, **extra_headers}
    return client.options(path, headers=headers)


def read_sse(manager, until, headers=None, query="", timeout=10.0, during=None):
    """Read SSE frames off ``GET /v1/jobs/events`` until ``until(frame)``.

    Drives the ASGI app directly instead of using ``TestClient``, which buffers
    the whole body and can never read a stream that stays open -- the defining
    property of the endpoint under test.

    Returns ``(status, response headers, frames)``, each frame a dict of
    ``{"id", "event", "data"}``. ``during`` is an optional coroutine run
    concurrently, for asserting live events arrive on an open stream.
    """
    import asyncio

    app = jobs_app(manager)

    class _Stop(Exception):
        pass

    async def run():
        state = {"status": None, "headers": {}, "buffer": "", "frames": []}
        blocked = asyncio.Event()

        async def receive():
            await blocked.wait()
            return {"type": "http.disconnect"}

        async def send(message):
            if message["type"] == "http.response.start":
                state["status"] = message["status"]
                state["headers"] = {
                    key.decode().lower(): value.decode() for key, value in message["headers"]
                }
                return
            state["buffer"] += message.get("body", b"").decode()
            while "\n\n" in state["buffer"]:
                raw, state["buffer"] = state["buffer"].split("\n\n", 1)
                frame = parse_sse_frame(raw)
                state["frames"].append(frame)
                if until(frame):
                    raise _Stop

        scope = {
            "type": "http",
            "asgi": {"version": "3.0", "spec_version": "2.3"},
            "http_version": "1.1",
            "method": "GET",
            "scheme": "http",
            "path": "/v1/jobs/events",
            "raw_path": b"/v1/jobs/events",
            "root_path": "",
            "query_string": query.encode(),
            "headers": [(b"host", b"testserver")] + list(headers or []),
            "client": ("127.0.0.1", 12345),
            "server": ("testserver", 80),
        }

        tasks = [asyncio.create_task(app(scope, receive, send))]
        if during is not None:
            tasks.append(asyncio.create_task(during()))
        try:
            await asyncio.wait_for(asyncio.gather(*tasks), timeout)
        except _Stop:
            pass
        except TimeoutError:
            raise AssertionError(
                f"timed out; frames so far: {[f['event'] for f in state['frames']]}"
            ) from None
        finally:
            for task in tasks:
                task.cancel()
        return state["status"], state["headers"], state["frames"]

    return asyncio.run(run())


def parse_sse_frame(raw):
    """Parse one SSE frame. Comment lines are a contract violation, not data."""
    frame = {}
    for line in raw.split("\n"):
        line = line.rstrip("\r")
        if not line:
            continue
        assert not line.startswith(":"), f"comment frames are forbidden: {line!r}"
        field, _, value = line.partition(": ")
        frame[field] = json.loads(value) if field == "data" else value
    return frame


def until_done(frame):
    """The ``until`` predicate almost every SSE test uses."""
    return frame["event"] == "done"


# ---------------------------------------------------------------------------
# Job phase vocabulary
# ---------------------------------------------------------------------------

#: The phase names the web UI knows how to render, per job kind.
#:
#: Transcribed from ``frontend/src/components/jobs/phases.ts``, which turns each
#: into a label and a "step 4 of 8" counter, and gives three of them explicit
#: "this is silent and slow" copy. A phase the backend emits that is not here
#: still renders (humanized), but loses the step counter and the reassurance --
#: so the two lists are kept in sync deliberately, not incidentally.
PHASE_VOCABULARY = {
    "pull": {
        "resolve", "query", "stage_lookup", "download",
        "verify", "unpack", "register", "symlink",
    },
    "build": {
        "resolve_recipe", "seek_keys", "validate", "run",
        "digest", "register", "stage",
    },
    "genome_init": {"resolve", "digest", "register", "build"},
}
