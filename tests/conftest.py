"""
Refgenie Unit Test Fixtures

Fixtures and pytest hooks only. Every plain helper function, class and constant
lives in ``tests/helpers.py``; this module imports what its fixtures need.

Lightweight fixtures for fast unit tests (<15s total). Tests that build real
genome folders, asset files or archives on disk are the ``component`` tier (see
``tests/README.md``); ``tests/integration/`` is the tier that needs live
services (Docker, PostgreSQL, bulker).
"""

import atexit
import os
import shutil
import tempfile
from pathlib import Path

import pytest

# Environment variables that can silently redirect a test run at real data.
#
# REFGENIE_DB_CONFIG_PATH is the dangerous one: refgenie resolves genome_folder
# and genome_stage_folder from the catalog's `configuration` row, NOT from the
# environment, so pointing this at a production catalog sends every write --
# genome rows, collections, aliases, .seq files -- into that catalog's store,
# wherever it lives. Nothing in the test suite overrides it, because the suite
# has no reason to expect it to be set.
#
# The guard below refuses to run when these are inherited from the environment.
_AMBIENT_ENV_VARS = (
    "REFGENIE_DB_CONFIG_PATH",
    "REFGENIE_GENOME_FOLDER",
    "REFGENIE_GENOME_STAGE_FOLDER",
    "REFGENIE_HOME_PATH",
    "REFGENIE_REFGET_STORE_URL",
)
_ALLOW_AMBIENT = "REFGENIE_ALLOW_AMBIENT_ENV"

# Snapshot the ambient environment BEFORE the redirect below injects a value of
# its own; pytest_configure inspects this snapshot, not the live environment.
_AMBIENT_SNAPSHOT = {var: os.environ.get(var) for var in _AMBIENT_ENV_VARS}

# Point refgenie's home at a disposable directory before anything imports
# refgenie. `refgenie.config` reads REFGENIE_HOME_PATH at import time and mkdirs
# it, and `Refgenie.init()` with no genome_folder falls back to
# config.genome_folder -- so without this, merely importing refgenie creates
# ~/.refgenie, and any bare init() writes fake genomes into the user's real
# data. setdefault keeps a deliberate override; if that override points at real
# data, pytest_configure below still aborts the run.
_TEST_HOME = Path(tempfile.mkdtemp(prefix="refgenie-test-home-"))
os.environ.setdefault("REFGENIE_HOME_PATH", str(_TEST_HOME))
if os.environ["REFGENIE_HOME_PATH"] == str(_TEST_HOME):
    atexit.register(lambda: shutil.rmtree(_TEST_HOME, ignore_errors=True))
else:
    shutil.rmtree(_TEST_HOME, ignore_errors=True)

from refgenie import Refgenie  # noqa: E402  (must follow the home redirect above)

from tests.helpers import (  # noqa: E402  (must follow the home redirect above)
    TESTS_DATA_DIR,
    build_rcrsd,
    make_engine,
    register_fasta,
    stage_rcrsd,
    web_assets_present,
)


def _is_disposable(value: str) -> bool:
    """Whether a path is somewhere a test may legitimately be pointed."""
    if value.startswith(("http://", "https://")):
        return False  # a remote store URL is never disposable
    try:
        resolved = Path(value).expanduser().resolve()
    except (OSError, RuntimeError):
        return False
    roots = [Path(tempfile.gettempdir()).resolve(), Path("/tmp"), Path("/scratch")]
    return any(resolved == r or r in resolved.parents for r in roots)


def _assert_default_home_is_not_real():
    """Fail loudly if refgenie's *default* paths still resolve under $HOME.

    The redirect at the top of this module only works if it runs before
    ``refgenie.config`` is imported. An import-ordering regression would silently
    restore the old behavior, so verify the outcome rather than the intent.
    """
    from refgenie.config import config as _refgenie_config

    genome_folder = Path(_refgenie_config.genome_folder)
    if Path.home() in genome_folder.parents:
        raise pytest.UsageError(
            f"refgenie's default genome_folder still resolves under $HOME "
            f"({genome_folder}); tests/conftest.py must set REFGENIE_HOME_PATH "
            "before importing refgenie."
        )


def pytest_configure(config):
    """Refuse to run against a catalog or store inherited from the environment.

    Set REFGENIE_ALLOW_AMBIENT_ENV=1 to override, but read the message first --
    the override exists for deliberate integration runs, not for silencing this.
    """
    if os.environ.get(_ALLOW_AMBIENT) == "1":
        return
    offenders = {
        var: value
        for var, value in _AMBIENT_SNAPSHOT.items()
        if value and not _is_disposable(value)
    }
    if not offenders:
        _assert_default_home_is_not_real()
        return
    lines = [
        "",
        "refgenie tests refuse to run: the environment points at real refgenie data.",
        "",
    ]
    lines += [f"    {var}={value}" for var, value in sorted(offenders.items())]
    lines += [
        "",
        "These are almost certainly inherited -- a sourced environment file or a",
        "shell/daemon that pre-exports them. refgenie reads genome_folder from the",
        "CATALOG, not the environment, so a stray REFGENIE_DB_CONFIG_PATH sends every",
        "test write into that catalog's genome store: fake genome rows, collections,",
        "aliases and .seq files land in production.",
        "",
        "Run the tests in a clean environment instead:",
        "",
        f"    env {' '.join('-u ' + v for v in sorted(offenders))} pytest",
        "",
        f"If you really mean to run against these paths, set {_ALLOW_AMBIENT}=1.",
        "",
    ]
    raise pytest.UsageError("\n".join(lines))


def pytest_sessionstart(session):
    """CI's anti-rot guard: a missing web UI bundle must fail loudly, not skip.

    ``tests/test_web.py`` marks most of its SPA tests ``requires_web_assets``
    because a fresh checkout has no built frontend. Left as a plain skip, that
    is exactly the failure mode this exists to prevent -- the old dash's browse
    path was a hard 500 for months and a green, silently-skipping suite is why
    nobody noticed. CI sets REFGENIE_REQUIRE_WEB_ASSETS=1 in every job that
    runs pytest (see .github/workflows/build-frontend.yml and its consumers);
    locally it is unset, so a skip stays a convenience.
    """
    if os.environ.get("REFGENIE_REQUIRE_WEB_ASSETS") and not web_assets_present():
        pytest.exit(
            "REFGENIE_REQUIRE_WEB_ASSETS is set but the web UI bundle "
            "(refgenie/server/webui/index.html) is missing. Build it first:\n"
            "    npm --prefix frontend ci && npm --prefix frontend run build",
            returncode=1,
        )


@pytest.fixture(autouse=True)
def _isolated_default_genome_folder(tmp_path, monkeypatch):
    """Any bare ``Refgenie.init()`` defaults into this test's tmp_path.

    ``Refgenie.__init__`` -> ``check_for_db_migrations`` can call a bare
    ``init()`` on its own (refgenie/core/lifecycle.py), so the default must be
    safe, not fatal. Tests should still pass an explicit ``genome_folder``; this
    is a net, not a license.
    """
    from refgenie.config import config

    monkeypatch.setattr(config, "genome_folder", tmp_path / "default_genomes")
    monkeypatch.setattr(config, "genome_stage_folder", tmp_path / "default_archives")


@pytest.fixture
def engine():
    """Function-scoped in-memory SQLite engine."""
    return make_engine()


@pytest.fixture(scope="session")
def fixtures_path():
    """Path to test data fixtures."""
    return TESTS_DATA_DIR


@pytest.fixture
def rCRSd_fasta_file_path(fixtures_path):
    """Path to rCRSd FASTA test file."""
    return fixtures_path / "rCRSd.fa"


@pytest.fixture
def bowtie2_index_asset_class_file_path(fixtures_path):
    """Path to bowtie2 index asset class definition."""
    return fixtures_path / "bowtie2_index_asset_class.yaml"


@pytest.fixture
def bowtie2_index_recipe_file_path(fixtures_path):
    """Path to bowtie2 index recipe definition."""
    return fixtures_path / "bowtie2_index_asset_recipe.yaml"


@pytest.fixture
def refgenie_minimal(engine, tmp_path):
    """
    Refgenie initialized without any assets built.
    Use for lightweight tests that don't need pre-built assets.

    Inits into tmp_path so tests never write to the real
    ``config.genome_folder`` default.
    """
    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.init(
        genome_folder=tmp_path / "genomes",
        genome_stage_folder=tmp_path / "archives",
    )
    return r


@pytest.fixture
def refgenie_with_fasta(engine, tmp_path, fixtures_path):
    """
    Refgenie initialized with fasta asset class and recipe registered
    (but no assets built). Use for tests that need fasta definitions.

    Inits into tmp_path so unit tests never write to the real
    ``config.genome_folder`` default.
    """
    r = Refgenie(database_engine=engine, suppress_migrations=True)
    r.init(
        genome_folder=tmp_path / "genomes",
        genome_stage_folder=tmp_path / "archives",
    )
    register_fasta(r, fixtures_path)
    return r


@pytest.fixture
def refgenie_with_genome(fixtures_path, tmp_path):
    """Refgenie with a genome initialized from FASTA (no built assets).

    Unlike ``refgenie_fs`` this one registers no fasta asset class or recipe,
    so nothing can be built against it.
    """
    r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    r.init(genome_folder=tmp_path / "genomes")
    r.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=["rCRSd"],
        description="rCRSd genome",
    )
    return r


@pytest.fixture
def refgenie_fs(refgenie_with_fasta, fixtures_path):
    """Refgenie on a real filesystem with the rCRSd genome initialized, nothing built."""
    refgenie_with_fasta.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=["rCRSd"],
        description="rCRSd genome",
    )
    return refgenie_with_fasta


@pytest.fixture
def refgenie_built(refgenie_fs):
    """Refgenie on a real filesystem with one fasta asset built as ``test``."""
    build_rcrsd(refgenie_fs)
    return refgenie_fs


@pytest.fixture
def staged_refgenie(refgenie_built):
    """
    Function-scoped refgenie with a built FASTA asset and staged records.

    Uses tmp_path for real filesystem paths so that StageManager.create() can
    write the .tgz tarball to disk. Function-scoped because unstage
    tests are destructive and each test needs a fresh staged asset.
    """
    stage_rcrsd(refgenie_built)
    return refgenie_built


@pytest.fixture(scope="session")
def refgenie_session(fixtures_path, tmp_path_factory):
    """
    Session-scoped refgenie with built FASTA asset.
    Use for read-only tests to avoid per-test setup overhead.

    The fasta recipe uses refgenie-build-fasta (pure Python, no samtools needed).

    Uses an isolated genome folder (not the default ``~/.refgenie``): content is
    digest-addressed, so a persistent genome folder paired with this in-memory
    DB would leave a build flag pointing at content the fresh DB doesn't know,
    and the skipped rebuild could not find its (relocated) output.
    """
    rCRSd_fasta = fixtures_path / "rCRSd.fa"
    genomes = tmp_path_factory.mktemp("refgenie_session_genomes")
    r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    r.init(genome_folder=genomes)
    register_fasta(r, fixtures_path)
    r.genome.initialize_genome(
        fasta_file_path=rCRSd_fasta,
        alias_names=["rCRSd"],
        description="rCRSd genome",
    )
    build_rcrsd(r)
    return r


@pytest.fixture(autouse=True)
def reset_source_cache():
    """Clear the remote-source cache around every test.

    The cache is a module-level dict keyed by URL, so a source built by one
    test is handed to the next one. Clearing it is one dict operation, so it
    is cheaper to do globally than to remember which classes need it.
    """
    from refgenie.managers.sources.genomes import clear_source_cache

    clear_source_cache()
    yield
    clear_source_cache()


@pytest.fixture
def server_client_world(engine, tmp_path, fixtures_path):
    """A served server + subscribed client pair: ``(client_rg, server_rg, url)``.

    Yields inside the live ``serve_refgenie`` context, so remote calls made by
    the test hit the real server routers over an in-process ASGI transport.
    """
    from tests.helpers import (
        MOCK_MIRROR_URL,
        MOCK_SERVER_URL,
        make_server_client_world,
        serve_refgenie,
    )

    server_rg, client_rg = make_server_client_world(engine, tmp_path, fixtures_path)
    # Both URLs are registered against the one server so multi-server tests
    # (which pass force_server_urls explicitly) need no separate fixture.
    with serve_refgenie(client_rg, server_rg, MOCK_SERVER_URL, MOCK_MIRROR_URL):
        yield client_rg, server_rg, MOCK_SERVER_URL


@pytest.fixture
def fake_subprocess():
    """``subprocess.run`` patched to succeed silently; yields the mock."""
    from unittest.mock import MagicMock, patch

    with patch("subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(stdout="", stderr="", returncode=0)
        yield mock_run


@pytest.fixture
def failing_subprocess():
    """Factory: ``failing_subprocess(exc)`` patches ``subprocess.run`` to raise."""
    from unittest.mock import patch

    with patch("subprocess.run") as mock_run:

        def configure(exc):
            mock_run.side_effect = exc
            return mock_run

        yield configure


@pytest.fixture
def job_manager_factory():
    """Build ``JobManager``s against fake runners and guarantee shutdown.

    A manager that outlives its test leaves two thread pools and a handler on
    the ``refgenie`` logger behind. Covers both the bare-manager tests (which
    pass ``runners=`` explicitly) and the router tests (default all-instant
    runners over a ``FakeRefgenie``).
    """
    from refgenie.server.jobs.manager import JobManager
    from refgenie.server.jobs.schemas import JobKind

    from tests.helpers import FakeRefgenie, instant_runner

    made = []

    def make(runners=None, refgenie=None, **kwargs):
        manager = JobManager(
            refgenie if refgenie is not None else FakeRefgenie(),
            runners=runners if runners is not None else {k: instant_runner for k in JobKind},
            **kwargs,
        )
        made.append(manager)
        return manager

    yield make
    for manager in made:
        manager.shutdown(wait=False)


# ---------------------------------------------------------------------------
# Test tiers
# ---------------------------------------------------------------------------
#
# Every test carries exactly one tier marker. Rather than decorating ~1000
# tests by hand, the tier is derived from location and can be overridden by
# marking a module explicitly (``pytestmark = pytest.mark.component``).
#
#   tests/e2e/**          -> e2e          (subprocess CLI, external binaries)
#   tests/integration/**  -> integration  (Docker services)
#   everything else       -> unit, unless the module says otherwise
#
# `addopts` in pyproject.toml deselects component and e2e so the bare `pytest`
# inner loop stays fast.

TIER_MARKERS = ("unit", "component", "e2e", "integration")

_TIER_BY_SUBDIR = {
    "e2e": "e2e",
    "integration": "integration",
}


def pytest_collection_modifyitems(config, items):
    """Assign a tier marker to every test that does not already declare one."""
    tests_root = Path(__file__).parent
    for item in items:
        if any(item.get_closest_marker(name) for name in TIER_MARKERS):
            continue
        try:
            relative = Path(str(item.path)).relative_to(tests_root)
        except ValueError:
            relative = None
        subdir = relative.parts[0] if relative and len(relative.parts) > 1 else ""
        item.add_marker(_TIER_BY_SUBDIR.get(subdir, "unit"))
