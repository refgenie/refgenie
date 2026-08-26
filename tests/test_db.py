"""
Tests for refgenie.db: ORM models/relationships and the alembic migration chain.

Manager-level write-path atomicity (add_from_path, write ordering) lives with
the managers tests, not here.

Also here: refgenie configuration (database config, migration state, the
legacy env-var warning, and server subscriptions), the server exception
classifier (refgenie.server.errors.classify_exception), and the
fixture-hygiene guards over the test tree itself.
"""

import ast
import importlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pytest
from alembic.autogenerate import compare_metadata
from alembic.config import Config
from alembic.migration import MigrationContext
from alembic.script import ScriptDirectory
from sqlalchemy import create_engine, inspect, text as sa_text
from sqlalchemy.engine.url import make_url
from sqlmodel import Session, SQLModel, select

import refgenie.config as config_module
import refgenie.db.tables  # noqa: F401  -- populates SQLModel.metadata
from refgenie import Refgenie
from refgenie.config.db import DatabaseType, create_default_db_config
from refgenie.const import TARGET_ALEMBIC_VERSION
from refgenie.db.migrations.utils import run_sql_migrations
from refgenie.utils.build import BUILD_DIGEST_SCHEME, build_level1_to_digest
from refgenie.exceptions import (
    AssetExistsError,
    MissingAliasError,
    MissingAssetError,
    MissingBuildInputError,
    MissingGenomeError,
    MissingRecipeError,
    NoArchiveError,
    PullFailedError,
    PullSkipped,
    RefgenieError,
)
from refgenie.server.errors import (
    _CLASSIFICATION,
    ErrorCode,
    WebError,
    action_error,
    classify_exception,
)


class TestDatabaseModels:
    """refgenie.db.tables: genome/asset persistence and ORM relationships."""

    def test_genome_round_trip_by_digest(self, refgenie_session):
        """alias.resolve yields a 32-char seqcol digest; genome.get round-trips it."""
        genome_digest = refgenie_session.alias.resolve("rCRSd")
        assert genome_digest is not None
        assert len(genome_digest) == 32  # seqcol digest length

        genome = refgenie_session.genome.get(genome_digest)
        assert genome is not None
        assert genome.digest == genome_digest

    def test_asset_relationship_chain(self, refgenie_session):
        """The asset → asset_group → genome ORM chain resolves both ways."""
        genome_digest = refgenie_session.alias.resolve("rCRSd")
        asset = refgenie_session.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name="test"
        )
        assert asset is not None
        assert asset.asset_group is not None
        assert asset.asset_group.name == "fasta"
        assert asset.asset_group.genome.digest == genome_digest


class TestManagerFacades:
    """Registry-facade smoke coverage.

    The only place configuration.list_all, stage.list_all, and
    sources.list_channels are exercised.
    """

    def test_facades_return_lists(self, refgenie_minimal):
        r = refgenie_minimal
        assert isinstance(r.configuration.list_all(), list)
        assert isinstance(r.stage.list_all(), list)
        assert isinstance(r.sources.list_channels(), list)


# ---------------------------------------------------------------------------
# Guards on the alembic migration chain.
#
# Both tests defend against the migration chain and the code that consumes it
# drifting apart silently. This has already happened once: TARGET_ALEMBIC_VERSION
# sat three revisions behind head, so every new database was stamped mid-chain and
# the last three migrations never ran on anything -- and nothing said so.
#
# If one of these fails after you add a migration, the migration is fine and the
# test is doing its job. Bump TARGET_ALEMBIC_VERSION in refgenie/const.py to the
# new head. Do not delete the test.
# ---------------------------------------------------------------------------

VERSION_TABLE = "alembic_version"  # alembic-managed; never part of a migration's DDL
MIGRATIONS_DIR = Path(refgenie.db.tables.__file__).parent / "migrations"


def _script_directory() -> ScriptDirectory:
    config = Config()
    config.set_main_option("script_location", str(MIGRATIONS_DIR))
    return ScriptDirectory.from_config(config)


def _downgrade(db_url, revision: str) -> None:
    """``run_sql_migrations`` in reverse; there is no shipped downgrade helper."""
    from alembic import command

    config = Config()
    config.set_main_option("script_location", str(MIGRATIONS_DIR))
    config.set_main_option(
        "sqlalchemy.url", db_url.render_as_string(hide_password=False).replace("%", "%%")
    )
    config.set_main_option("version_path_separator", "os")
    command.downgrade(config=config, revision=revision)


class TestAlembicChain:
    """The chain head and the schema it produces must stay aligned with the models."""

    def test_target_version_is_chain_head(self):
        """
        TARGET_ALEMBIC_VERSION must name the single head of the chain. If it names
        anything else, every revision after it is dead code: new databases are
        stamped to a version they are not and check_for_db_migrations() reports
        them up to date. Nothing raises; the schema is simply wrong.
        """
        heads = _script_directory().get_heads()
        assert len(heads) == 1, (
            f"expected a single migration head, found {heads}. "
            "A branched chain means someone added a revision without rebasing onto the head."
        )
        assert TARGET_ALEMBIC_VERSION == heads[0], (
            f"TARGET_ALEMBIC_VERSION is {TARGET_ALEMBIC_VERSION!r} but the chain head is "
            f"{heads[0]!r}. Update TARGET_ALEMBIC_VERSION in refgenie/const.py to the new head, "
            "or the migrations between the two will never run."
        )

    def test_upgrade_from_base_matches_models(self, tmp_path):
        """
        `alembic upgrade head` on an empty database must reproduce
        SQLModel.metadata — tables, columns, and (via compare_metadata) types,
        nullability, indexes, unique constraints and FKs. Also catches a chain
        that cannot apply from base at all. tmp_path is a throwaway.
        """
        # Scope the comparison to refgenie's own tables. SQLModel.metadata is a
        # process-global registry: gtars' refget-store models (sequence,
        # sequencecollection, namesattr, ...) register into it too once any test
        # opens a store, so comparing against the whole metadata would report
        # those foreign tables as drift depending on test order.
        model_table_names = {
            m.__tablename__
            for m in vars(refgenie.db.tables).values()
            if isinstance(m, type)
            and issubclass(m, SQLModel)
            and getattr(m, "__table__", None) is not None
        } - {VERSION_TABLE}

        db_path = tmp_path / "migrated.db"
        run_sql_migrations(db_url=make_url(f"sqlite:///{db_path}"))

        engine = create_engine(f"sqlite:///{db_path}")
        try:
            inspector = inspect(engine)

            db_tables = set(inspector.get_table_names()) - {VERSION_TABLE}
            assert model_table_names == db_tables, (
                f"tables only in the models: {sorted(model_table_names - db_tables)}; "
                f"tables only in the migrated database: {sorted(db_tables - model_table_names)}"
            )

            for table in sorted(model_table_names):
                model_columns = {c.name for c in SQLModel.metadata.tables[table].columns}
                db_columns = {c["name"] for c in inspector.get_columns(table)}
                assert model_columns == db_columns, (
                    f"{table}: columns only in the model {sorted(model_columns - db_columns)}; "
                    f"columns only in the migrated database {sorted(db_columns - model_columns)}"
                )

            # Beyond names: types, nullability, indexes, unique constraints, FKs.
            # include_object restricts the diff to refgenie's tables so foreign
            # SQLModel tables sharing the global metadata are ignored.
            def _include_object(obj, name, type_, reflected, compare_to):
                if type_ == "table":
                    return name in model_table_names or name == VERSION_TABLE
                return True

            with engine.connect() as conn:
                context = MigrationContext.configure(
                    conn, opts={"compare_type": True, "include_object": _include_object}
                )
                diffs = [
                    d
                    for d in compare_metadata(context, SQLModel.metadata)
                    if VERSION_TABLE not in repr(d)
                ]
            assert not diffs, "migrated schema differs from SQLModel.metadata:\n" + "\n".join(
                f"  {d}" for d in diffs
            )
        finally:
            engine.dispose()

    def test_check_for_db_migrations_does_not_leak_connections(self, tmp_path):
        """check_for_db_migrations must return its pooled connection every call;
        `MigrationContext.configure(engine.connect())` once never closed it."""
        engine = create_engine(f"sqlite:///{tmp_path / 'refgenie.db'}", echo=False)
        r = Refgenie(database_engine=engine, suppress_migrations=True)
        r.init(genome_folder=tmp_path / "genomes")

        baseline = engine.pool.checkedout()
        for _ in range(5):
            r.check_for_db_migrations()

        assert engine.pool.checkedout() == baseline


# ---------------------------------------------------------------------------
# Configuration: database config, migration state, the legacy env-var warning,
# and server subscriptions
#
# Kept to the real behavioral contracts: the default config, the instance's
# resolved config, the migration-needed check, the import-time warning fired
# when legacy ``$REFGENIE`` is set but ``$REFGENIE_DB_CONFIG_PATH`` is not,
# and subscribe/unsubscribe (pure database operations on the Configuration
# table). (Schema/migration mechanics are covered by TestAlembicChain above,
# including the migration-check connection-leak edge.)
# ---------------------------------------------------------------------------


class TestDatabaseConfig:
    """refgenie.config.db: default and instance database configuration."""

    def test_default_db_config_is_sqlite(self):
        """The default config is SQLite with a sqlite:///…/refgenie URL."""
        config = create_default_db_config()
        assert config.type == DatabaseType.SQLITE
        assert config.url.startswith("sqlite:///")
        assert config.url.endswith("/refgenie")

    def test_instance_database_config(self, refgenie_minimal):
        """get_database_config resolves to a SQLite config with a sqlite:/// URL."""
        config = refgenie_minimal.get_database_config()
        assert config.type == DatabaseType.SQLITE
        assert config.url.startswith("sqlite:///")


class TestMigrationCheck:
    """The catalog reports whether it needs a schema upgrade."""

    def test_freshly_initialized_catalog_needs_no_migration(self, refgenie_minimal):
        """A freshly initialized, head-stamped catalog reports no migration needed."""
        assert refgenie_minimal.check_for_db_migrations(log=False) is False


def test_warning_fires_when_REFGENIE_set_and_new_var_missing(monkeypatch, capsys):
    monkeypatch.setenv("REFGENIE", "/legacy/genome_config.yaml")
    monkeypatch.delenv("REFGENIE_DB_CONFIG_PATH", raising=False)
    importlib.reload(config_module)
    captured = capsys.readouterr()
    assert "REFGENIE_DB_CONFIG_PATH" in captured.err
    assert "/legacy/genome_config.yaml" in captured.err


@pytest.mark.parametrize(
    "legacy_set, new_var_set",
    [
        pytest.param(True, True, id="new_var_set"),
        pytest.param(False, False, id="REFGENIE_unset"),
    ],
)
def test_warning_silent(monkeypatch, capsys, tmp_path, legacy_set, new_var_set):
    """No warning when the new var is set, or when legacy $REFGENIE is unset."""
    if legacy_set:
        monkeypatch.setenv("REFGENIE", "/legacy/genome_config.yaml")
    else:
        monkeypatch.delenv("REFGENIE", raising=False)
    if new_var_set:
        monkeypatch.setenv("REFGENIE_DB_CONFIG_PATH", str(tmp_path / "db.yaml"))
    else:
        monkeypatch.delenv("REFGENIE_DB_CONFIG_PATH", raising=False)
    importlib.reload(config_module)
    captured = capsys.readouterr()
    assert "WARNING" not in captured.err


S1 = "http://server1.example.com"
S2 = "http://server2.example.com"
OLD = "http://old-server.example.com"
NEW = "http://new-server.example.com"


class TestSubscribe:
    """Test server subscription management."""

    @pytest.mark.parametrize(
        "ops, expected",
        [
            ([("subscribe", S1, {})], [S1]),
            ([("subscribe", [S1, S2], {})], [S1, S2]),
            ([("subscribe", S1, {}), ("subscribe", S1, {})], [S1]),
            ([("subscribe", OLD, {}), ("subscribe", NEW, {"reset": True})], [NEW]),
            ([("subscribe", [S1, S2], {}), ("unsubscribe", [S1], {})], [S2]),
            ([("subscribe", S1, {}), ("unsubscribe", ["http://nope.example.com"], {})], [S1]),
            ([], []),
        ],
        ids=[
            "single", "multiple", "deduplicates", "reset", "unsubscribe-one",
            "unsubscribe-unknown-is-noop", "none-subscribed",
        ],
    )
    def test_subscription_crud(self, refgenie_minimal, ops, expected):
        """Each row applies its operations in order, then pins the full subscription set."""
        for method, arg, kwargs in ops:
            getattr(refgenie_minimal.configuration, method)(arg, **kwargs)
        subs = list(refgenie_minimal.configuration.get_server_subscriptions())
        assert sorted(subs) == sorted(expected)

    def test_subscribe_when_configuration_id_is_not_one(self, refgenie_minimal):
        """The sole Configuration row need not have id 1 (restored dump, merged db)."""
        from refgenie.db.tables import Configuration

        with Session(refgenie_minimal.database_engine) as session:
            row = session.exec(select(Configuration)).one()
            session.delete(row)
            session.commit()
            session.add(
                Configuration(
                    id=42,
                    version=row.version,
                    genome_folder=row.genome_folder,
                    genome_stage_folder=row.genome_stage_folder,
                )
            )
            session.commit()

        refgenie_minimal.configuration.subscribe("http://server1.example.com")
        assert "http://server1.example.com" in refgenie_minimal.configuration.get_latest().servers
        refgenie_minimal.configuration.unsubscribe(["http://server1.example.com"])
        assert list(refgenie_minimal.configuration.get_server_subscriptions()) == []


# ---------------------------------------------------------------------------
# The one exception classifier: `refgenie.server.errors.classify_exception`
#
# The trap this section pins: `AssetExistsError` and `NoArchiveError` both
# subclass `PullFailedError`, so the classifier's isinstance scan must be
# ordered most specific first. A type-keyed dict or a mis-ordered `except`
# chain silently collapses the three into one code, destroying exactly the
# differentiation the pull UI exists to show.
# ---------------------------------------------------------------------------


class TestClassificationOrdering:
    """Most specific first, or the differentiation silently dies."""

    def test_asset_exists_is_not_collapsed_into_pull_failed(self):
        status, code, _ = classify_exception(AssetExistsError("asset already there"))
        assert (status, code) == (409, ErrorCode.ASSET_EXISTS)

    def test_no_archive_is_not_collapsed_into_pull_failed(self):
        status, code, _ = classify_exception(NoArchiveError("no archive on server"))
        assert (status, code) == (404, ErrorCode.NO_ARCHIVE)

    def test_plain_pull_failure_is_pull_failed(self):
        status, code, _ = classify_exception(PullFailedError("download died"))
        assert (status, code) == (502, ErrorCode.PULL_FAILED)

    def test_no_subclass_is_listed_after_its_parent(self):
        """The structural invariant behind the three tests above: a subclass
        appearing after its parent in the table would never be reached."""
        types = [entry[0] for entry in _CLASSIFICATION]
        for i, earlier in enumerate(types):
            for later in types[i + 1 :]:
                if later is earlier:
                    continue
                assert not issubclass(later, earlier), (
                    f"{later.__name__} is listed after its parent {earlier.__name__}; "
                    "it would never match"
                )


class TestSpecificCodes:
    """Each lookup failure keeps its own code despite subclassing RefgenieError."""

    @pytest.mark.parametrize(
        "exc,status,code",
        [
            (MissingGenomeError(genome="x"), 404, ErrorCode.GENOME_NOT_FOUND),
            (MissingAliasError("hg38"), 404, ErrorCode.ALIAS_NOT_FOUND),
            (MissingAssetError(digest="abc"), 404, ErrorCode.ASSET_NOT_FOUND),
            (MissingRecipeError("fasta"), 404, ErrorCode.RECIPE_NOT_FOUND),
            (MissingBuildInputError("need a file"), 400, ErrorCode.MISSING_BUILD_INPUT),
            (PullSkipped("declined"), 409, ErrorCode.PULL_SKIPPED),
        ],
    )
    def test_lookup_and_pull_codes(self, exc, status, code):
        got_status, got_code, _ = classify_exception(exc)
        assert (got_status, got_code) == (status, code)

    def test_bare_refgenie_error_is_the_500_catch_all(self):
        status, code, _ = classify_exception(RefgenieError("something domain-y"))
        assert (status, code) == (500, ErrorCode.REFGENIE_ERROR)

    def test_value_error_is_a_conflict(self):
        status, code, _ = classify_exception(ValueError("asset has children"))
        assert (status, code) == (409, ErrorCode.CONFLICT)

    def test_unknown_exception_is_internal_error(self):
        status, code, _ = classify_exception(KeyError("nope"))
        assert (status, code) == (500, ErrorCode.INTERNAL_ERROR)


class TestWebError:
    """Synthesized failures carry their own code through the same classifier."""

    def test_web_error_passes_its_code_and_status_through(self):
        exc = WebError("no servers", ErrorCode.NO_SUBSCRIPTIONS, 409)
        assert classify_exception(exc) == (409, ErrorCode.NO_SUBSCRIPTIONS, "no servers")

    def test_build_failed_is_synthesizable(self):
        exc = WebError("pipeline failed", ErrorCode.BUILD_FAILED, 500)
        status, code, _ = classify_exception(exc)
        assert (status, code) == (500, ErrorCode.BUILD_FAILED)


class TestActionError:
    """`action_error` wraps the classification into an envelope HTTPException."""

    def test_carries_classified_status_and_envelope_detail(self):
        exc = action_error(MissingGenomeError(genome="abc123"))
        assert exc.status_code == 404
        assert exc.detail["code"] == str(ErrorCode.GENOME_NOT_FOUND)
        assert "abc123" in exc.detail["message"]


# ---------------------------------------------------------------------------
# Fixture hygiene: guard against fixture shadowing and helpers leaking into
# conftest
#
# Two shadowed fixtures have already cost real debugging time: a module-local
# ``refgenie_built`` that silently changed the asset name for ten tests, and a
# duplicate ``fixtures_path`` in tests/integration/conftest.py. Both are the
# same mistake -- redefining a name the root conftest already owns -- and both
# are statically detectable.
# ---------------------------------------------------------------------------

TESTS_ROOT = Path(__file__).parent

#: Names allowed to shadow a root-conftest fixture, with the reason.
ALLOWED_SHADOWS: dict[str, str] = {}


def _fixture_names(path: Path) -> set[str]:
    """Every ``@pytest.fixture``-decorated function name defined in ``path``."""
    tree = ast.parse(path.read_text(), filename=str(path))
    names = set()
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for decorator in node.decorator_list:
            target = decorator.func if isinstance(decorator, ast.Call) else decorator
            attr = getattr(target, "attr", None) or getattr(target, "id", None)
            if attr == "fixture":
                names.add(node.name)
    return names


def _test_modules() -> list[Path]:
    return sorted(p for p in TESTS_ROOT.rglob("*.py") if p != TESTS_ROOT / "conftest.py")


def test_no_fixture_shadows_the_root_conftest():
    """No module may redefine a fixture the root conftest already provides."""
    root_fixtures = _fixture_names(TESTS_ROOT / "conftest.py")
    shadows = {}
    for module in _test_modules():
        collisions = (_fixture_names(module) & root_fixtures) - set(ALLOWED_SHADOWS)
        if collisions:
            shadows[str(module.relative_to(TESTS_ROOT))] = sorted(collisions)
    assert not shadows, (
        f"These modules redefine root-conftest fixtures: {shadows}. Rename the "
        "local fixture, or add it to ALLOWED_SHADOWS with a reason."
    )


def test_conftest_defines_only_fixtures_and_hooks():
    """tests/conftest.py holds fixtures and pytest hooks; helpers live in helpers.py."""
    tree = ast.parse((TESTS_ROOT / "conftest.py").read_text())
    fixtures = _fixture_names(TESTS_ROOT / "conftest.py")
    offenders = [
        node.name
        for node in tree.body
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef))
        and node.name not in fixtures
        and not node.name.startswith(("pytest_", "_"))
    ]
    assert not offenders, (
        f"tests/conftest.py defines plain helper function(s) {offenders}; "
        "move them to tests/helpers.py."
    )
