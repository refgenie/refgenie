"""
Genome alias managers. Three backends satisfy AliasBackend and mode.py selects
between them: AliasManager (SQL), StoreAliasManager (RefgetStore), and
FederatedAliasManager (the union of the two, used by local mode).
"""

from pathlib import Path
from collections.abc import Callable, Iterable
from typing import Any, Protocol

from refget.store import RefgetStore
from rich.table import Table
from sqlalchemy.engine import Engine
from sqlmodel import select

from refgenie.db.tables import Alias, Genome
from refgenie.exceptions import FederatedAliasError, MissingAliasError, MissingGenomeError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.utils.symlinks import remove_alias_files
from refgenie.utils.tables import build_table

NAMESPACE = "refgenie"


class AliasBackend(Protocol):
    """Protocol satisfied by every alias backend in this module."""

    def resolve(self, name: str) -> str: ...
    def add(self, name: str, genome_digest: str) -> Any: ...
    def remove(self, name: str) -> None: ...
    def list_all(self, genome_digest: str | None = None) -> Iterable: ...
    def get_for_genome(self, genome_digest: str) -> list[str]: ...
    def exists(self, name: str) -> bool: ...
    def invalidate(self) -> None: ...
    def table(
        self,
        aliases: list[str] | None = None,
        genome_digests: list[str] | None = None,
    ) -> Table: ...


def _alias_table(
    by_digest: dict[str, tuple[list[str], str | None]],
    aliases: list[str] | None,
    genome_digests: list[str] | None,
    extra_column: str | None = None,
) -> Table:
    """
    Build the "Genome aliases" Rich table shared by both alias backends.

    Args:
        by_digest: Mapping of genome digest to (alias names, extra value).
        aliases: If given, keep only digests carrying at least one of these names.
        genome_digests: If given, keep only these digests.
        extra_column: Optional third column header; its values come from the
            second element of each ``by_digest`` tuple.

    Returns:
        Table: A Rich table showing aliases.
    """
    rows = []
    for digest, (names, extra) in by_digest.items():
        if aliases is not None and not any(a in names for a in aliases):
            continue
        if genome_digests is not None and digest not in genome_digests:
            continue
        row = [", ".join(names), digest]
        if extra_column is not None:
            row.append(extra)
        rows.append(row)

    columns = ["Aliases", "Genome digest"]
    if extra_column is not None:
        columns.append(extra_column)
    return build_table("Genome aliases", columns, rows)


class AliasManager(ResourceManager):
    """
    Manager for genome alias operations.

    Handles CRUD for aliases - mappings between human-readable names
    and genome digests.
    """

    def __init__(self, database_engine: Engine):
        super().__init__(database_engine)

    def add(self, name: str, genome_digest: str, store_name: str | None = None) -> Alias:
        """
        Add an alias for a genome.

        Args:
            name: The alias name.
            genome_digest: The genome digest to alias.
            store_name: Owning store of the aliased genome (federation). Recorded
                on the row so collisions and qualified lookups are answerable
                without re-reading the stores.

        Returns:
            Alias: The created alias.

        Raises:
            MissingGenomeError: If genome_digest doesn't exist.
        """
        with self._database_session as session:
            genome = one_or_raise(
                session,
                select(Genome).where(Genome.digest == genome_digest),
                MissingGenomeError(genome=genome_digest),
                unique=True,
            )

            existing = session.exec(select(Alias).where(Alias.name == name)).unique().one_or_none()

            if existing:
                if existing.genome_digest == genome_digest:
                    logger.info(f"Alias '{name}' already exists for genome {genome_digest}")
                    if store_name is not None and existing.store_name != store_name:
                        existing.store_name = store_name
                        session.commit()
                        session.refresh(existing)
                    return existing
                # Update to point to new genome
                existing.genome_digest = genome_digest
                if store_name is not None:
                    existing.store_name = store_name
                session.commit()
                session.refresh(existing)
                logger.info(f"Updated alias '{name}' to point to {genome_digest}")
                return existing

            alias = Alias(name=name, genome=genome, store_name=store_name)
            session.add(alias)
            session.commit()
            session.refresh(alias)
            logger.info(f"Added alias: {name}")
            return alias

    def federated_sync(
        self,
        name: str,
        genome_digest: str,
        store_name: str,
        priorities: dict[str, int],
    ) -> str:
        """Register a store's ``(alias, digest)`` under the collision policy.

        The bare name maps to exactly one genome (``Alias.name`` is a global PK).
        When two stores offer the same name for different digests, the
        higher-priority store (lower ``priority`` integer) wins the bare name;
        the loser is demoted to a qualified ``store::name`` row, reachable by
        that qualified name or by digest. Every real collision emits a WARNING
        naming both stores, the alias, and both digests.

        Returns one of ``"inserted"``, ``"noop"``, ``"repointed"``, ``"kept"``.
        """

        def _priority(sname: str | None) -> float:
            if sname is None:
                return float("inf")
            return priorities.get(sname, float("inf"))

        with self._database_session as session:
            existing = session.exec(select(Alias).where(Alias.name == name)).unique().one_or_none()

            # No existing row -> insert.
            if existing is None:
                session.add(Alias(name=name, genome_digest=genome_digest, store_name=store_name))
                session.commit()
                return "inserted"

            # Same digest -> natural dedup; keep the row, backfill the owner.
            if existing.genome_digest == genome_digest:
                if existing.store_name is None:
                    existing.store_name = store_name
                    session.commit()
                return "noop"

            # Different digest -> a real collision resolved by store priority.
            existing_owner = existing.store_name
            if _priority(store_name) < _priority(existing_owner):
                # New store wins the bare name. Demote the previous winner.
                logger.warning(
                    f"Alias collision on '{name}': store '{store_name}' "
                    f"({genome_digest}) outranks '{existing_owner}' "
                    f"({existing.genome_digest}); '{existing_owner}' is now only "
                    f"reachable as '{existing_owner}::{name}' or by digest."
                )
                if existing_owner is not None:
                    self._upsert_qualified(
                        session, existing_owner, name, existing.genome_digest
                    )
                existing.genome_digest = genome_digest
                existing.store_name = store_name
                session.commit()
                return "repointed"

            # New store loses; store it under a qualified name.
            logger.warning(
                f"Alias collision on '{name}': store '{existing_owner}' "
                f"({existing.genome_digest}) outranks '{store_name}' "
                f"({genome_digest}); '{store_name}' is only reachable as "
                f"'{store_name}::{name}' or by digest."
            )
            self._upsert_qualified(session, store_name, name, genome_digest)
            session.commit()
            return "kept"

    @staticmethod
    def _upsert_qualified(session, store_name: str, name: str, genome_digest: str) -> None:
        """Insert/update a literal ``store::name`` row for a collision loser."""
        qualified = f"{store_name}::{name}"
        row = session.exec(select(Alias).where(Alias.name == qualified)).unique().one_or_none()
        if row is None:
            session.add(
                Alias(name=qualified, genome_digest=genome_digest, store_name=store_name)
            )
        else:
            row.genome_digest = genome_digest
            row.store_name = store_name

    def remove(self, name: str) -> None:
        """
        Remove an alias.

        Args:
            name: The alias name to remove.

        Raises:
            MissingAliasError: If alias doesn't exist.
        """
        with self._database_session as session:
            alias = one_or_raise(
                session,
                select(Alias).where(Alias.name == name),
                MissingAliasError(name),
                unique=True,
            )
            session.delete(alias)
            session.commit()
            logger.info(f"Removed alias: {name}")

    def resolve(self, name: str) -> str:
        """
        Resolve an alias to a genome digest.

        A ``store::alias`` qualified name resolves within that store: first via a
        literal qualified row (collision losers are stored this way), else via
        the bare winner if it belongs to that store. A bare name resolves to its
        global winner.

        Args:
            name: The alias name, optionally ``store::alias`` qualified.

        Returns:
            str: The genome digest.

        Raises:
            MissingAliasError: If alias doesn't exist.
        """
        with self._database_session as session:
            if "::" in name:
                store_name, _, bare = name.partition("::")
                # 1. A literal qualified row (how collision losers are stored).
                literal = session.exec(
                    select(Alias.genome_digest).where(Alias.name == name)
                ).one_or_none()
                if literal is not None:
                    return literal
                # 2. The bare winner, but only if it belongs to that store.
                row = session.exec(
                    select(Alias).where(Alias.name == bare)
                ).unique().one_or_none()
                if row is not None and row.store_name == store_name:
                    return row.genome_digest
                raise MissingAliasError(name)

            statement = select(Genome.digest).join(Alias).where(Alias.name == name)
            return one_or_raise(session, statement, MissingAliasError(name), unique=True)

    def list_all(self, genome_digest: str | None = None) -> Iterable[Alias]:
        """
        List all aliases, optionally filtered by genome.

        Args:
            genome_digest: Optional genome digest to filter by.

        Returns:
            Iterable[Alias]: All matching aliases.
        """
        with self._database_session as session:
            statement = select(Alias)
            if genome_digest:
                statement = statement.join(Genome).where(Genome.digest == genome_digest)
            return session.exec(statement).unique().all()

    def get_for_genome(self, genome_digest: str) -> list[str]:
        """
        Get alias names for a specific genome.

        Args:
            genome_digest: The genome digest.

        Returns:
            list[str]: Alias names for the genome.
        """
        with self._database_session as session:
            genome = one_or_raise(
                session,
                select(Genome).where(Genome.digest == genome_digest),
                MissingGenomeError(genome=genome_digest),
                unique=True,
            )
            return [alias.name for alias in genome.aliases]

    def exists(self, name: str) -> bool:
        """
        Check if an alias exists.

        Args:
            name: The alias name.

        Returns:
            bool: Whether the alias exists.
        """
        statement = select(Alias).where(Alias.name == name)
        with self._database_session as session:
            result = session.exec(statement)
            return bool(result.first())

    def invalidate(self) -> None:
        """No-op. Every read here goes straight to the database."""

    def table(
        self,
        aliases: list[str] | None = None,
        genome_digests: list[str] | None = None,
    ) -> Table:
        """
        Get a Rich table of aliases.

        Args:
            aliases: List of aliases to filter by.
            genome_digests: List of genome digests to filter by.

        Returns:
            Table: A Rich table showing aliases.
        """
        with self._database_session as session:
            genomes = session.exec(select(Genome)).all()

            aliases_by_genome = {}
            for genome in genomes:
                genome_aliases = session.exec(
                    select(Alias.name).where(Alias.genome_digest == genome.digest)
                ).all()
                aliases_by_genome[genome.digest] = (
                    list(genome_aliases),
                    genome.description,
                )

            return _alias_table(
                aliases_by_genome, aliases, genome_digests, extra_column="Genome description"
            )


class StoreAliasManager(AliasBackend):
    """Alias manager backed by RefgetStore 'refgenie' namespace.

    Maintains a Python-side mapping alongside the store because the store's
    get_collection_metadata_by_alias requires the collection to exist in
    the store, but aliases may point to digests not loaded as collections.
    """

    def __init__(
        self,
        refget_store_getter: Callable[[], RefgetStore] | None = None,
        genome_folder_getter: Callable[[], Path] | None = None,
    ):
        self._get_store = refget_store_getter
        self._get_genome_folder = genome_folder_getter
        self._cache: dict[str, str] = {}  # name -> digest

    def set_store_getter(self, getter: Callable[[], RefgetStore]):
        self._get_store = getter

    def set_genome_folder_getter(self, getter: Callable[[], Path]):
        self._get_genome_folder = getter

    @property
    def _store(self) -> RefgetStore:
        if self._get_store is None:
            raise RuntimeError("Store getter not set on StoreAliasManager")
        return self._get_store()

    def invalidate(self) -> None:
        """Drop the cache so the next read re-reads the store.

        The cache is keyed to one RefgetStore instance. Whenever the caller
        swaps that instance out, the cached names describe a store this manager
        no longer reads.
        """
        self._cache.clear()

    def _sync_from_store(self):
        """Load aliases from the store into the cache."""
        if self._cache:
            return
        aliases = self._store.list_collection_aliases(NAMESPACE) or []
        for alias_name in aliases:
            meta = self._store.get_collection_metadata_by_alias(NAMESPACE, alias_name)
            if meta is not None:
                self._cache[alias_name] = meta.digest

    def resolve(self, name: str) -> str:
        if name in self._cache:
            return self._cache[name]
        # Try from store metadata (works when collection is loaded)
        meta = self._store.get_collection_metadata_by_alias(NAMESPACE, name)
        if meta is not None:
            self._cache[name] = meta.digest
            return meta.digest
        raise MissingAliasError(name)

    def add(self, name: str, genome_digest: str) -> "_AliasRecord":
        self._store.add_collection_alias(NAMESPACE, name, genome_digest)
        self._cache[name] = genome_digest
        logger.info(f"Added alias: {name}")
        return _AliasRecord(name=name, genome_digest=genome_digest)

    def remove(self, name: str) -> None:
        """
        Remove an alias from the store, then the trees it owned.

        The store is this backend's catalog, so it decides first and the files
        follow -- the same order the SQL backend uses (its ``Alias``
        before_delete handler queues file removal for after the commit).

        The file removal runs even when the store has no such alias. That case is
        exactly what an interrupted removal leaves behind, and the trees are
        derived from the alias name alone, so cleaning them is safe whether or
        not the alias was there. Removing an alias that never existed is still
        an error, but the caller gets a clean filesystem with it.
        """
        missing = False
        if name not in self._cache:
            # Check store too
            meta = self._store.get_collection_metadata_by_alias(NAMESPACE, name)
            if meta is None:
                missing = True

        if not missing:
            self._store.remove_collection_alias(NAMESPACE, name)
            self._cache.pop(name, None)

        self._remove_alias_files(name)

        if missing:
            raise MissingAliasError(name)
        logger.info(f"Removed alias: {name}")

    def _remove_alias_files(self, name: str) -> None:
        """
        Remove the on-disk trees keyed by this alias.

        The SQL-backed AliasManager gets this from the ``Alias`` before_delete
        event handler, which never fires here because this manager stores
        aliases in the refget store rather than the database. Without it, local
        mode (the default) orphans both ``alias/<name>/`` and ``builds/<name>/``
        on every alias removal.
        """
        if self._get_genome_folder is None:
            logger.debug(f"No genome folder getter set; leaving files for alias '{name}' in place")
            return
        remove_alias_files(genome_folder=self._get_genome_folder(), alias_name=name)

    def list_all(self, genome_digest: str | None = None) -> list:
        self._sync_from_store()
        results = []
        for name, digest in self._cache.items():
            if genome_digest and digest != genome_digest:
                continue
            results.append(_AliasRecord(name=name, genome_digest=digest))
        return results

    def get_for_genome(self, genome_digest: str) -> list[str]:
        return [a.name for a in self.list_all(genome_digest=genome_digest)]

    def exists(self, name: str) -> bool:
        if name in self._cache:
            return True
        meta = self._store.get_collection_metadata_by_alias(NAMESPACE, name)
        if meta is not None:
            self._cache[name] = meta.digest
            return True
        return False

    def table(
        self,
        aliases: list[str] | None = None,
        genome_digests: list[str] | None = None,
    ) -> Table:
        all_aliases = self.list_all()

        by_digest: dict[str, tuple[list[str], str | None]] = {}
        for a in all_aliases:
            by_digest.setdefault(a.genome_digest, ([], None))[0].append(a.name)

        return _alias_table(by_digest, aliases, genome_digests)


class FederatedAliasManager(AliasBackend):
    """Local store aliases, with the SQL alias table behind them.

    A build node holds two halves of one alias space. Names for the genomes it
    built live in its RefgetStore. Names for the genomes it federates over live
    in the SQL ``alias`` table, written by ``refgenie store sync``. Neither half
    alone answers "what is this genome called", so reads consult both.

    Reads are local first. A genome built on this node owns its name and a
    federated store must never shadow it, which mirrors the store priority rule.
    With no stores registered the SQL half is empty and every read behaves
    exactly as the store-backed manager alone would.

    Writes go to the local store only. This node is not the owner of a federated
    store's alias space and must not mutate it.
    """

    def __init__(self, store_backend: "StoreAliasManager", sql_backend: AliasManager):
        self._local = store_backend
        self._sql = sql_backend

    def set_store_getter(self, getter: Callable[[], RefgetStore]):
        self._local.set_store_getter(getter)

    def resolve(self, name: str) -> str:
        try:
            return self._local.resolve(name)
        except MissingAliasError:
            return self._sql.resolve(name)

    def add(self, name: str, genome_digest: str) -> "_AliasRecord":
        """Add an alias to the local store.

        Refuses a name the SQL half already holds for a different genome.
        Local names win every read, so allowing this would leave the federated
        genome silently unreachable by the only name it has.
        """
        try:
            federated_digest = self._sql.resolve(name)
        except MissingAliasError:
            federated_digest = None
        if federated_digest is not None and federated_digest != genome_digest:
            raise FederatedAliasError(name, federated_digest)
        return self._local.add(name, genome_digest)

    def remove(self, name: str) -> None:
        self._local.remove(name)

    def list_all(self, genome_digest: str | None = None) -> list:
        seen: set[str] = set()
        results = []
        for backend in (self._local, self._sql):
            for alias in backend.list_all(genome_digest=genome_digest):
                if alias.name in seen:
                    continue
                seen.add(alias.name)
                results.append(_AliasRecord(name=alias.name, genome_digest=alias.genome_digest))
        return results

    def get_for_genome(self, genome_digest: str) -> list[str]:
        return [a.name for a in self.list_all(genome_digest=genome_digest)]

    def exists(self, name: str) -> bool:
        return self._local.exists(name) or self._sql.exists(name)

    def invalidate(self) -> None:
        self._local.invalidate()

    def table(
        self,
        aliases: list[str] | None = None,
        genome_digests: list[str] | None = None,
    ) -> Table:
        by_digest: dict[str, tuple[list[str], str | None]] = {}
        for a in self.list_all():
            by_digest.setdefault(a.genome_digest, ([], None))[0].append(a.name)
        return _alias_table(by_digest, aliases, genome_digests)


class _AliasRecord:
    """Lightweight record matching the interface of the SQL Alias model."""

    __slots__ = ("name", "genome_digest")

    def __init__(self, name: str, genome_digest: str):
        self.name = name
        self.genome_digest = genome_digest
