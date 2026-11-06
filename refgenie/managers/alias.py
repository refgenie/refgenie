"""
Genome alias managers. Two backends satisfy AliasBackend and mode.py selects
between them: AliasManager (SQL), and FederatedAliasManager (the local
RefgetStore's aliases, with the SQL ones behind them).
"""

from pathlib import Path
from collections.abc import Callable, Iterable
from typing import Any, Protocol

from refget.store import RefgetStore
from rich.table import Table

from refgenie.db.engine import Engine
from refgenie.db.queries import alias as queries
from refgenie.db.queries import genome as genome_queries
from refgenie.db.queries.search import validate_search_fields
from refgenie.db.tables import Alias
from refgenie.exceptions import FederatedAliasError, MissingAliasError
from refgenie.logger import logger
from refgenie.models import GenomeAlias, GenomeDigest, Page
from refgenie.managers.base import ResourceManager
from refgenie.plugins.events import NULL_EVENTS, EventSink, update_scope
from refgenie.plugins.hooks import Change
from refgenie.utils.search import Search, name_matches_operator
from refgenie.utils.symlinks import remove_alias_files
from refgenie.utils.tables import build_table

NAMESPACE = "refgenie"


class AliasBackend(Protocol):
    """Protocol satisfied by every alias backend in this module."""

    def resolve(self, name: GenomeAlias) -> GenomeDigest: ...
    def add(self, name: GenomeAlias, genome_digest: GenomeDigest) -> Any: ...
    def remove(self, name: GenomeAlias) -> None: ...
    def list_all(self, genome_digest: GenomeDigest | None = None) -> Iterable: ...
    def get_for_genome(self, genome_digest: GenomeDigest) -> list[str]: ...
    def get_for_genomes(
        self, genome_digests: list[GenomeDigest]
    ) -> dict[GenomeDigest, list[str]]: ...
    def exists(self, name: GenomeAlias) -> bool: ...
    def page(
        self,
        *,
        name: str | None = None,
        genome_digest: GenomeDigest | None = None,
        search: Search | None = None,
        offset: int,
        limit: int,
    ) -> Page: ...
    def invalidate(self) -> None: ...
    def table(
        self,
        aliases: list[GenomeAlias] | None = None,
        genome_digests: list[GenomeDigest] | None = None,
    ) -> Table: ...


#: Fields an alias search may match.
ALIAS_SEARCH_FIELDS = ["name"]


def _page_aliases(
    backend: AliasBackend,
    *,
    name: str | None,
    genome_digest: GenomeDigest | None,
    search: Search | None,
    offset: int,
    limit: int,
) -> Page:
    """
    One page of a backend's aliases, by name, filtered in memory.

    Shared by every backend: aliases may live in the RefgetStore, where SQL
    cannot reach them, so filtering and slicing happen over ``list_all``.

    Raises:
        ValueError: If the search names a field other than ``name``.
    """
    records = sorted(backend.list_all(genome_digest=genome_digest), key=lambda r: r.name)
    if name:
        records = [r for r in records if r.name == name]
    if search is not None:
        validate_search_fields(search.fields, ALIAS_SEARCH_FIELDS)
        if term := search.term.strip():
            records = [r for r in records if name_matches_operator(r.name, term, search.operator)]
    return Page(
        items=records[offset : offset + limit], offset=offset, limit=limit, total=len(records)
    )


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

    def __init__(self, database_engine: Engine, events: EventSink | None = None):
        """
        Args:
            database_engine: The database engine.
            events: Where alias changes are recorded for plugins. Given only
                when this manager is the top-level alias backend; inside
                ``FederatedAliasManager`` it keeps ``NULL_EVENTS``.
        """
        super().__init__(database_engine)
        self._events = events or NULL_EVENTS

    @update_scope
    def add(
        self, name: GenomeAlias, genome_digest: GenomeDigest, store_name: str | None = None
    ) -> Alias:
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
            genome = genome_queries.get(session, genome_digest)
            existing = queries.find(session, name)

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
                self._events.record(Change(action="alias_added", genome=genome_digest, alias=name))
                return existing

            alias = Alias(name=name, genome=genome, store_name=store_name)
            session.add(alias)
            session.commit()
            session.refresh(alias)
            logger.info(f"Added alias: {name}")
        self._events.record(Change(action="alias_added", genome=genome_digest, alias=name))
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
            existing = queries.find(session, name)

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
                    self._upsert_qualified(session, existing_owner, name, existing.genome_digest)
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
        row = queries.find(session, qualified)
        if row is None:
            session.add(Alias(name=qualified, genome_digest=genome_digest, store_name=store_name))
        else:
            row.genome_digest = genome_digest
            row.store_name = store_name

    @update_scope
    def remove(self, name: GenomeAlias) -> None:
        """
        Remove an alias.

        Args:
            name: The alias name to remove.

        Raises:
            MissingAliasError: If alias doesn't exist.
        """
        with self._database_session as session:
            alias = queries.get(session, name)
            genome_digest = alias.genome_digest
            session.delete(alias)
            session.commit()
            logger.info(f"Removed alias: {name}")
        self._events.record(Change(action="alias_removed", genome=genome_digest, alias=name))

    def resolve(self, name: GenomeAlias) -> GenomeDigest:
        """
        Resolve an alias to a genome digest.

        A ``store::alias`` qualified name resolves within that store: first via a
        literal qualified row (collision losers are stored this way), else via
        the bare winner if it belongs to that store. A bare name resolves to its
        global winner.

        Args:
            name: The alias name, optionally ``store::alias`` qualified.

        Returns:
            GenomeDigest: The genome digest.

        Raises:
            MissingAliasError: If alias doesn't exist.
        """
        with self._database_session as session:
            if "::" in name:
                store_name, _, bare = name.partition("::")
                # 1. A literal qualified row (how collision losers are stored).
                literal = queries.find(session, name)
                if literal is not None:
                    return GenomeDigest(literal.genome_digest)
                # 2. The bare winner, but only if it belongs to that store.
                row = queries.find(session, bare)
                if row is not None and row.store_name == store_name:
                    return GenomeDigest(row.genome_digest)
                raise MissingAliasError(name)

            return GenomeDigest(queries.get(session, name).genome_digest)

    def list_all(self, genome_digest: GenomeDigest | None = None) -> Iterable[Alias]:
        """
        List all aliases, optionally filtered by genome.

        Args:
            genome_digest: Optional genome digest to filter by.

        Returns:
            Iterable[Alias]: All matching aliases.
        """
        with self._database_session as session:
            return queries.list_all(session, genome_digest)

    def get_for_genome(self, genome_digest: GenomeDigest) -> list[str]:
        """
        Get alias names for a specific genome.

        Args:
            genome_digest: The genome digest.

        Returns:
            list[str]: Alias names for the genome.
        """
        with self._database_session as session:
            genome_queries.get(session, genome_digest)
            return queries.names_by_genome(session, [genome_digest]).get(genome_digest, [])

    def get_for_genomes(self, genome_digests: list[GenomeDigest]) -> dict[GenomeDigest, list[str]]:
        """
        Get alias names for many genomes at once, in a single query.

        Unlike :meth:`get_for_genome`, an unknown digest maps to an empty list
        rather than raising.

        Args:
            genome_digests: The genome digests to look up.

        Returns:
            dict[GenomeDigest, list[str]]: Alias names keyed by genome digest,
            with an entry for every requested digest.
        """
        by_digest: dict[GenomeDigest, list[str]] = {digest: [] for digest in genome_digests}
        if not genome_digests:
            return by_digest
        with self._database_session as session:
            by_digest.update(queries.names_by_genome(session, list(genome_digests)))
        return by_digest

    def exists(self, name: GenomeAlias) -> bool:
        """
        Check if an alias exists.

        Args:
            name: The alias name.

        Returns:
            bool: Whether the alias exists.
        """
        with self._database_session as session:
            return queries.exists(session, name)

    def page(self, *, name=None, genome_digest=None, search=None, offset: int, limit: int) -> Page:
        """One page of aliases, by name. See :func:`_page_aliases`."""
        return _page_aliases(
            self, name=name, genome_digest=genome_digest, search=search, offset=offset, limit=limit
        )

    def invalidate(self) -> None:
        """No-op. Every read here goes straight to the database."""

    def table(
        self,
        aliases: list[GenomeAlias] | None = None,
        genome_digests: list[GenomeDigest] | None = None,
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
            genomes = genome_queries.list_all(session)
            names = queries.names_by_genome(session)
        aliases_by_genome = {
            genome.digest: (names.get(genome.digest, []), genome.description) for genome in genomes
        }
        return _alias_table(
            aliases_by_genome, aliases, genome_digests, extra_column="Genome description"
        )


class _StoreAliases:
    """The local RefgetStore half of :class:`FederatedAliasManager`.

    Keeps a name -> digest cache beside the store, because the store's
    get_collection_metadata_by_alias needs the collection loaded, and an alias
    may point to a digest that is not.
    """

    def __init__(
        self,
        refget_store_getter: Callable[[], RefgetStore],
        genome_folder_getter: Callable[[], Path],
    ):
        self._get_store = refget_store_getter
        self._get_genome_folder = genome_folder_getter
        self._cache: dict[str, str] = {}  # name -> digest

    @property
    def _store(self) -> RefgetStore:
        return self._get_store()

    def invalidate(self) -> None:
        """Drop the cache; call it whenever the RefgetStore instance is swapped out."""
        self._cache.clear()

    def _sync_from_store(self):
        """Load aliases from the store into the cache."""
        if self._cache:
            return
        for alias_name in self._store.list_collection_aliases(NAMESPACE) or []:
            meta = self._store.get_collection_metadata_by_alias(NAMESPACE, alias_name)
            if meta is not None:
                self._cache[alias_name] = meta.digest

    def resolve(self, name: GenomeAlias) -> GenomeDigest:
        if name in self._cache:
            return GenomeDigest(self._cache[name])
        meta = self._store.get_collection_metadata_by_alias(NAMESPACE, name)
        if meta is None:
            raise MissingAliasError(name)
        self._cache[name] = meta.digest
        return GenomeDigest(meta.digest)

    def add(self, name: GenomeAlias, genome_digest: GenomeDigest) -> "_AliasRecord":
        self._store.add_collection_alias(NAMESPACE, name, genome_digest)
        self._cache[name] = genome_digest
        logger.info(f"Added alias: {name}")
        return _AliasRecord(name=name, genome_digest=genome_digest)

    def remove(self, name: GenomeAlias) -> None:
        """
        Remove an alias from the store, then ``alias/<name>/`` and ``builds/<name>/``.

        The files are removed even when the alias is missing, which cleans up
        after an interrupted removal; MissingAliasError is still raised.
        """
        try:
            self.resolve(name)
            missing = False
        except MissingAliasError:
            missing = True
        if not missing:
            self._store.remove_collection_alias(NAMESPACE, name)
            self._cache.pop(name, None)
        remove_alias_files(genome_folder=self._get_genome_folder(), alias_name=name)
        if missing:
            raise MissingAliasError(name)
        logger.info(f"Removed alias: {name}")

    def list_all(self, genome_digest: GenomeDigest | None = None) -> list["_AliasRecord"]:
        self._sync_from_store()
        return [
            _AliasRecord(name=name, genome_digest=digest)
            for name, digest in self._cache.items()
            if not genome_digest or digest == genome_digest
        ]

    def get_for_genomes(self, genome_digests: list[GenomeDigest]) -> dict[GenomeDigest, list[str]]:
        """Alias names for many genomes, from one pass over the cached store."""
        by_digest: dict[GenomeDigest, list[str]] = {digest: [] for digest in genome_digests}
        if not genome_digests:
            return by_digest
        self._sync_from_store()
        for name, digest in self._cache.items():
            if digest in by_digest:
                by_digest[digest].append(name)
        return by_digest


class FederatedAliasManager(AliasBackend):
    """Local store aliases, with the SQL alias table behind them.

    A build node holds two halves of one alias space. Names for the genomes it
    built live in its RefgetStore. Names for the genomes it federates over live
    in the SQL ``alias`` table, written by ``refgenie store sync``. Neither half
    alone answers "what is this genome called", so reads consult both.

    Reads are local first, so a federated store never shadows a name built
    here. Writes go to the local store only; this node does not own a
    federated store's alias space.
    """

    def __init__(
        self,
        store_backend: _StoreAliases,
        sql_backend: AliasManager,
        events: EventSink | None = None,
    ):
        self._local = store_backend
        self._sql = sql_backend
        # The outermost backend records alias changes; the two halves do not,
        # so one change is recorded once.
        self._events = events or NULL_EVENTS

    def resolve(self, name: GenomeAlias) -> GenomeDigest:
        try:
            return self._local.resolve(name)
        except MissingAliasError:
            return self._sql.resolve(name)

    @update_scope
    def add(self, name: GenomeAlias, genome_digest: GenomeDigest) -> "_AliasRecord":
        """Add an alias to the local store.

        Refuses a name the SQL half holds for a different genome, since the
        local name would hide the federated genome.
        """
        try:
            federated_digest = self._sql.resolve(name)
        except MissingAliasError:
            federated_digest = None
        if federated_digest is not None and federated_digest != genome_digest:
            raise FederatedAliasError(name, federated_digest)
        try:
            unchanged = self._local.resolve(name) == genome_digest
        except MissingAliasError:
            unchanged = False
        record = self._local.add(name, genome_digest)
        if not unchanged:
            self._events.record(Change(action="alias_added", genome=genome_digest, alias=name))
        return record

    @update_scope
    def remove(self, name: GenomeAlias) -> None:
        try:
            genome_digest = self._local.resolve(name)
        except MissingAliasError:
            genome_digest = None
        self._local.remove(name)
        self._events.record(Change(action="alias_removed", genome=genome_digest, alias=name))

    def list_all(self, genome_digest: GenomeDigest | None = None) -> list:
        seen: set[str] = set()
        results = []
        for backend in (self._local, self._sql):
            for alias in backend.list_all(genome_digest=genome_digest):
                if alias.name in seen:
                    continue
                seen.add(alias.name)
                results.append(_AliasRecord(name=alias.name, genome_digest=alias.genome_digest))
        return results

    def get_for_genome(self, genome_digest: GenomeDigest) -> list[str]:
        return [a.name for a in self.list_all(genome_digest=genome_digest)]

    def get_for_genomes(self, genome_digests: list[GenomeDigest]) -> dict[GenomeDigest, list[str]]:
        """Alias names for many genomes, merged local-first like :meth:`list_all`."""
        by_digest: dict[GenomeDigest, list[str]] = {digest: [] for digest in genome_digests}
        if not genome_digests:
            return by_digest
        seen: set[str] = set()
        for backend in (self._local, self._sql):
            for digest, names in backend.get_for_genomes(genome_digests).items():
                for name in names:
                    if name in seen:
                        continue
                    seen.add(name)
                    by_digest[digest].append(name)
        return by_digest

    def exists(self, name: GenomeAlias) -> bool:
        try:
            self.resolve(name)
        except MissingAliasError:
            return False
        return True

    def page(self, *, name=None, genome_digest=None, search=None, offset: int, limit: int) -> Page:
        return _page_aliases(
            self, name=name, genome_digest=genome_digest, search=search, offset=offset, limit=limit
        )

    def invalidate(self) -> None:
        self._local.invalidate()

    def table(
        self,
        aliases: list[GenomeAlias] | None = None,
        genome_digests: list[GenomeDigest] | None = None,
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
