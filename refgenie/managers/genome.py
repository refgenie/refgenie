"""
GenomeManager - handles genome CRUD operations.
"""

import tempfile
from dataclasses import dataclass
from pathlib import Path
from collections.abc import Callable, Iterable
from typing import TYPE_CHECKING, Any, Protocol

from rich.table import Table

from refgenie.db.engine import Engine
from refgenie.db.queries import genome as queries
from refgenie.db.tables import Genome
from refgenie.exceptions import FastaDigestMismatchError, MissingAliasError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.sequence import open_remote_store
from refgenie.models import GenomeAlias, GenomeDetail, GenomeDigest, GenomeListing, Page
from refgenie.utils.search import Search, alias_digests_matching
from refgenie.utils.species_file import SpeciesEntry
from refgenie.utils.tables import build_table
from refgenie.managers.sources import RemoteGenomeSource
from refgenie.plugins.events import NULL_EVENTS, EventSink, update_scope
from refgenie.plugins.hooks import Change

if TYPE_CHECKING:
    from gtars.refget import RefgetStore
    from refgenie.managers.alias import AliasBackend
    from refgenie.managers.sources.servers import ServerManager
    from refgenie.managers.species import SpeciesManager


#: Genome columns a pull copies from the server's record when unset locally.
SERVER_GENOME_COLUMNS = (
    "description",
    "species_name",
    "synonyms",
    "taxon_id",
    "assembly_source",
    "assembly_accession",
    "assembly_level",
)

#: The placeholder description a genome gets when none was known at creation.
NO_DESCRIPTION = "No description provided"


def _is_unset(value, column: str) -> bool:
    """Whether a local genome column holds no user-set value."""
    if value is None or value == "" or value == []:
        return True
    return column == "description" and value == NO_DESCRIPTION


def _clean_synonyms(names: list) -> list[str] | None:
    """FHR ``genomeSynonym`` as the ``synonyms`` column stores it.

    Keeps non-empty strings, stripped, each once ignoring case, in order. None
    if nothing is left.
    """
    seen = set()
    kept = []
    for name in names:
        if isinstance(name, str) and (name := name.strip()) and name.casefold() not in seen:
            seen.add(name.casefold())
            kept.append(name)
    return kept or None


@dataclass
class GenomeCreationResult:
    """What :meth:`GenomeManager.ensure_from_remote` did, so a pull can undo it."""

    success: bool
    genome_digest: GenomeDigest | None = None
    alias_name: GenomeAlias | None = None
    created_genome: bool = False
    created_alias: bool = False
    #: The server's record of the genome, if this call fetched it, so the
    #: caller need not ask the server again.
    genome_attrs: dict | None = None


def _compare_level2(level2_a: dict, level2_b: dict, digest_a: str, digest_b: str) -> dict:
    """Compare two level2 collection dicts locally using seqcol comparison logic."""
    attributes = {}
    array_elements = {}

    for attr in ("names", "lengths", "sequences"):
        arr_a = level2_a.get(attr, [])
        arr_b = level2_b.get(attr, [])
        set_a = set(arr_a) if isinstance(arr_a, list) else set()
        set_b = set(arr_b) if isinstance(arr_b, list) else set()

        attributes[attr] = {
            "a_count": len(arr_a),
            "b_count": len(arr_b),
            "a_and_b": len(set_a & set_b),
            "a_only": len(set_a - set_b),
            "b_only": len(set_b - set_a),
        }

        array_elements[attr] = {
            "a_and_b": sorted(set_a & set_b),
            "a_only": sorted(set_a - set_b),
            "b_only": sorted(set_b - set_a),
        }

    return {
        "digests": {"a": digest_a, "b": digest_b},
        "attributes": attributes,
        "array_elements": array_elements,
    }


class GenomeStoreRouter(Protocol):
    """What genome reads need from the root's store router.

    ``refgenie.core.store_router.RefgetStoreRouter`` satisfies it. Described
    structurally so this manager does not import up into ``refgenie.core``.
    """

    def store_for_genome(self, genome_digest: str, store_name: str | None = None) -> Any: ...
    def get_fhr_metadata(self, genome_digest: str, store_name: str | None = None) -> Any: ...


class GenomeManager(ResourceManager):
    """Manager for genome operations."""

    def __init__(
        self,
        database_engine: Engine,
        refget_store_getter: Callable[[], "RefgetStore"],
        alias_manager_getter: "Callable[[], AliasBackend]",
        servers_getter: "Callable[[], ServerManager]",
        events: EventSink | None = None,
        *,
        store_router_getter: Callable[[], GenomeStoreRouter],
        species_manager: "SpeciesManager",
    ):
        """
        Initialize the GenomeManager.

        The alias manager and the server manager are passed as getters because
        ``Refgenie`` builds them lazily, after this manager.

        Args:
            database_engine: The database engine.
            refget_store_getter: Returns the writable RefgetStore.
            alias_manager_getter: Returns the alias manager.
            servers_getter: Returns the ServerManager, for genomes created from servers.
            events: Where genome changes are recorded for plugins.
            store_router_getter: Returns the store router. Genome reads go
                through it, so each genome is read from the store that owns it.
            species_manager: The species list, for listings' species and a
                pulled genome's species row.
        """
        super().__init__(database_engine)
        self._refget_store_getter = refget_store_getter
        self._store_router_getter = store_router_getter
        self._alias_manager_getter = alias_manager_getter
        self._servers_getter = servers_getter
        self._species = species_manager
        self._events = events or NULL_EVENTS

    @update_scope
    def add(
        self,
        digest: GenomeDigest,
        description: str,
        alias_names: list[GenomeAlias],
        species_name: str | None = None,
        remote_url: str | None = None,
        store_name: str | None = None,
    ) -> Genome:
        """Add a genome."""
        with self._database_session as session:
            genome = Genome(
                digest=digest,
                description=description,
                species_name=species_name,
                remote_url=remote_url,
                store_name=store_name,
            )
            session.add(genome)
            session.commit()
            session.refresh(genome)
            logger.info(f"Added genome: {digest}")
        self._events.record(Change(action="genome_added", genome=digest))

        for name in alias_names:
            self._alias_manager_getter().add(name, digest)

        return genome

    def apply_fhr(self, digest: str, fhr_dict: dict) -> None:
        """Apply a normalized FHR record to the genome row and the store sidecar.

        Writes the ``genome`` columns (see :meth:`apply_fhr_columns`) and the
        RefgetStore FHR sidecar. Every FHR write goes through here so the row
        and the sidecar cannot diverge. The genome row must already exist.
        """
        self.apply_fhr_columns(digest, fhr_dict)
        self._refget_store_getter().set_fhr_metadata(digest, self._fhr_metadata_from_dict(fhr_dict))
        logger.info(f"Applied FHR metadata to genome {digest}")

    def apply_fhr_columns(self, digest: str, fhr_dict: dict) -> None:
        """Write the queryable ``genome`` columns from an FHR record, no sidecar.

        Maps ``documentation``, ``taxon.name``, ``taxon.uri``, ``genomeSynonym``,
        ``assemblySource``, ``accessionID.name`` and ``assemblyLevel``. Only
        fields present in the record are written, so a partial record never
        nulls existing columns. An empty ``genomeSynonym`` clears the
        synonyms. Federated ``store sync`` uses this directly, since the
        owning store already holds the sidecar.
        """
        taxon = fhr_dict.get("taxon") or {}
        synonyms = fhr_dict.get("genomeSynonym")
        values = {
            "description": fhr_dict.get("documentation"),
            "species_name": taxon.get("name"),
            "synonyms": None if synonyms is None else _clean_synonyms(synonyms),
            "taxon_id": self._taxon_id_from_uri(taxon.get("uri")),
            "assembly_source": fhr_dict.get("assemblySource"),
            "assembly_accession": (fhr_dict.get("accessionID") or {}).get("name"),
            "assembly_level": fhr_dict.get("assemblyLevel"),
        }
        updates = {column: value for column, value in values.items() if value is not None}
        if synonyms is not None:
            # A record that lists genomeSynonym, even empty, states the full set.
            updates["synonyms"] = values["synonyms"]

        with self._database_session as session:
            genome = queries.get(session, digest)
            for column, value in updates.items():
                setattr(genome, column, value)
            session.add(genome)
            session.commit()

    @update_scope
    def fill_from_server(self, digest: str, attrs: dict) -> None:
        """Copy a server's genome record onto the local genome, where it is unset.

        Fills ``description``, ``species_name``, ``synonyms``, ``taxon_id``,
        ``assembly_source``, ``assembly_accession`` and ``assembly_level`` from
        ``attrs`` (the server's ``/genomes/{digest}`` body), and stores the
        server's FHR record as the local FHR sidecar once the collection is in
        the local store. A value the user already set locally is never
        overwritten, and neither is an existing sidecar.

        The server's ``species`` row is stored every time, replacing a local
        copy of that taxon. It carries the common name: genomes have none of
        their own.

        Args:
            digest: The local genome to fill. It must already exist.
            attrs: The server's genome attributes; may be empty.
        """
        if not attrs:
            return
        if attrs.get("species"):
            try:
                self._species.replace([SpeciesEntry.model_validate(attrs["species"])])
            except ValueError as exc:  # a malformed record from the server
                logger.warning(f"Could not store the server's species record for {digest}: {exc}")
        with self._database_session as session:
            genome = queries.get(session, digest)
            changed = False
            for column in SERVER_GENOME_COLUMNS:
                value = attrs.get(column)
                if value in (None, "", []) or not _is_unset(getattr(genome, column), column):
                    continue
                setattr(genome, column, value)
                changed = True
            if changed:
                session.add(genome)
                session.commit()

        fhr = attrs.get("fhr")
        if not isinstance(fhr, dict) or not fhr:
            return
        store = self._refget_store_getter()
        try:
            store.get_collection_metadata(digest)
        except Exception:  # noqa: BLE001 - not in the local store yet
            # The sidecar hangs off the collection; a pull that brings the
            # fasta adds the collection and calls this again.
            logger.debug(f"Collection {digest} not in the local store; FHR record not stored")
            return
        try:
            if store.get_fhr_metadata(digest) is not None:
                return
            store.set_fhr_metadata(digest, self._fhr_metadata_from_dict(fhr))
        except Exception as exc:  # noqa: BLE001 - the sidecar is best-effort
            logger.warning(f"Could not store the server's FHR record for {digest}: {exc}")

    def set_store_name(self, digest: str, store_name: str | None) -> None:
        """Set the owning store of a genome (federation)."""
        with self._database_session as session:
            genome = queries.get(session, digest)
            genome.store_name = store_name
            session.add(genome)
            session.commit()

    @staticmethod
    def _taxon_id_from_uri(uri: str | None) -> int | None:
        """Parse the id from ``https://identifiers.org/taxonomy:<id>``, or None."""
        if not uri:
            return None
        try:
            return int(uri.rstrip("/").split(":")[-1])
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _fhr_metadata_from_dict(fhr_dict: dict):
        """Build a gtars ``FhrMetadata`` from a camelCase FHR dict.

        gtars only reads FHR from a file, so the dict goes through a temp file.
        """
        import json
        import os

        from gtars.refget import FhrMetadata

        fd, tmp = tempfile.mkstemp(suffix=".fhr.json")
        try:
            with os.fdopen(fd, "w") as fh:
                json.dump(fhr_dict, fh)
            return FhrMetadata.from_json(tmp)
        finally:
            os.unlink(tmp)

    def get(self, digest: GenomeDigest) -> Genome:
        """Get a genome by its digest, with its SQL alias rows loaded."""
        with self._database_session as session:
            return queries.get(session, digest, with_aliases=True)

    def page(
        self,
        *,
        digest: GenomeDigest | None = None,
        alias: GenomeAlias | None = None,
        search: Search | None = None,
        has_assets: bool | None = None,
        taxon_id: int | None = None,
        offset: int,
        limit: int,
    ) -> Page[GenomeListing]:
        """
        One page of genomes with their aliases, asset counts and species rows.

        Genomes with assets come first, then by species name, then digest.
        Aliases come from the mode-selected alias backend (the RefgetStore in
        local mode), read for the whole page at once, and so does an
        ``alias`` filter and the ``aliases`` search field.

        Args:
            digest: Only this genome.
            alias: Only the genome this alias names. An unknown alias gives an
                empty page, not an error.
            search: Match the genome search fields (``queries.SEARCH_FIELDS``).
            has_assets: Only genomes with (True) or without (False) assets.
            taxon_id: Only genomes of this NCBI taxon.
            offset: Rows to skip.
            limit: Most rows to return.

        Raises:
            ValueError: If the search names a field that is not searchable.
        """
        alias_mgr = self._alias_manager_getter()
        if alias is not None:
            try:
                resolved = alias_mgr.resolve(alias)
            except MissingAliasError:
                return Page(items=[], offset=offset, limit=limit, total=0)
            if digest is not None and digest != resolved:
                return Page(items=[], offset=offset, limit=limit, total=0)
            digest = resolved
        alias_matches = []
        if queries.searches_aliases(search):
            alias_matches = alias_digests_matching(alias_mgr, search.term.strip(), search.operator)
        with self._database_session as session:
            genomes, total = queries.list_page(
                session,
                digest=digest,
                search=search,
                alias_matches=alias_matches,
                has_assets=has_assets,
                taxon_id=taxon_id,
                offset=offset,
                limit=limit,
            )
            digests = [genome.digest for genome in genomes]
            counts = queries.count_asset_groups(session, digests)
        aliases = alias_mgr.get_for_genomes(digests)
        species = self._species.get_many({g.taxon_id for g in genomes if g.taxon_id is not None})
        listings = [
            GenomeListing(
                genome=genome,
                aliases=aliases.get(genome.digest, []),
                asset_count=counts.get(genome.digest, 0),
                species=species.get(genome.taxon_id),
            )
            for genome in genomes
        ]
        return Page(items=listings, offset=offset, limit=limit, total=total)

    def detail(self, digest: GenomeDigest) -> GenomeDetail:
        """
        One genome as its detail page shows it: listing fields (species row included) plus FHR.

        The session is closed before the store is read for the FHR record, so
        no connection is held during the store lookup.

        Raises:
            MissingGenomeError: If there is no genome with this digest.
        """
        with self._database_session as session:
            genome = queries.get(session, digest)
            count = queries.count_asset_groups(session, [digest]).get(digest, 0)
        aliases = self._alias_manager_getter().get_for_genomes([digest]).get(digest, [])
        species = self._species.get(genome.taxon_id) if genome.taxon_id is not None else None
        fhr = self._store_router_getter().get_fhr_metadata(digest, store_name=genome.store_name)
        return GenomeDetail(
            genome=genome,
            aliases=aliases,
            asset_count=count,
            species=species,
            fhr=fhr.to_dict() if fhr is not None else None,
        )

    @update_scope
    def remove(self, digest: GenomeDigest) -> None:
        """
        Remove a genome and all its aliases, asset groups, and assets.

        The catalog commits first, then aliases are removed. If interrupted in
        between, an alias can be left naming a missing genome; ``refgenie alias
        remove <name>`` clears it.
        """
        # Read before the delete, which leaves nothing to enumerate them from.
        alias_mgr = self._alias_manager_getter()
        alias_names = list(alias_mgr.get_for_genome(digest))

        with self._database_session as session:
            genome = queries.get(session, digest)
            session.delete(genome)
            session.commit()
        logger.info(
            f"Removed genome '{digest}' and all corresponding aliases, asset groups and assets"
        )
        self._events.record(Change(action="genome_removed", genome=digest))

        # One failed alias removal must not stop the rest.
        for alias_name in alias_names:
            try:
                alias_mgr.remove(alias_name)
            except MissingAliasError:
                logger.debug(f"Alias '{alias_name}' was already gone")
            except Exception:
                logger.exception(
                    f"Could not remove alias '{alias_name}' for the removed genome "
                    f"'{digest}'; it now names a genome that is gone. "
                    f"Remove it with `refgenie alias remove {alias_name}`."
                )

    def list_all(self) -> Iterable[Genome]:
        """List all genomes, assets-first, then by species and digest."""
        with self._database_session as session:
            return queries.list_all(session)

    def exists(self, digest: GenomeDigest) -> bool:
        """Check if a genome exists."""
        with self._database_session as session:
            return queries.exists(session, digest)

    def has_collection(self, digest: GenomeDigest) -> bool:
        """Whether the genome's sequence collection is in the local RefgetStore."""
        try:
            return self._refget_store_getter().get_collection_metadata(digest) is not None
        except Exception:  # noqa: BLE001 - a missing collection can raise
            return False

    def table(self) -> Table:
        """Create a Rich table of all genomes."""
        genomes = self.list_all()
        aliases = self._alias_manager_getter().get_for_genomes([g.digest for g in genomes])
        rows = [
            (
                genome.digest,
                ", ".join(aliases[genome.digest]),
                f"remote:{genome.remote_url}" if genome.remote_url else "local",
                genome.species_name or "",
                genome.description or "",
            )
            for genome in genomes
        ]
        return build_table(
            "Genomes", ["Digest", "Aliases", "Source", "Species", "Description"], rows
        )

    def _owning_store(self, genome: Genome):
        """The store that owns ``genome``, with its collection loaded."""
        return self._store_router_getter().store_for_genome(genome.digest, genome.store_name)

    def get_metadata(self, digest: GenomeDigest) -> dict:
        """Get metadata for a genome from the RefgetStore that owns it."""
        genome = self.get(digest)
        level2 = self._owning_store(genome).get_collection_level2(digest)

        if level2 is None:
            raise ValueError(f"No collection metadata found for {digest}")

        total_length = sum(level2["lengths"])

        result = {
            "digest": digest,
            "n_sequences": len(level2["names"]),
            "total_length": total_length,
            "source": "local",
        }

        if genome.remote_url:
            result["source"] = "remote"
            result["remote_url"] = genome.remote_url

        return result

    def compare(self, digest_a: GenomeDigest, digest_b: GenomeDigest) -> dict:
        """Compare two genomes by their digests."""
        genome_a = self.get(digest_a)
        genome_b = self.get(digest_b)

        if genome_a.remote_url or genome_b.remote_url:
            return self._compare_with_remote(genome_a, genome_b)
        # Each genome may live in a different store of the router; the store's
        # own compare needs both collections in the same store.
        store_a = self._owning_store(genome_a)
        store_b = self._owning_store(genome_b)
        if store_a is store_b:
            return store_a.compare(digest_a, digest_b)
        return _compare_level2(
            store_a.get_collection_level2(digest_a),
            store_b.get_collection_level2(digest_b),
            digest_a,
            digest_b,
        )

    def _compare_with_remote(self, genome_a: Genome, genome_b: Genome) -> dict:
        """Compare genomes where at least one is remote."""
        if genome_a.remote_url and genome_a.remote_url == genome_b.remote_url:
            with open_remote_store(genome_a.remote_url) as remote_store:
                return remote_store.compare(genome_a.digest, genome_b.digest)
        level2_a, level2_b = (
            self._get_remote_level2(g.remote_url, g.digest)
            if g.remote_url
            else self._owning_store(g).get_collection_level2(g.digest)
            for g in (genome_a, genome_b)
        )
        return _compare_level2(level2_a, level2_b, genome_a.digest, genome_b.digest)

    @staticmethod
    def _get_remote_level2(remote_url: str, digest: str) -> dict:
        """Fetch level2 data from a remote RefgetStore."""
        with open_remote_store(remote_url) as remote_store:
            return remote_store.get_collection_level2(digest)

    @update_scope
    def initialize_genome(
        self,
        description: str,
        alias_names: list[GenomeAlias],
        fasta_file_path: Path | None = None,
        source: RemoteGenomeSource | None = None,
        digest: GenomeDigest | None = None,
        use_existing: bool = False,
        species_name: str | None = None,
        expected_digest: GenomeDigest | None = None,
    ) -> tuple[GenomeDigest, bool]:
        """Initialize a genome from a local FASTA file or a remote source.

        ``expected_digest`` (local FASTA only) makes a FASTA whose digest
        differs fail before any genome row or alias is written.
        """
        has_local = fasta_file_path is not None
        has_remote = source is not None or digest is not None
        if has_local and has_remote:
            raise ValueError("Cannot specify both fasta_file_path and source/digest.")
        if not has_local and not has_remote:
            raise ValueError("Must specify either fasta_file_path or source + digest.")
        if has_remote and source is None:
            raise ValueError("A RemoteGenomeSource is required for remote init.")

        if has_local:
            return self._initialize_genome_local(
                fasta_file_path=fasta_file_path,
                description=description,
                alias_names=alias_names,
                use_existing=use_existing,
                species_name=species_name,
                expected_digest=expected_digest,
            )
        else:
            return self._initialize_genome_remote(
                source=source,
                digest=digest,
                description=description,
                alias_names=alias_names,
                use_existing=use_existing,
                species_name=species_name,
            )

    def _initialize_genome_local(
        self,
        fasta_file_path: Path,
        description: str,
        alias_names: list[GenomeAlias],
        use_existing: bool = False,
        species_name: str | None = None,
        expected_digest: GenomeDigest | None = None,
    ) -> tuple[GenomeDigest, bool]:
        """Initialize a genome from a local FASTA file."""
        if not fasta_file_path.exists():
            raise FileNotFoundError(
                f"FASTA file not found: {fasta_file_path}. Cannot initialize genome."
            )
        logger.info(f"Initializing genome from FASTA file: {fasta_file_path}")

        store = self._refget_store_getter()
        metadata, was_new = store.add_sequence_collection_from_fasta(fasta_file_path)
        digest = GenomeDigest(metadata.digest)
        logger.debug(f"{fasta_file_path} digest: {digest}")
        if expected_digest is not None and digest != expected_digest:
            if was_new:
                store.remove_collection(digest)
            raise FastaDigestMismatchError(fasta_file_path.name, expected_digest, digest)

        return digest, self._register(
            digest, description, alias_names, use_existing, species_name=species_name
        )

    def _register(
        self,
        digest: GenomeDigest,
        description: str,
        alias_names: list[GenomeAlias],
        use_existing: bool,
        **columns,
    ) -> bool:
        """Add the genome row, or reuse an existing one if allowed; True if added."""
        if not self.exists(digest):
            self.add(digest=digest, description=description, alias_names=alias_names, **columns)
            return True
        if not use_existing:
            raise ValueError(
                f"Genome already exists with digest {digest}. Can't reinitialize genomes."
            )
        logger.debug(f"Genome already exists: {digest}. Using existing.")
        self._repair_missing_aliases(digest, alias_names)
        return False

    def _repair_missing_aliases(self, digest: str, alias_names: list[str]) -> None:
        """Re-add any of ``alias_names`` that no longer resolve to ``digest``.

        Genome rows live in SQL but aliases live in the RefgetStore, whose
        unlocked index can lose an alias, so re-init repairs them.
        """
        if not alias_names:
            return
        alias_mgr = self._alias_manager_getter()
        for name in alias_names:
            try:
                if alias_mgr.resolve(name) == digest:
                    continue
            except Exception:  # noqa: BLE001 - unresolvable reads as missing
                pass
            logger.warning(
                f"Alias '{name}' did not resolve to existing genome {digest}; re-adding it."
            )
            alias_mgr.add(name, digest)

    def _initialize_genome_remote(
        self,
        source: RemoteGenomeSource,
        description: str,
        alias_names: list[GenomeAlias],
        digest: GenomeDigest | None = None,
        use_existing: bool = False,
        species_name: str | None = None,
    ) -> tuple[GenomeDigest, bool]:
        """Initialize a genome from a remote source."""
        if digest is None:
            if alias_names:
                resolved = source.resolve_alias(alias_names[0])
                if resolved is not None:
                    digest = resolved
                    logger.info(f"Resolved alias '{alias_names[0]}' to digest: {digest}")
            if digest is None:
                raise ValueError(
                    "No digest provided and could not resolve alias. "
                    "Provide --digest or use a source that supports alias resolution."
                )

        try:
            collection_data = source.verify_collection(digest)
        except ConnectionError as e:
            raise ConnectionError(f"Cannot reach source at {source.url}: {e}") from e

        if collection_data is None:
            raise ValueError(
                f"Collection {digest} not found on source {source.url}. "
                f"Use 'refgenie genome browse' to see available collections."
            )

        n_sequences = len(collection_data.get("names", []))
        logger.info(f"Verified collection {digest} on {source.url} ({n_sequences} sequences)")

        # Import the collection into the local store so aliases resolve.
        # TODO: Replace with store.add_remote(url); store.load_collection(digest)
        # once gtars supports multi-remote.
        remote_url = source.store_url
        if remote_url is not None:
            self._import_remote_collection(remote_url, digest)
        else:
            logger.warning(
                f"Source {source.url} does not advertise a store URL. "
                f"Collection not imported to local store; aliases may not resolve."
            )

        return digest, self._register(
            digest,
            description,
            alias_names,
            use_existing,
            species_name=species_name,
            remote_url=remote_url,
        )

    def _import_remote_collection(self, remote_url: str, digest: str) -> None:
        """Import a collection from a remote store into the local store."""
        try:
            with open_remote_store(remote_url) as remote_store:
                collection = remote_store.get_collection(digest)
                local_store = self._refget_store_getter()
                local_store.add_sequence_collection(collection)
                logger.info(f"Imported collection {digest} into local store")
        except Exception as e:
            logger.warning(
                f"Could not import collection {digest} from {remote_url}: {e}. "
                f"Aliases may not resolve until the collection is available locally."
            )

    # =========================================================================
    # Bootstrapping a genome from subscribed servers
    # =========================================================================

    @update_scope
    def ensure_from_remote(
        self,
        alias_name: GenomeAlias,
        genome_digest: GenomeDigest | None = None,
        genome_description: str | None = None,
        server_urls: list[str] | None = None,
        genome_attrs: dict | None = None,
    ) -> GenomeCreationResult:
        """
        Make a genome and its alias exist locally, asking the servers if needed.

        Resolves the alias against the servers, creates what is missing, and
        reports what it created so a failed pull can roll it back. Where
        possible the genome is created through a ``RemoteGenomeSource``, so
        ``remote_url`` is set and sequences can be fetched lazily.

        Args:
            alias_name: The name of the alias.
            genome_digest: The digest of the genome, if already known.
            genome_description: The description of the genome.
            server_urls: Servers to ask. Defaults to the subscriptions.
            genome_attrs: The server's genome record, if the caller already
                has it. Fetched here when needed and not given.

        Returns:
            GenomeCreationResult: Success flag, what this call created, and
            the server's genome record if it was fetched or given.
        """
        servers = self._servers_getter()
        alias_manager = self._alias_manager_getter()

        if genome_digest is None:
            genome_digest = servers.resolve_alias(alias_name, server_urls)
            if genome_digest is None:
                return GenomeCreationResult(success=False)
            logger.info(f"Determined digest for {alias_name}: {genome_digest}")

        # Checked before any writes, so a rollback never deletes rows that
        # existed before the pull.
        genome_existed = self.exists(genome_digest)
        alias_existed = alias_manager.exists(alias_name)

        if genome_existed:
            alias_manager.add(alias_name, genome_digest)
            return GenomeCreationResult(
                success=True,
                genome_digest=genome_digest,
                alias_name=alias_name,
                created_genome=False,
                created_alias=not alias_existed,
                genome_attrs=genome_attrs,
            )

        if genome_attrs is None:
            genome_attrs = servers.genome_attrs(genome_digest, server_urls)
        genome_description = genome_description or genome_attrs.get("description") or NO_DESCRIPTION

        created = True
        source = servers.genome_source(server_urls)
        initialized = False
        if source is not None:
            try:
                _digest, created = self.initialize_genome(
                    source=source,
                    digest=genome_digest,
                    description=genome_description,
                    alias_names=[alias_name],
                    use_existing=True,
                )
                initialized = True
            except Exception as exc:
                logger.debug(f"Source-based genome init failed, falling back: {exc}")
        if not initialized:
            self.add(genome_digest, genome_description, [alias_name])
        return GenomeCreationResult(
            success=True,
            genome_digest=genome_digest,
            alias_name=alias_name,
            created_genome=created,
            created_alias=not alias_existed,
            genome_attrs=genome_attrs,
        )

    @update_scope
    def init_from_remote(
        self,
        alias_name: GenomeAlias,
        genome_digest: GenomeDigest | None = None,
        genome_description: str | None = None,
        server_urls: list[str] | None = None,
        require_collection: bool = False,
    ) -> bool:
        """
        Register a genome locally from remote metadata, downloading no assets.

        Args:
            alias_name: The alias name to look up on remote servers.
            genome_digest: The digest of the genome, if already known.
            genome_description: The description of the genome.
            server_urls: Optional list of server URLs to query.
            require_collection: If True, undo the registration and return False
                when the genome's sequence collection is not in the local store.

        Returns:
            True if the alias exists afterwards, False if no server knew it, or
            the collection was required and is missing.
        """
        if self._alias_manager_getter().exists(alias_name):
            logger.info(f"Alias '{alias_name}' already exists locally")
            return True

        result = self.ensure_from_remote(
            alias_name=alias_name,
            genome_digest=genome_digest,
            genome_description=genome_description,
            server_urls=server_urls,
        )
        if result.success:
            if require_collection and not self.has_collection(result.genome_digest):
                if result.created_genome:
                    self.remove(result.genome_digest)  # also removes its aliases
                elif result.created_alias:
                    self._alias_manager_getter().remove(alias_name)
                logger.error(
                    f"Could not register '{alias_name}': its sequence collection is not in "
                    f"the server's sequence store. Run 'refgenie pull {alias_name}/fasta' instead."
                )
                return False
            if result.genome_attrs is not None:
                self.fill_from_server(result.genome_digest, result.genome_attrs)
            logger.info(f"Registered genome for alias '{alias_name}': {result.genome_digest}")
            return True
        logger.error(f"Could not find alias '{alias_name}' on remote servers")
        return False
