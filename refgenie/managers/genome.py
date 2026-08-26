"""
GenomeManager - handles genome CRUD operations.
"""

from pathlib import Path
from collections.abc import Callable, Iterable
from typing import TYPE_CHECKING

from rich.table import Table
from sqlalchemy.engine import Engine
from sqlalchemy.orm import selectinload
from sqlmodel import select

from refgenie.db.tables import AssetGroup, Genome
from refgenie.exceptions import MissingAliasError, MissingGenomeError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.utils.tables import build_table
from refgenie.managers.sources import RemoteGenomeSource

if TYPE_CHECKING:
    from gtars.refget import RefgetStore


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


class GenomeManager(ResourceManager):
    """Manager for genome operations."""

    def __init__(
        self,
        database_engine: Engine,
        refget_store_getter: Callable[[], "RefgetStore"],
        alias_manager_getter: Callable | None = None,
    ):
        super().__init__(database_engine)
        self._refget_store_getter = refget_store_getter
        self._alias_manager_getter = alias_manager_getter

    def add(
        self,
        digest: str,
        description: str,
        alias_names: list[str],
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

        # Add aliases via the alias manager (works for both SQL and store backends)
        if alias_names and self._alias_manager_getter:
            alias_mgr = self._alias_manager_getter()
            for name in alias_names:
                alias_mgr.add(name, digest)

        return genome

    def apply_fhr(self, digest: str, fhr_dict: dict) -> None:
        """Apply a normalized FHR record to the genome row and the store sidecar.

        The single source of truth for the FHR -> catalog mapping. It:

          * upserts the ``genome`` row's metadata columns from the FHR record
            (``description`` <- ``documentation``; ``species_name`` <- ``genome``,
            i.e. the scientific name; plus the queryable assembly/organism core:
            ``common_name`` <- ``commonName``, ``taxon_id`` <- the id in
            ``taxon.uri``, ``assembly_source`` <- ``assemblySource``,
            ``assembly_accession`` <- ``accessionID.name``, ``assembly_level`` <-
            ``assemblyLevel``); and
          * writes the RefgetStore FHR sidecar via ``set_fhr_metadata`` so
            ``get_fhr_metadata`` (the server) and the S3 mirror carry the record.

        All FHR write paths (including ``genome init --fhr`` and
        ``genome set-metadata --fhr``) funnel through here so the row and the
        sidecar cannot diverge.

        Idempotent; the genome row must already exist (``genome init``/``add``
        creates it) -- ``apply_fhr`` updates, it never inserts. Only fields
        present in the record are written, so a partial record never nulls
        existing columns.
        """
        self.apply_fhr_columns(digest, fhr_dict)
        self._refget_store_getter().set_fhr_metadata(
            digest, self._fhr_metadata_from_dict(fhr_dict)
        )
        logger.info(f"Applied FHR metadata to genome {digest}")

    def apply_fhr_columns(self, digest: str, fhr_dict: dict) -> None:
        """Write the queryable ``genome`` columns from an FHR record.

        The column half of :meth:`apply_fhr`, with no RefgetStore sidecar write.
        Federated ``store sync`` uses this: the owning store already holds the
        FHR sidecar, so the server reads it through the router at serve time and
        only the faceted columns need to land on the ``genome`` row here. Only
        fields present in the record are written; a partial record never nulls
        existing columns.
        """
        updates: dict = {}
        documentation = fhr_dict.get("documentation")
        if documentation is not None:
            updates["description"] = documentation
        species_name = fhr_dict.get("genome")
        if species_name is not None:
            updates["species_name"] = species_name
        common_name = fhr_dict.get("commonName")
        if common_name is not None:
            updates["common_name"] = common_name
        taxon_id = self._taxon_id_from_uri((fhr_dict.get("taxon") or {}).get("uri"))
        if taxon_id is not None:
            updates["taxon_id"] = taxon_id
        assembly_source = fhr_dict.get("assemblySource")
        if assembly_source is not None:
            updates["assembly_source"] = assembly_source
        assembly_accession = (fhr_dict.get("accessionID") or {}).get("name")
        if assembly_accession is not None:
            updates["assembly_accession"] = assembly_accession
        assembly_level = fhr_dict.get("assemblyLevel")
        if assembly_level is not None:
            updates["assembly_level"] = assembly_level

        with self._database_session as session:
            genome = one_or_raise(
                session,
                select(Genome).where(Genome.digest == digest),
                MissingGenomeError(genome=digest),
            )
            for column, value in updates.items():
                setattr(genome, column, value)
            session.add(genome)
            session.commit()

    def set_store_name(self, digest: str, store_name: str | None) -> None:
        """Set the owning store of a genome (federation)."""
        with self._database_session as session:
            genome = one_or_raise(
                session,
                select(Genome).where(Genome.digest == digest),
                MissingGenomeError(genome=digest),
            )
            genome.store_name = store_name
            session.add(genome)
            session.commit()

    @staticmethod
    def _taxon_id_from_uri(uri: str | None) -> int | None:
        """Parse the NCBI taxon id out of an identifiers.org taxonomy URI.

        Expects an identifiers.org taxonomy URI
        (``https://identifiers.org/taxonomy:<id>``) and extracts the trailing
        id into the queryable ``taxon_id`` column (the URI is re-derived from
        the id for display). Returns None for a missing or unparseable URI.
        """
        if not uri:
            return None
        try:
            return int(uri.rstrip("/").split(":")[-1])
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _fhr_metadata_from_dict(fhr_dict: dict):
        """Build a gtars ``FhrMetadata`` from a camelCase FHR dict.

        gtars only deserializes FHR from a file (``FhrMetadata.from_json(path)``),
        so the dict is round-tripped through a temp file. The file is the exact
        camelCase shape gtars serializes, so nothing is lost.
        """
        import json
        import os
        import tempfile

        from gtars.refget import FhrMetadata

        fd, tmp = tempfile.mkstemp(suffix=".fhr.json")
        try:
            with os.fdopen(fd, "w") as fh:
                json.dump(fhr_dict, fh)
            return FhrMetadata.from_json(tmp)
        finally:
            os.unlink(tmp)

    def get(self, digest: str) -> Genome:
        """Get a genome by its digest."""
        statement = (
            select(Genome)
            .options(
                selectinload(Genome.aliases),
                selectinload(Genome.asset_groups).selectinload(AssetGroup.assets),
            )
            .where(Genome.digest == digest)
        )
        with self._database_session as session:
            return one_or_raise(session, statement, MissingGenomeError(genome=digest), unique=True)

    def remove(self, digest: str) -> None:
        """
        Remove a genome and all its aliases, asset groups, and assets.

        The catalog commits first. If interrupted afterward, store aliases may
        point at a gone genome -- drift a repair pass can clean and no read
        path mistakes for a working genome.
        """
        # Resolved before the delete: after it, the genome the aliases belong to
        # is gone and there is nothing left to enumerate them from.
        alias_mgr = self._alias_manager_getter() if self._alias_manager_getter else None
        alias_names = list(alias_mgr.get_for_genome(digest)) if alias_mgr else []

        with self._database_session as session:
            # Queried through this session rather than self.get(): _database_session
            # must not be re-entered from inside an open block.
            genome = one_or_raise(
                session,
                select(Genome).where(Genome.digest == digest),
                MissingGenomeError(genome=digest),
            )
            session.delete(genome)
            session.commit()
        logger.info(
            f"Removed genome '{digest}' and all corresponding aliases, asset groups and assets"
        )

        # Cleanup, after the catalog has decided. Each alias removal is
        # independent, so one failure must not strand the rest.
        for alias_name in alias_names:
            try:
                alias_mgr.remove(alias_name)
            except MissingAliasError:
                logger.debug(f"Alias '{alias_name}' was already gone")

    def list_all(self) -> Iterable[Genome]:
        """List all genomes."""
        with self._database_session as session:
            return (
                session.exec(
                    select(Genome).options(
                        selectinload(Genome.aliases),
                        selectinload(Genome.asset_groups).selectinload(AssetGroup.assets),
                    )
                )
                .unique()
                .all()
            )

    def exists(self, digest: str) -> bool:
        """Check if a genome exists."""
        statement = select(Genome).where(Genome.digest == digest)
        with self._database_session as session:
            result = session.exec(statement)
            return bool(result.first())

    def table(self) -> Table:
        """Create a Rich table of all genomes."""
        alias_mgr = self._alias_manager_getter() if self._alias_manager_getter else None
        rows = []
        with self._database_session as session:
            genomes = session.exec(select(Genome)).unique().all()
            for genome in genomes:
                alias_names = alias_mgr.get_for_genome(genome.digest) if alias_mgr else []
                rows.append(
                    (
                        genome.digest,
                        ", ".join(alias_names) if alias_names else "",
                        f"remote:{genome.remote_url}" if genome.remote_url else "local",
                        genome.species_name or "",
                        genome.description or "",
                    )
                )
        return build_table(
            "Genomes", ["Digest", "Aliases", "Source", "Species", "Description"], rows
        )

    def get_metadata(self, digest: str) -> dict:
        """Get metadata for a genome from the RefgetStore."""
        genome = self.get(digest)
        store = self._refget_store_getter()
        level2 = store.get_collection_level2(digest)

        if level2 is None:
            raise ValueError(f"No collection metadata found for {digest}")

        total_length = sum(level2["lengths"])

        result = {
            "digest": digest,
            "n_sequences": len(level2["names"]),
            "total_length": total_length,
            "source": "local",
        }

        if hasattr(genome, "remote_url") and genome.remote_url:
            result["source"] = "remote"
            result["remote_url"] = genome.remote_url

        return result

    def compare(self, digest_a: str, digest_b: str) -> dict:
        """Compare two genomes by their digests."""
        genome_a = self.get(digest_a)
        genome_b = self.get(digest_b)

        if genome_a.remote_url or genome_b.remote_url:
            return self._compare_with_remote(genome_a, genome_b)
        else:
            store = self._refget_store_getter()
            return store.compare(digest_a, digest_b)

    def _compare_with_remote(self, genome_a: Genome, genome_b: Genome) -> dict:
        """Compare genomes where at least one is remote."""
        store = self._refget_store_getter()

        if genome_a.remote_url:
            level2_a = self._get_remote_level2(genome_a.remote_url, genome_a.digest)
        else:
            level2_a = store.get_collection_level2(genome_a.digest)

        if genome_b.remote_url:
            level2_b = self._get_remote_level2(genome_b.remote_url, genome_b.digest)
        else:
            level2_b = store.get_collection_level2(genome_b.digest)

        if (
            genome_a.remote_url
            and genome_b.remote_url
            and genome_a.remote_url == genome_b.remote_url
        ):
            import tempfile

            from gtars.refget import RefgetStore as GtarsRefgetStore

            with tempfile.TemporaryDirectory() as tmpdir:
                remote_store = GtarsRefgetStore.open_remote(tmpdir, genome_a.remote_url)
                return remote_store.compare(genome_a.digest, genome_b.digest)

        return _compare_level2(level2_a, level2_b, genome_a.digest, genome_b.digest)

    @staticmethod
    def _get_remote_level2(remote_url: str, digest: str) -> dict:
        """Fetch level2 data from a remote RefgetStore."""
        import tempfile

        from gtars.refget import RefgetStore as GtarsRefgetStore

        with tempfile.TemporaryDirectory() as tmpdir:
            remote_store = GtarsRefgetStore.open_remote(tmpdir, remote_url)
            return remote_store.get_collection_level2(digest)

    def initialize_genome(
        self,
        description: str,
        alias_names: list[str],
        fasta_file_path: Path | None = None,
        source: RemoteGenomeSource | None = None,
        digest: str | None = None,
        use_existing: bool = False,
        species_name: str | None = None,
    ) -> tuple[str, bool]:
        """Initialize a genome from a local FASTA file or a remote source."""
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
        alias_names: list[str],
        use_existing: bool = False,
        species_name: str | None = None,
    ) -> tuple[str, bool]:
        """Initialize a genome from a local FASTA file."""
        if not fasta_file_path.exists():
            raise FileNotFoundError(
                f"FASTA file not found: {fasta_file_path}. Cannot initialize genome."
            )
        logger.info(f"Initializing genome from FASTA file: {fasta_file_path}")

        store = self._refget_store_getter()
        metadata, _was_new = store.add_sequence_collection_from_fasta(fasta_file_path)
        digest = metadata.digest
        logger.debug(f"{fasta_file_path} digest: {digest}")

        try:
            self.get(digest)
            created = False
        except MissingGenomeError:
            self.add(
                digest=digest,
                description=description,
                alias_names=alias_names,
                species_name=species_name,
            )
            created = True
        else:
            if not use_existing:
                raise ValueError(
                    f"Genome already exists with digest {digest}. Can't reinitialize genomes."
                )
            logger.debug(f"Genome already exists: {digest}. Using existing.")
            self._repair_missing_aliases(digest, alias_names)

        return digest, created

    def _repair_missing_aliases(self, digest: str, alias_names: list[str]) -> None:
        """Re-add any of ``alias_names`` that no longer resolve to ``digest``.

        Genome rows live in the SQL catalog but aliases live in the RefgetStore,
        so the two can disagree -- and the store's index is rewritten wholesale
        with no locking, so a concurrent writer can drop an alias while leaving
        the genome row untouched.

        Re-init is a no-op when the digest row exists, so dropped aliases must
        be repaired explicitly here.
        """
        if not alias_names or not self._alias_manager_getter:
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
        alias_names: list[str],
        digest: str | None = None,
        use_existing: bool = False,
        species_name: str | None = None,
    ) -> tuple[str, bool]:
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

        try:
            self.get(digest)
            created = False
        except MissingGenomeError:
            self.add(
                digest=digest,
                description=description,
                alias_names=alias_names,
                species_name=species_name,
                remote_url=remote_url,
            )
            created = True
        else:
            if not use_existing:
                raise ValueError(f"Genome already exists with digest {digest}.")
            logger.debug(f"Genome already exists: {digest}. Using existing.")
            # Same repair the local path does: without it, a re-sync can never
            # attach aliases that appeared at the source after first registration.
            self._repair_missing_aliases(digest, alias_names)

        return digest, created

    def _import_remote_collection(self, remote_url: str, digest: str) -> None:
        """Import a collection from a remote store into the local store."""
        import tempfile

        from gtars.refget import RefgetStore as GtarsRefgetStore

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                remote_store = GtarsRefgetStore.open_remote(tmpdir, remote_url)
                collection = remote_store.get_collection(digest)
                local_store = self._refget_store_getter()
                local_store.add_sequence_collection(collection)
                logger.info(f"Imported collection {digest} into local store")
        except Exception as e:
            logger.warning(
                f"Could not import collection {digest} from {remote_url}: {e}. "
                f"Aliases may not resolve until the collection is available locally."
            )

    def export_fasta(self, digest: str, output_path: Path) -> Path:
        """Export a genome's sequences to a FASTA file."""
        genome = self.get(digest)
        store = self._refget_store_getter()

        if genome.remote_url:
            import tempfile

            from gtars.refget import RefgetStore

            try:
                with tempfile.TemporaryDirectory() as tmpdir:
                    remote_store = RefgetStore.open_remote(tmpdir, genome.remote_url)
                    remote_store.export_fasta(digest, str(output_path))
            except (IOError, ConnectionError) as e:
                raise ConnectionError(f"Cannot export FASTA from {genome.remote_url}: {e}") from e
        else:
            store.export_fasta(digest, str(output_path))

        return output_path
