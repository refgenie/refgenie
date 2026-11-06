"""Sequence access: locus parsing, retrieval, remote fetch and FASTA export."""

import re
import tempfile
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol

from refgenie.exceptions import MissingGenomeError
from refgenie.logger import logger
from refgenie.models import GenomeDigest

if TYPE_CHECKING:
    from gtars.refget import RefgetStore

    from refgenie.managers.genome import GenomeManager


_LOCUS_REGEX = re.compile(r"^(?P<name>[A-Za-z0-9_.|-]+)(?::(?P<start>\d+)(?:-(?P<end>\d+))?)?$")


def parse_locus(locus: str) -> tuple[str, int | None, int | None]:
    """
    Parse a genomic locus string into (name, start, end).

    Coordinates are 0-based, half-open [start, end).

    Args:
        locus: Locus string, e.g. 'chr1', 'chr1:0-1000', 'V01146.1:100-200'.

    Returns:
        Tuple of (name, start, end). start and end are None if not specified.

    Raises:
        ValueError: If the locus format is invalid.
    """
    if not (match := _LOCUS_REGEX.match(locus)):
        raise ValueError(
            f"Invalid locus format: '{locus}'. Expected 'name', 'name:start', or 'name:start-end'."
        )
    name = match.group("name")
    start_str = match.group("start")
    end_str = match.group("end")
    start = int(start_str) if start_str is not None else None
    end = int(end_str) if end_str is not None else None
    return name, start, end


@contextmanager
def open_remote_store(url: str) -> Iterator["RefgetStore"]:
    """A remote RefgetStore, opened over a temporary cache directory."""
    from gtars.refget import RefgetStore

    with tempfile.TemporaryDirectory() as tmpdir:
        yield RefgetStore.open_remote(tmpdir, url)


class StoreRouter(Protocol):
    """What sequence retrieval needs from the root's store router.

    ``refgenie.core.store_router.RefgetStoreRouter`` satisfies it. Described
    structurally so this manager does not import up into ``refgenie.core``.
    """

    names: Any
    default_store: Any

    def get_store(self, name: str) -> Any: ...


class SequenceManager:
    """Sequence retrieval through the RefgetStore router, with remote fallback."""

    def __init__(
        self,
        genome_manager: "GenomeManager",
        store_router_getter: Callable[[], StoreRouter],
        sequences_enabled: bool,
    ):
        """
        Args:
            genome_manager: Genome lookups (store owner, remote_url).
            store_router_getter: Returns the current router. A getter, not a
                router, because ``Refgenie.reload_store_router`` replaces it.
            sequences_enabled: Whether the mode may fetch sequence bytes from a
                remote (local mode only).
        """
        self._genome = genome_manager
        self._router = store_router_getter
        self._sequences_enabled = sequences_enabled

    def get(self, genome_digest: GenomeDigest, locus: str) -> str:
        """
        Return the sequence found in a selected range and chromosome.

        Uses RefgetStore for sequence retrieval (content-addressable, no FASTA file needed).
        For remote genomes, sequences are fetched on demand and cached locally.
        Coordinates are 0-based, half-open [start, end).

        Args:
            genome_digest: The genome digest.
            locus: Coordinates of desired sequence, e.g. 'chr1:0-1000'.
                Format: 'name', 'name:start', or 'name:start-end'.

        Returns:
            str: The requested sequence.

        Raises:
            ValueError: If locus format is invalid or sequence not found.
        """
        name, start, end = parse_locus(locus)
        router = self._router()
        # Route to the store that owns this genome, else the default store
        # (a genome may resolve from the store's alias index with no SQL row).
        store_name = None
        try:
            store_name = self._genome.get(genome_digest).store_name
        except MissingGenomeError:
            pass
        if store_name and store_name in router.names:
            store = router.get_store(store_name)
        else:
            store = router.default_store

        record = None
        try:
            record = store.get_sequence_by_name(genome_digest, name)
        except KeyError:
            record = None

        # Records come back without bytes; load them from disk, else fall back
        # to the genome's remote_url (a metadata-only genome has none locally).
        local_bytes_available = False
        if record is not None:
            try:
                store.load_sequence(record.metadata.sha512t24u)
                # Re-fetch so the record carries the freshly loaded bytes;
                # loading does not mutate an already-fetched record.
                record = store.get_sequence_by_name(genome_digest, name)
                local_bytes_available = True
            except (OSError, IOError):
                local_bytes_available = False

        if not local_bytes_available:
            if not self._sequences_enabled:
                raise RuntimeError("Sequence retrieval is not available in server mode")
            genome = self._genome.get(genome_digest)
            if not genome.remote_url:
                raise ValueError(
                    f"Sequence '{name}' not found in genome '{genome}' (digest: {genome_digest}). "
                    f"If this genome was initialized before RefgetStore integration, "
                    f"reinitialize it with: refgenie genome init --fasta <path> --name {genome}"
                )
            record = self._fetch_remote(genome.remote_url, genome_digest, name)

        if start is not None and end is not None:
            seq_digest = record.metadata.sha512t24u
            store.load_sequence(seq_digest)
            return store.get_substring(seq_digest, start, end)
        elif start is not None:
            # Single position: return from start to end of sequence
            seq_digest = record.metadata.sha512t24u
            store.load_sequence(seq_digest)
            return store.get_substring(seq_digest, start, record.metadata.length)
        else:
            # Whole chromosome
            result = record.decode()
            if result is None:
                raise ValueError(
                    f"Sequence data not available for '{name}' in genome '{genome_digest}'"
                )
            return result

    def fetch(self, digest: str) -> None:
        """Download a remote genome's missing sequences into the local store.

        Remote init imports only the collection's metadata, so anything that
        reads sequence bytes from the local store (the fasta recipe) needs this
        first. A no-op for a genome initialized from a local FASTA.
        """
        remote_url = self._genome.get(digest).remote_url
        if not remote_url:
            return
        local_store = self._router().default_store
        missing = []
        for seq_digest in local_store.get_collection_level2(digest)["sequences"]:
            seq_digest = seq_digest.removeprefix("SQ.")
            try:
                local_store.load_sequence(seq_digest)
            except OSError:
                missing.append(seq_digest)
        if not missing:
            return
        logger.info(f"Downloading sequences for {digest} from {remote_url}")
        with open_remote_store(remote_url) as remote_store:
            for seq_digest in missing:
                # force: the local store already holds a bytes-less stub of this sequence.
                local_store.add_sequence(remote_store.get_sequence(seq_digest), force=True)

    def export_fasta(self, digest: str, output_path: Path) -> Path:
        """Export a genome's sequences to a FASTA file."""
        genome = self._genome.get(digest)
        if not genome.remote_url:
            self._router().default_store.export_fasta(digest, str(output_path))
            return output_path
        try:
            with open_remote_store(genome.remote_url) as remote_store:
                remote_store.export_fasta(digest, str(output_path))
        except (IOError, ConnectionError) as e:
            raise ConnectionError(f"Cannot export FASTA from {genome.remote_url}: {e}") from e
        return output_path

    def _fetch_remote(self, remote_url: str, collection_digest: str, sequence_name: str):
        """Fetch a single sequence from a remote store and cache it in the default store."""
        from gtars.refget import RefgetStore as GtarsRefgetStore
        import tempfile

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                remote_store = GtarsRefgetStore.open_remote(tmpdir, remote_url)
                record = remote_store.get_sequence_by_name(collection_digest, sequence_name)

                # Cache ONLY this sequence; force, because the store holds a bytes-less stub.
                self._router().default_store.add_sequence(record, force=True)
        except (IOError, ConnectionError) as e:
            raise ConnectionError(
                f"Cannot fetch sequence '{sequence_name}' from {remote_url}. "
                f"Server may be unreachable: {e}"
            ) from e

        return record
