"""Sequence access: locus parsing, sequence retrieval, and remote sequence fallback."""

import re
from collections.abc import Callable
from typing import TYPE_CHECKING, Any, Protocol

from refgenie.exceptions import MissingGenomeError
from refgenie.models import GenomeDigest

if TYPE_CHECKING:
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
        # Route to the store that owns this genome. Local mode (the only mode
        # that serves sequences) has one store; a genome may resolve from the
        # store's alias index with no SQL row yet, so fall back to the default
        # (writable) store when there is no recorded owner.
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

        # Ensure the sequence bytes are available locally before use. The store
        # opens lazily: get_sequence_by_name returns a record with no bytes
        # loaded (healthy --fasta genome), and a metadata-only genome (from
        # `genome init --store`) either returns a bytes-less record or nothing at
        # all. Try to load the bytes from disk; if they are absent locally, fall
        # back to the genome's recorded remote_url. This is what keeps a
        # store-initialized genome usable instead of raising OSError.
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
            # Local mode: check if genome has a remote_url for fallback
            genome = self._genome.get(genome_digest)
            if not genome.remote_url:
                raise ValueError(
                    f"Sequence '{name}' not found in genome '{genome}' (digest: {genome_digest}). "
                    f"If this genome was initialized before RefgetStore integration, "
                    f"reinitialize it with: refgenie genome init --fasta <path> --name {genome}"
                )
            # Fetch from remote and cache locally
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

    def _fetch_remote(self, remote_url: str, collection_digest: str, sequence_name: str):
        """Fetch a single sequence from a remote store and cache it in the default store."""
        from gtars.refget import RefgetStore as GtarsRefgetStore
        import tempfile

        try:
            with tempfile.TemporaryDirectory() as tmpdir:
                remote_store = GtarsRefgetStore.open_remote(tmpdir, remote_url)
                record = remote_store.get_sequence_by_name(collection_digest, sequence_name)

                # Cache ONLY this individual sequence in the global store
                self._router().default_store.add_sequence(record)
        except (IOError, ConnectionError) as e:
            raise ConnectionError(
                f"Cannot fetch sequence '{sequence_name}' from {remote_url}. "
                f"Server may be unreachable: {e}"
            ) from e

        return record
