"""Sequence access: locus parsing, retrieval, remote fetch and FASTA export."""

import re
import tempfile
from collections.abc import Callable, Iterator
from contextlib import ExitStack, contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Any, Protocol

from refgenie.exceptions import MissingGenomeError, MissingSequenceError, UnservedStoreError
from refgenie.logger import logger
from refgenie.models import GenomeDigest

if TYPE_CHECKING:
    from gtars.refget import RefgetStore

    from refgenie.managers.genome import GenomeManager


_LOCUS_REGEX = re.compile(
    r"^(?P<name>[A-Za-z0-9_.|-]+)(?::(?P<start>\d[\d,]*)(?:-(?P<end>\d[\d,]*))?)?$"
)


def parse_locus(locus: str, *, zero_based: bool = False) -> tuple[str, int | None, int | None]:
    """
    Parse a region written as text into (name, start, end), 0-based half-open.

    Text regions are 1-based and inclusive, as in samtools, UCSC and IGV:
    'chr17:43,044,295-43,044,394' is 100 bases. With ``zero_based`` the text is
    read as 0-based half-open instead (BED style). Commas in numbers are fine.
    Either way the result is 0-based half-open, as the GA4GH refget API and
    :meth:`SequenceManager.get` take it.

    Args:
        locus: e.g. 'chr1', 'chr1:1-1000', 'V01146.1:100-200'.
        zero_based: Read the numbers as 0-based, half-open.

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
    start_str, end_str = match.group("start"), match.group("end")
    start = int(start_str.replace(",", "")) if start_str is not None else None
    end = int(end_str.replace(",", "")) if end_str is not None else None
    if not zero_based and start is not None:
        if start < 1:
            raise ValueError(
                "Text regions are 1-based, so they start at 1; "
                "use --zero-based for 0-based, half-open numbers."
            )
        start -= 1
    if start is not None and end is not None and end <= start:
        raise ValueError(f"The region '{locus}' is empty: its end is not after its start.")
    return name, start, end


@contextmanager
def open_remote_store(url: str) -> Iterator["RefgetStore"]:
    """A remote RefgetStore, opened over a temporary cache directory."""
    from gtars.refget import RefgetStore

    with tempfile.TemporaryDirectory() as tmpdir:
        store = RefgetStore.open_remote(tmpdir, url)
        store.set_quiet(True)
        yield store


class StoreRouter(Protocol):
    """What sequence retrieval needs from the root's store router.

    ``refgenie.core.store_router.RefgetStoreRouter`` satisfies it. Described
    structurally so this manager does not import up into ``refgenie.core``.
    """

    default_store: Any

    def store_for_genome(self, genome_digest: str, store_name: str | None) -> Any: ...


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

    def get(
        self,
        genome_digest: GenomeDigest,
        name: str,
        start: int | None = None,
        end: int | None = None,
    ) -> str:
        """Bases [start, end) of sequence `name` (0-based, half-open).

        None start means 0; None end means the sequence's length. Reads only
        the bytes it needs: from the genome's store when it holds them, else by
        HTTP byte range from the genome's remote_url. Caches nothing.

        Raises:
            MissingSequenceError: If the genome has no sequence called `name`.
            ValueError: If the range is not within the sequence, or a
                metadata-only genome has no remote_url.
            RuntimeError: If the bytes are not local and the mode cannot fetch them.
        """
        router = self._router()
        # The genome's store_name picks its store. A genome may resolve from
        # the local store's alias index with no SQL row, which reads as None.
        genome = None
        try:
            genome = self._genome.get(genome_digest)
        except MissingGenomeError:
            pass
        try:
            store = router.store_for_genome(genome_digest, genome.store_name if genome else None)
        except UnservedStoreError:
            raise
        except MissingGenomeError:
            # Not in the store yet (a metadata-only pull); the remote has it.
            store = router.default_store

        with ExitStack() as stack:
            remote_store = None
            try:
                level2 = store.get_collection_level2(genome_digest)
            except (KeyError, OSError):
                if genome is None or not genome.remote_url or not self._sequences_enabled:
                    raise MissingGenomeError(genome_digest) from None
                remote_store = stack.enter_context(open_remote_store(genome.remote_url))
                level2 = remote_store.get_collection_level2(genome_digest)
            names = level2["names"]
            if name not in names:
                shown = ", ".join(names[:3])
                raise MissingSequenceError(
                    f"Sequence '{name}' is not in genome {genome_digest}; "
                    f"it has {len(names)} sequences, for example {shown}."
                )
            i = names.index(name)
            length = level2["lengths"][i]
            seq_digest = level2["sequences"][i].removeprefix("SQ.")
            start = 0 if start is None else start
            end = length if end is None else end
            if not 0 <= start < end <= length:
                raise ValueError(
                    f"{name} has {length:,} bases; the region {start:,}-{end:,} is outside it."
                )
            if remote_store is None:
                try:
                    return store.get_substring(seq_digest, start, end)
                except (KeyError, OSError):
                    pass  # this store has no bytes for the sequence
            if not self._sequences_enabled:
                raise RuntimeError("Sequence retrieval is not available in server mode")
            if genome is None or not genome.remote_url:
                raise ValueError(
                    f"Sequence '{name}' has no bytes in genome '{genome_digest}' here, and the "
                    f"genome has no remote_url. Reinitialize it with: "
                    f"refgenie genome init --fasta <path> --name <alias>"
                )
            if remote_store is None:
                remote_store = stack.enter_context(open_remote_store(genome.remote_url))
            return remote_store.get_substring(seq_digest, start, end)

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
