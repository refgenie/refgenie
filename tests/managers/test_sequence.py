"""
Tests for the sequence manager (``rgc.sequence``): getseq against a RefgetStore
backend, the remote fallback for metadata-only genomes, and locus parsing.
"""

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from refgenie import Refgenie
from refgenie.managers.sequence import parse_locus
from tests.helpers import make_engine


class TestGenomeSequence:
    """getseq coordinate semantics and error paths (unit tier)."""

    def test_getseq_missing_chromosome_raises(self, refgenie_session):
        """A locus on a chromosome absent from the genome raises.

        Coordinate parsing and substring semantics are covered by the
        RefgetStore tests below.
        """
        with pytest.raises(ValueError):
            refgenie_session.sequence.get(refgenie_session.alias.resolve("rCRSd"), "chr999:0-10")


# ---------------------------------------------------------------------------
# getseq() using the RefgetStore backend (unit tier)
#
# Tests use in-memory RefgetStore to avoid needing samtools or FASTA files on
# disk.
# ---------------------------------------------------------------------------

RCRSD_LENGTH = 33138


def _local_refgenie_with_bytes(tmp_path: Path, fasta_path: Path):
    """A LocalMode Refgenie whose on-disk store has full sequence bytes.

    The store is written to disk *before* Refgenie opens it, so Refgenie opens
    it lazily (no eager sequence load) -- exercising the on-demand load path.

    Returns (refgenie, collection_digest, sequence_name).
    """
    from gtars.refget import RefgetStore

    genome_folder = tmp_path / "genomes"
    store_path = genome_folder / ".refget_store"
    store_path.parent.mkdir(parents=True, exist_ok=True)

    store = RefgetStore.on_disk(str(store_path))
    meta, _ = store.add_sequence_collection_from_fasta(str(fasta_path))
    digest = meta.digest
    store.add_collection_alias("refgenie", "rCRSd", digest)

    r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    r.database.init(genome_folder=genome_folder)
    return r, digest, "rCRSd"


def test_getseq_whole_chromosome_lazy_store(tmp_path, fixtures_path):
    """Whole-chromosome getseq loads bytes on demand from a lazily-opened store.

    With eager load_all_sequences() gone, the record returned by
    get_sequence_by_name has no bytes until explicitly loaded, so record.decode()
    initially returns None. getseq must load the sequence and re-fetch.
    """
    r, _digest, name = _local_refgenie_with_bytes(tmp_path, fixtures_path / "rCRSd.fa")

    seq = r.sequence.get(r.alias.resolve("rCRSd"), name)
    assert isinstance(seq, str)
    assert len(seq) == RCRSD_LENGTH


def test_getseq_substring_lazy_store(tmp_path, fixtures_path):
    """Substring getseq works against a lazily-opened store (explicit load path)."""
    r, _digest, name = _local_refgenie_with_bytes(tmp_path, fixtures_path / "rCRSd.fa")

    digest = r.alias.resolve("rCRSd")
    whole = r.sequence.get(digest, name)
    sub = r.sequence.get(digest, f"{name}:0-10")
    assert len(sub) == 10
    assert whole.startswith(sub)


def test_getseq_metadata_only_routes_to_remote(tmp_path, fixtures_path):
    """getseq on a metadata-only (store-initialized) genome routes to remote fallback.

    A metadata-only collection is present in the store, so get_sequence_by_name
    returns a bytes-less record (NOT a KeyError). The fix must recognize the
    missing local bytes and fall back to the genome's remote_url rather than
    raising OSError -- this is the getseq facet of the store-init poisoning bug.
    """
    from gtars.refget import RefgetStore

    from refgenie.db.tables import Genome, Alias

    genome_folder = tmp_path / "genomes"
    store_path = genome_folder / ".refget_store"
    store_path.parent.mkdir(parents=True, exist_ok=True)

    # Source store with real bytes, used only to produce a valid collection.
    src = RefgetStore.on_disk(str(tmp_path / "src"))
    meta, _ = src.add_sequence_collection_from_fasta(str(fixtures_path / "rCRSd.fa"))
    digest = meta.digest

    # Local store: metadata only (no sequence bytes) -- the store-init state.
    tgt = RefgetStore.on_disk(str(store_path))
    tgt.add_sequence_collection(src.get_collection(digest))
    tgt.add_collection_alias("refgenie", "meta_only", digest)

    r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    r.database.init(genome_folder=genome_folder)

    # Genome row with a remote_url so the fallback has a URL (as store-init records).
    with r.database._database_session as session:
        session.add(Genome(digest=digest, description="meta", remote_url="https://x"))
        session.add(Alias(name="meta_only", genome_digest=digest))
        session.commit()

    mock_record = MagicMock()
    mock_record.decode.return_value = "ACGTACGT"
    mock_record.metadata.sha512t24u = "fake_digest"
    mock_record.metadata.length = 8

    with patch.object(r.sequence, "_fetch_remote", return_value=mock_record) as mock_fetch:
        result = r.sequence.get(r.alias.resolve("meta_only"), "rCRSd")

    # The bytes-less local record must NOT raise OSError; it routes to remote.
    mock_fetch.assert_called_once()
    assert result == "ACGTACGT"


class TestParseLocus:
    """Tests for the parse_locus helper function."""

    @pytest.mark.parametrize(
        "locus, expected",
        [
            ("chr1", ("chr1", None, None)),  # name only
            ("chr1:0-1000", ("chr1", 0, 1000)),  # name with range
            ("chr1:500", ("chr1", 500, None)),  # name with start only
            ("V01146.1:0-10", ("V01146.1", 0, 10)),  # dot in name
            ("ENST00000-1:0-10", ("ENST00000-1", 0, 10)),  # hyphen in name
            ("gi|12345|ref|NC_001:0-100", ("gi|12345|ref|NC_001", 0, 100)),  # pipes in name
            ("invalid locus", ValueError),  # space is invalid
            ("", ValueError),  # empty is invalid
        ],
    )
    def test_parse_locus(self, locus, expected):
        """parse_locus splits name/start/end across name char-classes; rejects bad input."""
        if expected is ValueError:
            with pytest.raises(ValueError, match="Invalid locus format"):
                parse_locus(locus)
        else:
            assert parse_locus(locus) == expected
