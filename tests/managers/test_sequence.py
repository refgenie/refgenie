"""
Tests for the sequence manager (``rgc.sequence``): getseq against a RefgetStore
backend, the remote fallback for metadata-only genomes, and locus parsing.
"""

from pathlib import Path
from unittest.mock import patch

import pytest

from refgenie import Refgenie
from refgenie.exceptions import MissingSequenceError
from refgenie.managers.sequence import parse_locus
from tests.support.factories import make_engine


class TestGenomeSequence:
    """getseq coordinate semantics and error paths (unit tier)."""

    def test_getseq_missing_chromosome_raises(self, refgenie_session):
        """A locus on a chromosome absent from the genome raises.

        Coordinate parsing and substring semantics are covered by the
        RefgetStore tests below.
        """
        with pytest.raises(MissingSequenceError, match="chr999"):
            refgenie_session.sequence.get(refgenie_session.alias.resolve("rCRSd"), "chr999", 0, 10)


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
    """Whole-chromosome getseq loads bytes on demand from a lazily-opened store."""
    r, _digest, name = _local_refgenie_with_bytes(tmp_path, fixtures_path / "rCRSd.fa")

    seq = r.sequence.get(r.alias.resolve("rCRSd"), name)
    assert isinstance(seq, str)
    assert len(seq) == RCRSD_LENGTH


def test_getseq_substring_lazy_store(tmp_path, fixtures_path):
    """Substring getseq works against a lazily-opened store (explicit load path)."""
    r, _digest, name = _local_refgenie_with_bytes(tmp_path, fixtures_path / "rCRSd.fa")

    digest = r.alias.resolve("rCRSd")
    whole = r.sequence.get(digest, name)
    sub = r.sequence.get(digest, name, 0, 10)
    assert len(sub) == 10
    assert whole.startswith(sub)


def test_getseq_metadata_only_routes_to_remote(tmp_path, fixtures_path):
    """getseq on a metadata-only genome falls back to the genome's remote_url.

    The store has the collection but no sequence bytes.
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

    with patch("refgenie.managers.sequence.open_remote_store") as opener:
        opener.return_value.__enter__.return_value = src
        result = r.sequence.get(r.alias.resolve("meta_only"), "rCRSd", 0, 10)

    opener.assert_called_once_with("https://x")
    assert result == src.get_substring(src.get_collection_level2(digest)["sequences"][0][3:], 0, 10)
    # Nothing was cached: the local store still has no bytes for the sequence.
    with pytest.raises(KeyError):
        tgt2 = RefgetStore.on_disk(str(store_path))
        tgt2.get_substring(src.get_collection_level2(digest)["sequences"][0][3:], 0, 10)


def test_getseq_server_mode(tmp_path, fixtures_path):
    """Server mode serves a genome whose own store has the bytes, not a metadata-only one."""
    r, digest, name = _local_refgenie_with_bytes(tmp_path, fixtures_path / "rCRSd.fa")
    r.sequence._sequences_enabled = False
    assert len(r.sequence.get(digest, name, 0, 10)) == 10


def test_getseq_errors(tmp_path, fixtures_path):
    r, digest, name = _local_refgenie_with_bytes(tmp_path, fixtures_path / "rCRSd.fa")
    with pytest.raises(MissingSequenceError, match="chr999.*rCRSd"):
        r.sequence.get(digest, "chr999")
    with pytest.raises(ValueError, match="33,138 bases"):
        r.sequence.get(digest, name, 0, RCRSD_LENGTH + 1)
    with pytest.raises(ValueError):
        r.sequence.get(digest, name, 10, 10)


class TestParseLocus:
    """Tests for the parse_locus helper function."""

    @pytest.mark.parametrize(
        "locus, expected",
        [
            ("chr1", ("chr1", None, None)),  # name only
            ("chr1:1-1000", ("chr1", 0, 1000)),  # 1-based inclusive -> 0-based half-open
            ("chr1:500", ("chr1", 499, None)),  # name with start only
            ("chr17:43,044,295-43,044,394", ("chr17", 43044294, 43044394)),  # commas
            ("V01146.1:1-10", ("V01146.1", 0, 10)),  # dot in name
            ("ENST00000-1:1-10", ("ENST00000-1", 0, 10)),  # hyphen in name
            ("gi|12345|ref|NC_001:1-100", ("gi|12345|ref|NC_001", 0, 100)),  # pipes in name
            ("invalid locus", ValueError),  # space is invalid
            ("", ValueError),  # empty is invalid
        ],
    )
    def test_parse_locus(self, locus, expected):
        """Text regions are 1-based and inclusive, like samtools; the result is
        0-based half-open."""
        if expected is ValueError:
            with pytest.raises(ValueError, match="Invalid locus format"):
                parse_locus(locus)
        else:
            assert parse_locus(locus) == expected

    def test_a_1_based_region_starts_at_1(self):
        with pytest.raises(ValueError, match="start at 1"):
            parse_locus("chr1:0-10")

    def test_an_empty_region_is_refused(self):
        with pytest.raises(ValueError, match="empty"):
            parse_locus("chr1:10-9")

    def test_zero_based(self):
        assert parse_locus("chr1:0-10", zero_based=True) == ("chr1", 0, 10)
        assert parse_locus("chr1:5", zero_based=True) == ("chr1", 5, None)
