"""refgenie's copies of refget 0.13.0's canonical_str and compare_seqcols.

They are copied so a plain install does not need refget. Build digests already
stored in catalogs, and the GA4GH comparison every compare path returns, depend
on them, so they must give exactly what refget gives.
"""

from pathlib import Path

import pytest
from gtars.refget import RefgetStore, sha512t24u_digest

from refgenie.utils.digest import _canonical_str, _canonicalize_value, build_level2
from refgenie.utils.seqcol import _compare_seqcols, seqcol_level2

FASTA = Path(__file__).resolve().parents[1] / "data" / "refget_test_fasta"


@pytest.fixture(scope="module")
def level2s():
    """The level 2 of each refget test FASTA, with every attribute the store holds."""
    store = RefgetStore.in_memory()
    store.set_quiet(True)
    level2s = {}
    for name in ("base", "different_names", "subset"):
        meta, _ = store.add_sequence_collection_from_fasta(str(FASTA / f"{name}.fa"))
        level2s[name] = seqcol_level2(store, meta.digest)
    return level2s


def test_canonical_digest_is_unchanged():
    """Computed with refget 0.13.0."""
    item = {"names": ["chrM"], "lengths": [16569], "sequences": ["SQ.abc"]}
    assert sha512t24u_digest(_canonical_str(item)) == "owW5hNjn8p7wX7-9p4HePlrnO6SWtg0l"


def test_canonical_str_rejects_floats():
    with pytest.raises(ValueError, match="float"):
        _canonical_str({"lengths": [16569.0]})


def test_copies_match_refget(level2s):
    pytest.importorskip("refget")
    from refget.utils import canonical_str, compare_seqcols

    build = _canonicalize_value(
        build_level2(
            genome_digest="genome-A",
            recipe_name="fasta",
            recipe_version="0.0.1",
            input_asset_digests={"parent": "digest-1"},
            input_file_digests={"source": "sha-1"},
            params={"threads": "4"},
        )
    )
    items = [build, {"names": ["chrÄ", "染色体1"]}, {"b": {"z": [], "a": 1}, "a": [-1, 2**40]}]
    for level2 in level2s.values():
        items.append(level2)
        items.extend(level2.values())
    for item in items:
        assert _canonical_str(item) == canonical_str(item)

    for a in level2s.values():
        for b in level2s.values():
            assert _compare_seqcols(a, b) == compare_seqcols(a, b)
