"""summarize_comparison on hand-built level-2 collections."""

from refgenie.utils.seqcol import compare_level2, summarize_comparison


def _level2(names, seqs, lengths=None):
    return {"names": names, "lengths": lengths or [10] * len(names), "sequences": seqs}


def _summary(a, b):
    comparison = compare_level2(a, b, "A" * 32, "B" * 32)
    return summarize_comparison(a, b, comparison, "a", "b")


ROMAN = ["chrI", "chrII", "chrIII", "chrIV", "chrM"]


def test_identical():
    a = _level2(["x", "y"], ["SQ.1", "SQ.2"])
    assert _summary(a, dict(a)) == ("Identical.", {})


def test_all_renamed():
    names = [f"chr{i}" for i in range(1, 17)] + ["chrM"]
    new = [str(i) for i in range(1, 17)] + ["MT"]
    seqs = [f"SQ.{i}" for i in range(17)]
    summary, name_map = _summary(_level2(names, seqs), _level2(new, seqs))
    assert summary == (
        "Same 17 sequences in the same order; all 17 names differ: "
        "chr1->1, chr2->2, chr3->3, ..., chr16->16, chrM->MT."
    )
    assert name_map["chrM"] == "MT"


def test_irregular_rename_is_not_hidden():
    """sacCer3 vs R64-1-1: chrM->Mito breaks the drop-chr pattern, so it is listed."""
    names = ROMAN[:3] + ["chrM"] + [f"chr{i}" for i in "VWXYZ"]
    new = [n.removeprefix("chr") for n in names[:3]] + ["Mito"] + list("VWXYZ")
    seqs = [f"SQ.{i}" for i in range(len(names))]
    summary, _ = _summary(_level2(names, seqs), _level2(new, seqs))
    assert summary.endswith("chrI->I, chrII->II, chrIII->III, ..., chrZ->Z, chrM->Mito.")


def test_some_renamed():
    seqs = ["SQ.1", "SQ.2", "SQ.3"]
    summary, name_map = _summary(_level2(["a", "b", "M"], seqs), _level2(["a", "b", "MT"], seqs))
    assert summary == "Same 3 sequences in the same order; 1 name differs: M->MT."
    assert name_map == {"M": "MT"}


def test_reordered():
    a = _level2(["a", "b"], ["SQ.1", "SQ.2"])
    b = _level2(["b", "a"], ["SQ.2", "SQ.1"])
    summary, name_map = _summary(a, b)
    assert summary == "Same 2 sequences in a different order; same names."
    assert name_map == {}


def test_same_names_and_lengths_but_content_differs():
    names = [f"c{i}" for i in range(25)]
    a = _level2(names, [f"SQ.a{i}" for i in range(25)])
    b = _level2(names, [f"SQ.a{i}" if i < 7 else f"SQ.b{i}" for i in range(25)])
    assert _summary(a, b)[0] == (
        "Same names and lengths, but 18 of 25 sequences differ in content: c7, c8, c9, ..., c24."
    )


def test_same_names_with_other_sequences_warn():
    """R64-1-1 vs EF2: one sequence is shared, but the chromosome names in both
    hold different sequences."""
    names = ["I", "II", "III", "IV", "V"]
    a = _level2(names + ["Mito"], [f"SQ.a{i}" for i in range(5)] + ["SQ.m"])
    b = _level2(names + ["MT"], [f"SQ.b{i}" for i in range(5)] + ["SQ.m"])
    assert _summary(a, b)[0] == (
        "1 sequence is shared (a has 6, b has 6); 1 name differs: Mito->MT. "
        "Warning: 5 names in both hold different sequences: I, II, III, ..., V."
    )
    c = _level2(["I"], ["SQ.c"], [99])
    assert _summary(a, c)[0] == (
        "No sequences in common. Warning: 1 name in both holds different sequences: I."
    )


def test_hg38_vs_primary_shape():
    a = _level2([f"c{i}" for i in range(195)], [f"SQ.{i}" for i in range(195)])
    b = _level2(
        [f"c{i}" for i in range(7)] + [f"p{i}" for i in range(18)],
        [f"SQ.{i}" for i in range(7)] + [f"SQ.p{i}" for i in range(18)],
    )
    assert _summary(a, b)[0] == "7 sequences are shared (a has 195, b has 25); same names."


def test_subset():
    a = _level2(["x", "y", "z"], ["SQ.1", "SQ.2", "SQ.3"])
    b = _level2(["x", "y"], ["SQ.1", "SQ.2"])
    assert _summary(a, b)[0] == "All 2 sequences of b are in a, which has 1 more; same names."
    assert _summary(a, _level2(["x"], ["SQ.1"]))[0] == (
        "The one sequence of b is in a, which has 2 more; same names."
    )


def test_overlap_and_disjoint():
    a = _level2(["x", "y"], ["SQ.1", "SQ.2"])
    b = _level2(["x", "w", "v"], ["SQ.1", "SQ.8", "SQ.9"])
    assert _summary(a, b)[0].startswith("1 sequence is shared (a has 2, b has 3)")
    c = _level2(["p"], ["SQ.7"], [99])
    assert _summary(a, c) == ("No sequences in common.", {})


def test_duplicate_digests_do_not_crash():
    a = _level2(["x", "y"], ["SQ.1", "SQ.1"])
    b = _level2(["x2", "y2"], ["SQ.1", "SQ.1"])
    summary, name_map = _summary(a, b)
    assert name_map == {"x": "x2", "y": "y2"}
    assert summary.startswith("Same 2 sequences")


def test_comparison_is_unchanged():
    """The GA4GH seqcol comparison, as refget 0.13.0 computes it: attributes in
    one side only, and each order outcome (same, different, undefined)."""
    a = {
        "names": ["chr1", "chr2", "chrM"],
        "lengths": [10, 20, 16569],
        "sequences": ["SQ.1", "SQ.2", "SQ.m"],
        "sorted_sequences": ["SQ.1", "SQ.2", "SQ.m"],
    }
    b = {
        "names": ["chr2", "chr1", "MT"],
        "lengths": [20, 20, 10],
        "sequences": ["SQ.1", "SQ.2", "SQ.x"],
        "name_length_pairs": [{"length": 20, "name": "chr2"}],
    }
    assert compare_level2(a, b, "A" * 32, "B" * 32) == {
        "attributes": {
            "a_only": ["sorted_sequences"],
            "b_only": ["name_length_pairs"],
            "a_and_b": ["lengths", "names", "sequences"],
        },
        "array_elements": {
            "a_count": {"lengths": 3, "names": 3, "sequences": 3, "sorted_sequences": 3},
            "b_count": {"lengths": 3, "name_length_pairs": 1, "names": 3, "sequences": 3},
            "a_and_b_count": {"lengths": 2, "names": 2, "sequences": 2},
            "a_and_b_same_order": {"lengths": None, "names": False, "sequences": True},
        },
        "digests": {"a": "A" * 32, "b": "B" * 32},
    }
