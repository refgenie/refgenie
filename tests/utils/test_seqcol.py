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
        "chr1->1, chr2->2, chr3->3, ..., chrM->MT."
    )
    assert name_map["chrM"] == "MT"


def test_some_renamed():
    seqs = ["SQ.1", "SQ.2", "SQ.3"]
    summary, name_map = _summary(_level2(["a", "b", "M"], seqs), _level2(["a", "b", "MT"], seqs))
    assert summary == "Same 3 sequences in the same order; 1 names differ: M->MT."
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
    assert (
        _summary(a, b)[0] == "Same 25 names and lengths, but 18 of 25 sequences differ in content."
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


def test_overlap_and_disjoint():
    a = _level2(["x", "y"], ["SQ.1", "SQ.2"])
    b = _level2(["x", "w", "v"], ["SQ.1", "SQ.8", "SQ.9"])
    assert _summary(a, b)[0].startswith("1 sequences are shared (a has 2, b has 3)")
    c = _level2(["p"], ["SQ.7"], [99])
    assert _summary(a, c) == ("No sequences in common.", {})


def test_duplicate_digests_do_not_crash():
    a = _level2(["x", "y"], ["SQ.1", "SQ.1"])
    b = _level2(["x2", "y2"], ["SQ.1", "SQ.1"])
    summary, name_map = _summary(a, b)
    assert name_map == {"x": "x2", "y": "y2"}
    assert summary.startswith("Same 2 sequences")
