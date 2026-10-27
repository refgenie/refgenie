"""
Tests for GenomeManager.search / search_condition, the one genome search shared
by the REST route and the MCP tool.

Uses `refgenie_minimal` (local, store-backed mode) seeded with two genomes, so
alias matching exercises the mode-selected alias manager.
"""

import pytest

from refgenie.managers.genome import GENOME_SEARCH_FIELDS
from refgenie.utils.search import SearchOperator
from tests.helpers import seed_search_genomes

A = "a" * 32
B = "b" * 32


@pytest.fixture
def rgc(refgenie_minimal):
    return seed_search_genomes(refgenie_minimal)


def _digests(rgc, *args, **kwargs):
    return sorted(g.digest for g in rgc.genome.search(*args, **kwargs))


@pytest.mark.parametrize(
    ["term", "field"],
    [
        ("sapiens", "species_name"),
        ("Human", "common_name"),
        ("Peculiar", "description"),
        ("UCSC", "assembly_source"),
        ("GCA_000001405", "assembly_accession"),
        ("hg38", "aliases"),
        ("aaaa", "digest"),
    ],
)
def test_search_matches_each_field(rgc, term, field):
    assert field in GENOME_SEARCH_FIELDS
    assert _digests(rgc, term) == [A]
    assert _digests(rgc, term, fields=[field]) == [A]


def test_contains_is_case_insensitive(rgc):
    assert _digests(rgc, "SAPIENS") == [A]


def test_equals_is_case_sensitive(rgc):
    assert _digests(rgc, "human", SearchOperator.EQUALS) == []
    assert _digests(rgc, "Human", SearchOperator.EQUALS) == [A]


def test_aliases_field_matches_alias_only(rgc):
    assert _digests(rgc, "hg38", fields=["aliases"]) == [A]
    assert _digests(rgc, "sapiens", fields=["aliases"]) == []


def test_term_shared_by_both_genomes(rgc):
    assert _digests(rgc, "m") == [A, B]


@pytest.mark.parametrize("term", ["", "   "])
def test_empty_term_returns_nothing(rgc, term):
    assert rgc.genome.search(term) == []
    assert rgc.genome.search_condition(term) is None


def test_no_fields_gives_no_condition(rgc):
    assert rgc.genome.search_condition("human", fields=[]) is None


def test_unknown_field_raises(rgc):
    with pytest.raises(ValueError, match="bogus.*Valid fields are"):
        rgc.genome.search("x", fields=["bogus"])
