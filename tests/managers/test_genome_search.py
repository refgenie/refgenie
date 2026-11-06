"""
Tests for the search in GenomeManager.page, the one genome search shared by
the REST route and the MCP tool.

Uses `refgenie_minimal` (local, store-backed mode) seeded with two genomes, so
alias matching exercises the mode-selected alias manager.
"""

import pytest

from sqlmodel import Session

from refgenie.db.queries.genome import SEARCH_FIELDS
from refgenie.db.tables import Genome
from refgenie.utils.search import Search, SearchOperator
from tests.support.factories import seed_search_genomes, seed_species

A = "a" * 32
B = "b" * 32


@pytest.fixture
def rgc(refgenie_minimal):
    return seed_search_genomes(refgenie_minimal)


def _page(rgc, term, operator=SearchOperator.CONTAINS, fields=None):
    search = Search(term=term, fields=fields, operator=operator)
    return rgc.genome.page(search=search, offset=0, limit=100)


def _digests(rgc, *args, **kwargs):
    return sorted(listing.genome.digest for listing in _page(rgc, *args, **kwargs).items)


@pytest.mark.parametrize(
    ["term", "field"],
    [
        ("sapiens", "species"),
        ("Peculiar", "description"),
        ("UCSC", "assembly_source"),
        ("GCA_000001405", "assembly_accession"),
        ("hg38", "aliases"),
        ("aaaa", "digest"),
    ],
)
def test_search_matches_each_field(rgc, term, field):
    assert field in SEARCH_FIELDS
    assert _digests(rgc, term) == [A]
    assert _digests(rgc, term, fields=[field]) == [A]


def test_contains_is_case_insensitive(rgc):
    assert _digests(rgc, "SAPIENS") == [A]


def test_equals_is_case_sensitive(rgc):
    assert _digests(rgc, "homo sapiens", SearchOperator.EQUALS) == []
    assert _digests(rgc, "Homo sapiens", SearchOperator.EQUALS) == [A]


def test_aliases_field_matches_alias_only(rgc):
    assert _digests(rgc, "hg38", fields=["aliases"]) == [A]
    assert _digests(rgc, "sapiens", fields=["aliases"]) == []


def test_term_shared_by_both_genomes(rgc):
    assert _digests(rgc, "m") == [A, B]


@pytest.mark.parametrize("term", ["", "   "])
def test_empty_term_filters_nothing(rgc, term):
    assert _digests(rgc, term) == [A, B]


@pytest.mark.parametrize(
    ["term", "operator", "found"],
    [
        ("maize", SearchOperator.CONTAINS, True),
        ("aiz", SearchOperator.CONTAINS, True),
        ("maize", SearchOperator.EQUALS, True),
        ("mai", SearchOperator.EQUALS, False),
        ("mai", SearchOperator.STARTS_WITH, True),
        ("MAI", SearchOperator.STARTS_WITH, True),
        ("orn", SearchOperator.STARTS_WITH, False),
        ("orn", SearchOperator.ENDS_WITH, True),
        ("mai", SearchOperator.ENDS_WITH, False),
    ],
)
def test_synonyms_match_each_name_on_its_own(rgc, term, operator, found):
    """The JSON quotes around each synonym bound eq, starts_with and ends_with."""
    maize = "c" * 32
    with Session(rgc.database_engine) as session:
        session.add(Genome(digest=maize, description="", synonyms=["maize", "corn"]))
        session.commit()
    expected = [maize] if found else []
    assert _digests(rgc, term, operator) == expected
    assert _digests(rgc, term, operator, fields=["synonyms"]) == expected


def test_non_ascii_synonym_is_found(rgc):
    """JSON text escapes non-ASCII names; the search escapes the term to match."""
    deer = "d" * 32
    with Session(rgc.database_engine) as session:
        session.add(Genome(digest=deer, description="", synonyms=["Père David's deer"]))
        session.commit()
    assert _digests(rgc, "Père") == [deer]
    assert _digests(rgc, "Père David's deer", SearchOperator.EQUALS) == [deer]


@pytest.fixture
def mouse_rgc(rgc):
    """The two genomes, plus a species-list row for B's taxon (10090).

    ``Hausmaus`` and ``field mouse`` are on no genome column, so a hit on
    either proves the species list was searched.
    """
    seed_species(rgc, 10090, "Mus musculus", "house mouse", ["Hausmaus", "field mouse"])
    return rgc


def _add_genome(rgc, digest, **columns):
    with Session(rgc.database_engine) as session:
        session.add(Genome(digest=digest, description="", **columns))
        session.commit()


@pytest.mark.parametrize("fields", [None, ["species"]], ids=["all-fields", "species"])
@pytest.mark.parametrize(
    ["term", "operator", "found"],
    [
        ("usma", SearchOperator.CONTAINS, True),
        ("field mouse", SearchOperator.CONTAINS, True),
        ("Hausmaus", SearchOperator.EQUALS, True),
        ("hausmaus", SearchOperator.EQUALS, False),
        ("Haus", SearchOperator.EQUALS, False),
        ("HAUS", SearchOperator.STARTS_WITH, True),
        ("maus", SearchOperator.STARTS_WITH, False),
        ("maus", SearchOperator.ENDS_WITH, True),
        ("Haus", SearchOperator.ENDS_WITH, False),
    ],
)
def test_species_list_names_match_each_name_on_its_own(mouse_rgc, term, operator, found, fields):
    assert _digests(mouse_rgc, term, operator, fields=fields) == ([B] if found else [])


def test_genome_without_a_species_row_is_found_by_its_species_name(mouse_rgc):
    assert _digests(mouse_rgc, "sapiens", fields=["species"]) == [A]


def test_genome_without_a_taxon_is_not_found_by_species_list_names(mouse_rgc):
    _add_genome(mouse_rgc, "c" * 32, species_name="Mus spretus")
    assert _digests(mouse_rgc, "Hausmaus") == [B]


def test_common_name_comes_from_the_species_list(mouse_rgc):
    """Every genome of a taxon shows its species row's common name; one without a row has none."""
    mouse = "c" * 32
    _add_genome(mouse_rgc, mouse, species_name="Mus musculus", taxon_id=10090)
    listings = {g.genome.digest: g for g in _page(mouse_rgc, "").items}
    assert listings[mouse].common_name == "house mouse"
    assert listings[B].common_name == "house mouse"
    assert listings[A].common_name is None
    assert mouse_rgc.genome.detail(mouse).common_name == "house mouse"


def test_unknown_field_raises(rgc):
    with pytest.raises(ValueError, match="bogus.*Valid fields are"):
        _page(rgc, "x", fields=["bogus"])
