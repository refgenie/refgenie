"""The word rules of the ``?q=`` search: the SQL and Python twins agree, and
description fields match only at word starts."""

import pytest
from sqlalchemy import event
from sqlmodel import Session, select

from refgenie.db.queries import genome as genome_queries
from refgenie.db.queries.search import rank_condition
from refgenie.db.tables import Genome, Species, SpeciesName
from refgenie.utils.search import MatchRank, Search, match_rank
from tests.support.factories import create_asset_class

NAMES = [
    "rat",
    "Norway rat",
    "Rattus norvegicus",
    "mRatBN7.2",
    "naked mole-rat",
    "Carassius auratus",
    "rat_strain (BN)",
    "100% rat",
]


@pytest.mark.parametrize("term", ["rat", "RAT", "mole-rat", "mole rat", "ratt", "BN", "100%"])
def test_sql_ranks_like_python(refgenie_minimal, term):
    """``rank_condition`` puts each name in the tier ``match_rank`` gives it."""
    with Session(refgenie_minimal.database_engine) as session:
        session.add_all(
            Species(taxon_id=i, scientific_name=name, names=[SpeciesName(taxon_id=i, name=name)])
            for i, name in enumerate(NAMES)
        )
        session.commit()
        for rank in (MatchRank.EXACT, MatchRank.WHOLE_WORD, MatchRank.WORD_START):
            condition = rank_condition(SpeciesName.name, term, rank)
            found = set(session.exec(select(SpeciesName.name).where(condition)).all())
            expected = {name for name in NAMES if match_rank(name, term) <= rank}
            assert found == expected, (term, rank)


def test_asset_class_description_matches_at_word_starts(refgenie_minimal):
    with Session(refgenie_minimal.database_engine) as session:
        create_asset_class(session, "bowtie2_index", ["file"])

    def names(term):
        page = refgenie_minimal.asset_class.page(
            search=Search(term=term, fields=["description"]), offset=0, limit=10
        )
        return [row.name for row in page.items]

    # The description is "Test bowtie2_index asset class".
    assert names("bowtie2") == ["bowtie2_index"]
    assert names("index") == ["bowtie2_index"]
    assert names("owtie") == []
    assert names("lass") == []


def test_alias_matches_add_no_bind_parameters(refgenie_minimal):
    """Alias matches travel as one parameter, so a short term matching every
    alias cannot pass the driver's parameter limit (65,535 on Postgres)."""
    engine = refgenie_minimal.database_engine
    with Session(engine) as session:
        session.add(Genome(digest="a" * 32, description=""))
        session.commit()

    def parameter_counts(matches):
        counts = []

        def count(conn, cursor, statement, parameters, context, executemany):
            counts.append(len(parameters))

        event.listen(engine, "before_cursor_execute", count)
        try:
            with Session(engine) as session:
                rows, total = genome_queries.list_page(
                    session, search=Search(term="a"), alias_matches=matches, offset=0, limit=10
                )
        finally:
            event.remove(engine, "before_cursor_execute", count)
        assert [row.digest for row in rows] == ["a" * 32]
        assert total == 1
        return counts

    one = {"a" * 32: MatchRank.EXACT}
    # Far more digests than SQLite (32,766) or Postgres allows as parameters.
    many = {f"{i:032d}": MatchRank.WHOLE_WORD for i in range(70_000)} | one
    assert parameter_counts(many) == parameter_counts(one)
