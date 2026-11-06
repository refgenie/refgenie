"""SpeciesManager: per-taxon replace, name matching and loading from stores."""

import json
import logging

import pytest

from refgenie.db.tables import Store, StoreType
from refgenie.utils.search import SearchOperator
from refgenie.utils.species_file import SpeciesEntry

RAT = SpeciesEntry(
    taxon_id=10116,
    scientific_name="Rattus norvegicus",
    common_name="Norway rat",
    synonyms=["brown rat", "rat"],
)
MAIZE = SpeciesEntry(
    taxon_id=4577, scientific_name="Zea mays", common_name="maize", synonyms=["corn"]
)
DEER = SpeciesEntry(
    taxon_id=59904, scientific_name="Elaphurus davidianus", common_name="Père David's deer"
)


def test_replace_and_get(refgenie_minimal):
    sp = refgenie_minimal.species
    assert sp.replace([RAT, MAIZE]) == 2
    rat = sp.get(10116)
    assert (rat.scientific_name, rat.common_name) == ("Rattus norvegicus", "Norway rat")
    assert rat.synonyms == ["brown rat", "rat"]
    assert sp.get(1) is None
    assert set(sp.get_many([10116, 4577, 999])) == {10116, 4577}


def test_replace_is_per_taxon(refgenie_minimal):
    sp = refgenie_minimal.species
    sp.replace([RAT, MAIZE])
    sp.replace([RAT.model_copy(update={"synonyms": ["brown rat"]})])
    assert sp.get(10116).synonyms == ["brown rat"]
    assert sp.get(4577).synonyms == ["corn"]
    sp.replace([DEER])  # the pull case: one new taxon, nothing else changes
    assert set(sp.get_many([10116, 4577, 59904])) == {10116, 4577, 59904}
    assert sp.get(4577).synonyms == ["corn"]


@pytest.mark.parametrize(
    "term, operator, expected",
    [
        ("brown", SearchOperator.CONTAINS, [10116]),
        ("MAIZE", SearchOperator.CONTAINS, [4577]),
        ("Norway", SearchOperator.STARTS_WITH, [10116]),
        ("orn", SearchOperator.ENDS_WITH, [4577]),
        ("brown rat", SearchOperator.EQUALS, [10116]),
        ("brown", SearchOperator.EQUALS, []),
        ("Brown rat", SearchOperator.EQUALS, []),
        ("Père", SearchOperator.CONTAINS, [59904]),
        ("  ", SearchOperator.CONTAINS, []),
    ],
)
def test_matching_taxon_ids(refgenie_minimal, term, operator, expected):
    refgenie_minimal.species.replace([RAT, MAIZE, DEER])
    assert refgenie_minimal.species.matching_taxon_ids(term, operator) == expected


def _store(tmp_path, name, priority, species=None, raw=None):
    path = tmp_path / name
    path.mkdir(exist_ok=True)
    if species is not None:
        (path / "species.json").write_text(json.dumps({"species": species}))
    if raw is not None:
        (path / "species.json").write_text(raw)
    return Store(name=name, url=str(path), type=StoreType.on_disk, priority=priority)


def _warnings(caplog):
    return [r for r in caplog.records if r.levelno == logging.WARNING]


class TestLoadFromStores:
    def test_store_without_file_is_fine(self, refgenie_minimal, tmp_path, caplog):
        stores = [_store(tmp_path, "a", 1, [RAT.model_dump()]), _store(tmp_path, "b", 2)]
        with caplog.at_level(logging.DEBUG):
            assert refgenie_minimal.species.load_from_stores(stores) == 1
        assert refgenie_minimal.species.get(10116) is not None
        assert _warnings(caplog) == []

    def test_merges_by_priority_not_argument_order(self, refgenie_minimal, tmp_path):
        high = _store(tmp_path, "high", 1, [{"taxon_id": 1, "scientific_name": "High"}])
        low = _store(tmp_path, "low", 2, [{"taxon_id": 1, "scientific_name": "Low"}])
        refgenie_minimal.species.load_from_stores([low, high])
        sp = refgenie_minimal.species.get(1)
        assert sp.scientific_name == "High"
        assert sp.synonyms == ["Low"]

    def test_unreadable_file_keeps_old_rows(self, refgenie_minimal, tmp_path, caplog):
        good = _store(tmp_path, "good", 1, [RAT.model_dump()])
        bad = _store(tmp_path, "bad", 2, [MAIZE.model_dump()])
        sp = refgenie_minimal.species
        sp.load_from_stores([good, bad])
        _store(tmp_path, "bad", 2, raw="{nope")
        with caplog.at_level(logging.DEBUG):
            assert sp.load_from_stores([good, bad]) is None
        assert len(_warnings(caplog)) == 1
        assert set(sp.get_many([10116, 4577])) == {10116, 4577}

    def test_pulled_rows_survive_and_dropped_synonyms_go(self, refgenie_minimal, tmp_path):
        sp = refgenie_minimal.species
        store = _store(tmp_path, "a", 1, [RAT.model_dump()])
        sp.load_from_stores([store])
        sp.replace([DEER])  # a pulled species no file lists
        _store(tmp_path, "a", 1, [{**RAT.model_dump(), "synonyms": ["brown rat"]}])
        sp.load_from_stores([store])
        assert sp.get(10116).synonyms == ["brown rat"]
        assert sp.get(59904) is not None
