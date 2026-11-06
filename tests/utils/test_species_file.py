"""Reading, validating and merging stores' species.json files."""

import json
import logging

import httpx
import pytest

from refgenie.exceptions import SpeciesFileError
from refgenie.utils.species_file import (
    SpeciesEntry,
    merge_species,
    parse_species,
    read_species_file,
)

RAT = {
    "taxon_id": 10116,
    "scientific_name": "Rattus norvegicus",
    "common_name": "Norway rat",
    "synonyms": ["brown rat", "rat"],
}


def entry(taxon_id, scientific_name, common_name=None, synonyms=None):
    return SpeciesEntry(
        taxon_id=taxon_id,
        scientific_name=scientific_name,
        common_name=common_name,
        synonyms=synonyms,
    )


class TestParse:
    def test_keeps_valid_entries_strips_and_ignores_unknown_keys(self):
        data = {"species": [{**RAT, "scientific_name": "  Rattus norvegicus ", "extra": 1}]}
        (e,) = parse_species(data, "f")
        assert e.scientific_name == "Rattus norvegicus"
        assert e.names() == ["Rattus norvegicus", "Norway rat", "brown rat", "rat"]

    @pytest.mark.parametrize(
        "bad",
        [
            {"scientific_name": "x"},
            {"taxon_id": "10116", "scientific_name": "x"},
            {"taxon_id": True, "scientific_name": "x"},
            {"taxon_id": 0, "scientific_name": "x"},
            {"taxon_id": 5, "scientific_name": "  "},
            {"taxon_id": 5, "scientific_name": "x", "synonyms": [1]},
            "not an object",
        ],
    )
    def test_bad_entries_are_skipped_and_logged(self, bad, caplog):
        with caplog.at_level(logging.DEBUG):
            kept = parse_species({"species": [bad, RAT]}, "f")
        assert [e.taxon_id for e in kept] == [10116]
        assert "skipped species entry 0" in caplog.text

    def test_repeated_taxon_keeps_first(self):
        other = {**RAT, "scientific_name": "Other"}
        (e,) = parse_species({"species": [RAT, other]}, "f")
        assert e.scientific_name == "Rattus norvegicus"

    @pytest.mark.parametrize("data", [[], {}, {"species": {}}])
    def test_bad_top_level_raises(self, data):
        with pytest.raises(SpeciesFileError):
            parse_species(data, "f")


class TestReadOnDisk:
    def test_no_file(self, tmp_path):
        assert read_species_file(str(tmp_path), on_disk=True) is None

    def test_file(self, tmp_path):
        (tmp_path / "species.json").write_text(json.dumps({"species": [RAT]}))
        (e,) = read_species_file(str(tmp_path), on_disk=True)
        assert e.taxon_id == 10116

    def test_malformed_json_raises(self, tmp_path):
        (tmp_path / "species.json").write_text("{nope")
        with pytest.raises(SpeciesFileError):
            read_species_file(str(tmp_path), on_disk=True)


class TestReadHttp:
    @pytest.fixture
    def respond(self, monkeypatch):
        """Install a handler behind every httpx.Client; return the requested URLs."""
        urls = []

        def install(handler):
            def wrapped(request):
                urls.append(str(request.url))
                return handler(request)

            real_client = httpx.Client

            def factory(*args, **kwargs):
                kwargs["transport"] = httpx.MockTransport(wrapped)
                return real_client(*args, **kwargs)

            monkeypatch.setattr(httpx, "Client", factory)
            return urls

        return install

    @pytest.mark.parametrize("root", ["https://x.org/jungle", "https://x.org/jungle/"])
    def test_200(self, respond, root):
        urls = respond(lambda r: httpx.Response(200, json={"species": [RAT]}))
        (e,) = read_species_file(root, on_disk=False)
        assert e.taxon_id == 10116
        assert urls == ["https://x.org/jungle/species.json"]

    @pytest.mark.parametrize("status", [403, 404])
    def test_missing_file(self, respond, status):
        respond(lambda r: httpx.Response(status))
        assert read_species_file("https://x.org/s", on_disk=False) is None

    def test_server_error_raises(self, respond):
        respond(lambda r: httpx.Response(500))
        with pytest.raises(SpeciesFileError, match="HTTP 500"):
            read_species_file("https://x.org/s", on_disk=False)

    def test_connect_error_raises(self, respond):
        def boom(request):
            raise httpx.ConnectError("down")

        respond(boom)
        with pytest.raises(SpeciesFileError):
            read_species_file("https://x.org/s", on_disk=False)

    def test_s3_root_is_not_read(self):
        assert read_species_file("s3://bucket/store", on_disk=False) is None


class TestMerge:
    def test_priority_names_and_union(self):
        high = [entry(1, "Sci A", "Common A", ["alpha"])]
        low = [entry(1, "Sci B", "Common B", ["ALPHA", "beta"])]
        (m,) = merge_species([high, low])
        assert m.scientific_name == "Sci A"
        assert m.common_name == "Common A"
        # case-insensitive dedup keeps the first spelling; the lower store's
        # different scientific and common names stay searchable
        assert m.synonyms == ["alpha", "Sci B", "Common B", "beta"]

    def test_blank_common_name_falls_through(self):
        (m,) = merge_species([[entry(1, "Sci", "  ")], [entry(1, "Sci", "Low common")]])
        assert m.common_name == "Low common"

    def test_sorted_by_taxon_id(self):
        merged = merge_species([[entry(9, "N")], [entry(2, "M")]])
        assert [m.taxon_id for m in merged] == [2, 9]
