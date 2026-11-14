"""display_alias: the one name a person sees for a genome."""

import pytest

from refgenie.models import display_alias

pytestmark = pytest.mark.unit


def test_the_shortest_plain_alias_wins():
    names = ["GRCh38-igenomes-ncbi", "hg38_analysis_set", "GRCh38", "hg38", "jungle::hg"]
    assert display_alias(names) == "hg38"


def test_paths_and_store_qualified_names_are_skipped():
    assert display_alias(["Homo_sapiens/NCBI/GRCh38", "jungle::mm39"]) is None


def test_ties_are_alphabetical():
    assert display_alias(["mm10", "hg19"]) == "hg19"
