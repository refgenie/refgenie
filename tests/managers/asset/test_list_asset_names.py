"""`list` shows every name of content that carries several (#295)."""

import pytest

from tests.support.factories import ASSET, GENOME, GROUP, stage_copy
from tests.support.snapshots import fasta_asset

pytestmark = pytest.mark.component


def _add_identical_copy(r, name):
    _, rel = stage_copy(r, fasta_asset(r), f"staging_{name}")  # identical bytes
    r.asset.content.add(
        asset_class_name="fasta",
        path=rel,
        asset_group_name=GROUP,
        genome_digest=r.alias.resolve(GENOME),
        asset_name=name,
    )


def test_list_shows_every_name_of_shared_content(refgenie_built):
    r = refgenie_built
    _add_identical_copy(r, "peer")
    digest = r.alias.resolve(GENOME)
    listed, _ = r.asset.list_all(genome_digests=[digest])
    assert f"{GROUP}:{ASSET}" in listed[digest]
    assert f"{GROUP}:peer" in listed[digest]


def test_list_with_seek_keys_covers_every_name(refgenie_built):
    r = refgenie_built
    _add_identical_copy(r, "peer")
    digest = r.alias.resolve(GENOME)
    listed, _ = r.asset.list_all(genome_digests=[digest], include_seek_keys=True)
    assert f"{GROUP}.fasta:{ASSET}" in listed[digest]
    assert f"{GROUP}.fasta:peer" in listed[digest]
