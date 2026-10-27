"""
Tests for asset groups and their default asset (refgenie.managers.asset.group).

``AssetGroupManager`` is reached as ``rgc.asset.group``. It is built here from
the database and a GenomeManager alone, which proves its constructor names
everything it needs. How the write path promotes a default (build, pull,
``content.add(set_default=...)``) lives in test_asset_content.py.
"""

import pytest
from sqlmodel import Session, select

from refgenie.db.tables import AssetName
from refgenie.exceptions import MissingAssetError, MissingAssetGroupError
from refgenie.managers.asset.group import AssetGroupManager
from tests.helpers import build_rcrsd

OLD_NAME = "0.7.17"
NEW_NAME = "0.7.19"


class TestAssetGroupManager:
    """Group lookups and the default read, on a built catalog (unit tier)."""

    def test_works_without_an_asset_manager(self, refgenie_session):
        """Built from the engine and GenomeManager only, every read works."""
        r = refgenie_session
        groups = AssetGroupManager(r.database_engine, r.genome)
        digest = r.alias.resolve("rCRSd")

        assert groups.exists("fasta", genome_digest=digest)
        assert not groups.exists("no_such_group", genome_digest=digest)
        assert groups.get("fasta", genome_digest=digest).name == "fasta"
        assert "fasta" in {g.name for g in groups.list_all(genome_digests=[digest])}
        assert groups.get_default("fasta", genome_digest=digest) == "test"

    def test_missing_group_raises(self, refgenie_session):
        """Neither the group nor its default can be read for an unknown group."""
        digest = refgenie_session.alias.resolve("rCRSd")
        with pytest.raises(MissingAssetGroupError):
            refgenie_session.asset.group.get("no_such_group", genome_digest=digest)
        with pytest.raises(MissingAssetGroupError):
            refgenie_session.asset.group.get_default("no_such_group", genome_digest=digest)

    def test_genome_is_required(self, refgenie_session):
        """A lookup without a genome is a caller error."""
        with pytest.raises(TypeError, match="genome_digest"):
            refgenie_session.asset.group.get("fasta")


class TestAssetGroupDefaults:
    """The one-default-per-group invariant (component tier)."""

    pytestmark = pytest.mark.component

    def test_set_default_records_exact_name_and_keeps_one_flag(self, refgenie_fs):
        """set_default records the name the caller passed (not the canonical name) and
        flipping it clears the previous flag -- at most one default per group."""
        r = refgenie_fs
        build_rcrsd(r, asset_name=OLD_NAME)
        build_rcrsd(r, asset_name=NEW_NAME)
        genome_digest = r.alias.resolve("rCRSd")
        for name in (NEW_NAME, OLD_NAME):
            r.asset.group.set_default(
                genome_digest=genome_digest, asset_group_name="fasta", asset_name=name
            )
            assert r.asset.group.get_default("fasta", genome_digest=genome_digest) == name
            with Session(r.database_engine) as session:
                defaults = session.exec(select(AssetName).where(AssetName.is_default)).all()
            assert [d.name for d in defaults] == [name]

    def test_set_default_from_different_group_raises(self, refgenie_fs):
        """A name that belongs to no group cannot be this group's default."""
        r = refgenie_fs
        build_rcrsd(r, asset_name=OLD_NAME)
        with pytest.raises(MissingAssetError):
            r.asset.group.set_default(
                genome_digest=r.alias.resolve("rCRSd"),
                asset_group_name="fasta",
                asset_name="a-name-from-nowhere",
            )
