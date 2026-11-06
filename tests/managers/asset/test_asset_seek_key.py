"""
Tests for reading an asset's seek keys (refgenie.managers.asset.seek_key).

``SeekKeyManager`` is reached as ``rgc.asset.seek_key``. It is built here from
the database, a GenomeManager and an AssetGroupManager alone, which proves its
constructor names everything it needs. Binding seek keys during a write is
covered with the write path (test_asset.py, test_asset_content.py).
"""

import pytest

from refgenie.exceptions import MissingAssetError, MissingSeekKeyError
from refgenie.managers.asset.group import AssetGroupManager
from refgenie.managers.asset.seek_key import SeekKeyManager


class TestSeekKeyManager:
    """Seek-key lookups on a built catalog (unit tier)."""

    def test_works_without_an_asset_manager(self, refgenie_session):
        """Built without AssetManager, every read works."""
        r = refgenie_session
        seek_keys = SeekKeyManager(
            r.database_engine,
            AssetGroupManager(r.database_engine, r.genome),
        )
        digest = r.alias.resolve("rCRSd")

        seek_key = seek_keys.get("fasta", "test", "fasta", genome_digest=digest)
        assert seek_key.name == "fasta"
        assert seek_key.asset.digest
        assert seek_keys.exists("fasta", "test", "fasta", genome_digest=digest)
        assert not seek_keys.exists("fasta", "test", "no_such_key", genome_digest=digest)
        assert seek_keys.get_default("fasta", "test", genome_digest=digest) == "fasta"
        assert "fasta" in seek_keys.list_all(digest, "fasta")

    def test_list_all(self, refgenie_session):
        """list_all returns the fasta seek-key names; the default asset
        resolves to the same keys as naming it explicitly."""
        digest = refgenie_session.alias.resolve("rCRSd")
        keys = refgenie_session.asset.seek_key.list_all(digest, "fasta")
        assert keys
        assert all(isinstance(k, str) for k in keys)
        assert "fasta" in keys
        assert keys == refgenie_session.asset.seek_key.list_all(digest, "fasta", asset_name="test")

    def test_missing_asset_or_key_raises(self, refgenie_session):
        """An unknown asset or seek key is reported as such."""
        seek_keys = refgenie_session.asset.seek_key
        digest = refgenie_session.alias.resolve("rCRSd")
        with pytest.raises(MissingSeekKeyError):
            seek_keys.get("fasta", "test", "no_such_key", genome_digest=digest)
        with pytest.raises(MissingAssetError):
            seek_keys.get_default("fasta", "no_such_asset", genome_digest=digest)
