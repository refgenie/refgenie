"""Tests for parse_registry_path and get_asset_vars."""

from __future__ import annotations

import os

from refgenie.refgenie import get_asset_vars, parse_registry_path


class TestParseRegistryPath:
    """Tests for parse_registry_path."""

    def test_full_path(self):
        result = parse_registry_path("hg38/fasta:default")
        assert result["genome"] == "hg38"
        assert result["asset"] == "fasta"
        assert result["tag"] == "default"
        assert result["seek_key"] is None
        assert result["protocol"] is None

    def test_genome_and_asset_only(self):
        result = parse_registry_path("hg38/fasta")
        assert result["genome"] == "hg38"
        assert result["asset"] == "fasta"
        assert result["tag"] is None

    def test_with_seek_key(self):
        result = parse_registry_path("hg38/fasta.fai:default")
        assert result["genome"] == "hg38"
        assert result["asset"] == "fasta"
        assert result["seek_key"] == "fai"
        assert result["tag"] == "default"

    def test_with_protocol(self):
        result = parse_registry_path("refgenie://hg38/fasta:default")
        assert result["protocol"] == "refgenie"
        assert result["genome"] == "hg38"
        assert result["asset"] == "fasta"
        assert result["tag"] == "default"

    def test_asset_only(self):
        result = parse_registry_path("fasta")
        assert result["genome"] is None
        assert result["asset"] == "fasta"
        assert result["tag"] is None

    def test_returns_dict(self):
        result = parse_registry_path("hg38/fasta")
        assert isinstance(result, dict)
        assert set(result.keys()) == {"protocol", "genome", "asset", "seek_key", "tag"}


class TestGetAssetVars:
    """Tests for get_asset_vars."""

    def test_basic_vars(self, tmp_path):
        outfolder = str(tmp_path)
        result = get_asset_vars("hg38", "fasta", "default", outfolder)
        assert result["genome"] == "hg38"
        assert result["asset"] == "fasta"
        assert result["tag"] == "default"
        assert result["asset_outfolder"] == os.path.join(outfolder, "fasta", "default")

    def test_with_specific_args(self, tmp_path):
        outfolder = str(tmp_path)
        args = {"fasta": "/path/to/file.fa"}
        result = get_asset_vars(
            "hg38", "fasta", "default", outfolder, specific_args=args
        )
        assert result["fasta"] == "/path/to/file.fa"

    def test_with_specific_params(self, tmp_path):
        outfolder = str(tmp_path)
        params = {"threads": "4"}
        result = get_asset_vars(
            "hg38", "fasta", "default", outfolder, specific_params=params
        )
        assert result["threads"] == "4"

    def test_extra_kwargs_passed_through(self, tmp_path):
        outfolder = str(tmp_path)
        result = get_asset_vars(
            "hg38", "fasta", "default", outfolder, custom_key="custom_val"
        )
        assert result["custom_key"] == "custom_val"
