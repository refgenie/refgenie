"""Tests for CLI entry points (refgenie command)."""

from __future__ import annotations

import subprocess
import sys


def _run_refgenie(*args):
    """Run refgenie as a subprocess and return the result."""
    return subprocess.run(
        [sys.executable, "-m", "refgenie", *args],
        capture_output=True,
        text=True,
    )


class TestCLIEntryPoint:
    """Tests for the refgenie command-line interface."""

    def test_help_exits_zero(self):
        result = _run_refgenie("--help")
        assert result.returncode == 0
        assert "refgenie" in result.stdout

    def test_version_in_help(self):
        result = _run_refgenie("--help")
        assert "version" in result.stdout.lower()

    def test_no_command_exits_nonzero(self):
        result = _run_refgenie()
        assert result.returncode != 0

    def test_init_help(self):
        result = _run_refgenie("init", "--help")
        assert result.returncode == 0
        assert "genome-config" in result.stdout

    def test_build_help(self):
        result = _run_refgenie("build", "--help")
        assert result.returncode == 0
        assert "asset-registry-paths" in result.stdout

    def test_list_help(self):
        result = _run_refgenie("list", "--help")
        assert result.returncode == 0

    def test_init_creates_config(self, tmp_path):
        cfg_path = str(tmp_path / "genome_config.yaml")
        result = _run_refgenie("init", "-c", cfg_path)
        assert result.returncode == 0
        assert (tmp_path / "genome_config.yaml").exists()

    def test_init_config_contains_expected_keys(self, tmp_path):
        cfg_path = str(tmp_path / "genome_config.yaml")
        result = _run_refgenie("init", "-c", cfg_path)
        assert result.returncode == 0
        content = (tmp_path / "genome_config.yaml").read_text()
        assert "genome_folder" in content
        assert "genome_server" in content
        assert "config_version" in content

    def test_list_with_config(self, tmp_path):
        cfg_path = str(tmp_path / "g.yaml")
        _run_refgenie("init", "-c", cfg_path)
        result = _run_refgenie("list", "-c", cfg_path)
        assert result.returncode == 0

    def test_seek_missing_asset_exits_nonzero(self, tmp_path):
        cfg_path = str(tmp_path / "g.yaml")
        _run_refgenie("init", "-c", cfg_path)
        result = _run_refgenie("seek", "-c", cfg_path, "hg38/fasta")
        assert result.returncode != 0
