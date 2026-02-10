"""Shared fixtures for refgenie tests."""

from __future__ import annotations

import gzip
import os

import pytest


@pytest.fixture
def tmp_dir(tmp_path):
    """A temporary directory for test outputs."""
    return str(tmp_path)


@pytest.fixture
def gzipped_file(tmp_path):
    """Create a gzipped file and return its path."""
    p = tmp_path / "test.gz"
    with gzip.open(str(p), "wb") as f:
        f.write(b"test content")
    return str(p)


@pytest.fixture
def plain_file(tmp_path):
    """Create a plain text file and return its path."""
    p = tmp_path / "test.txt"
    p.write_text("test content")
    return str(p)


@pytest.fixture
def test_data_dir():
    """Return the path to the tests/data directory."""
    return os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture
def test_fasta_gz(test_data_dir):
    """Return path to the test t7.fa.gz file."""
    return os.path.join(test_data_dir, "t7.fa.gz")
