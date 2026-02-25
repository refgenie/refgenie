"""Tests for refgenie.helpers module."""

from __future__ import annotations

import os

import pytest

from refgenie.helpers import (
    _parse_user_build_input,
    _single_folder_writeable,
    _writeable,
    make_sure_path_exists,
)


class TestParseUserBuildInput:
    """Tests for _parse_user_build_input."""

    def test_none_input(self):
        result = _parse_user_build_input(None)
        assert result == {}

    def test_empty_list(self):
        result = _parse_user_build_input([])
        assert result == {}

    def test_single_pair(self):
        result = _parse_user_build_input([["fasta=/path/to/file"]])
        assert result == {"fasta": "/path/to/file"}

    def test_multiple_pairs_single_list(self):
        result = _parse_user_build_input([["fasta=/path/to/file", "gtf=/path/to/gtf"]])
        assert result == {"fasta": "/path/to/file", "gtf": "/path/to/gtf"}

    def test_multiple_lists(self):
        result = _parse_user_build_input(
            [["fasta=/path/to/file"], ["gtf=/path/to/gtf"]]
        )
        assert result == {"fasta": "/path/to/file", "gtf": "/path/to/gtf"}

    def test_item_without_equals_ignored(self):
        result = _parse_user_build_input([["no_equals_here", "key=val"]])
        assert result == {"key": "val"}


class TestWriteable:
    """Tests for _writeable and _single_folder_writeable."""

    def test_existing_writable_dir(self, tmp_dir):
        assert _writeable(tmp_dir) is True

    def test_single_folder_writeable(self, tmp_dir):
        assert _single_folder_writeable(tmp_dir) is True

    def test_nonexistent_path_under_writable_parent(self, tmp_dir):
        nonexistent = os.path.join(tmp_dir, "a", "b", "c")
        assert _writeable(nonexistent) is True

    def test_strict_exists_raises_for_missing(self, tmp_dir):
        from refgenie.exceptions import MissingFolderError

        nonexistent = os.path.join(tmp_dir, "nope")
        with pytest.raises(MissingFolderError):
            _writeable(nonexistent, strict_exists=True)

    def test_none_defaults_to_cwd(self):
        assert _writeable(None) is True


class TestMakeSurePathExists:
    """Tests for make_sure_path_exists."""

    def test_creates_nested_dirs(self, tmp_dir):
        nested = os.path.join(tmp_dir, "a", "b", "c")
        make_sure_path_exists(nested)
        assert os.path.isdir(nested)

    def test_existing_dir_no_error(self, tmp_dir):
        make_sure_path_exists(tmp_dir)
        assert os.path.isdir(tmp_dir)
