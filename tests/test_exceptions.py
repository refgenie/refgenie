"""Tests for refgenie.exceptions module."""

from __future__ import annotations

import pytest

from refgenie.exceptions import (
    MissingFolderError,
    MissingGenomeConfigError,
    RefgenieError,
)


class TestRefgenieError:
    """Tests for the base exception."""

    def test_is_exception(self):
        assert issubclass(RefgenieError, Exception)

    def test_can_raise(self):
        with pytest.raises(RefgenieError):
            raise RefgenieError("test")


class TestMissingGenomeConfigError:
    """Tests for MissingGenomeConfigError."""

    def test_inherits_refgenie_error(self):
        assert issubclass(MissingGenomeConfigError, RefgenieError)

    def test_message_without_filepath(self):
        err = MissingGenomeConfigError()
        assert "environment variable" in str(err)

    def test_message_with_filepath(self):
        err = MissingGenomeConfigError("/path/to/config.yaml")
        assert "Not a file" in str(err)
        assert "/path/to/config.yaml" in str(err)


class TestMissingFolderError:
    """Tests for MissingFolderError."""

    def test_inherits_refgenie_error(self):
        assert issubclass(MissingFolderError, RefgenieError)

    def test_message_contains_folder_path(self):
        err = MissingFolderError("/some/folder")
        assert "/some/folder" in str(err)
