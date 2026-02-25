from __future__ import annotations

from refgenconf import CFG_ENV_VARS

__all__ = ["RefgenieError", "MissingGenomeConfigError", "MissingFolderError"]


class RefgenieError(Exception):
    """Base refgenie exception type."""

    pass


class MissingGenomeConfigError(RefgenieError):
    """Exception for when a genome config filepath doesn't point to a file."""

    def __init__(self, conf_file: str | None = None) -> None:
        """Create the error message, using optionally an attempt filepath.

        Args:
            conf_file: Path attempted to be used as genome config file.
        """
        msg = "You must provide a config file either as an argument or via an environment variable: {}".format(
            ", ".join(CFG_ENV_VARS)
        )
        if conf_file:
            msg = "Not a file {} -- {}.".format(conf_file, msg)
        super(MissingGenomeConfigError, self).__init__(msg)


class MissingFolderError(RefgenieError):
    """Exception for when a folder path doesn't exist."""

    def __init__(self, folder: str) -> None:
        """Create the error message.

        Args:
            folder: Path attempted to be used as folder to save a file to.
        """
        super(MissingFolderError, self).__init__(folder)
