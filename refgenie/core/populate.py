"""Resolution of ``refgenie://`` registry paths embedded in arbitrary input.

``Refgenie.populate`` and ``Refgenie.populater`` call :func:`populate_registry_paths`.
:mod:`refgenie.integrations.looper` is the thin looper-facing wrapper on top of
``Refgenie.populate``.
"""

import re
from typing import TYPE_CHECKING

from refgenie.exceptions import MissingAssetError, MissingAssetGroupError, MissingSeekKeyError
from refgenie.logger import logger

if TYPE_CHECKING:  # pragma: no cover - annotation only; importing would cycle
    from refgenie.core.root import Refgenie

Populateable = str | list[str, dict[str, str]]


_REGISTRY_PATH = re.compile(r"refgenie://([A-Za-z0-9_\-/\.\:]+)?")


def _populate_string(
    rgc: "Refgenie",
    input: str,
    remote: bool = False,
    server_urls: list[str] | None = None,
) -> str:
    """
    Populate refgenie registry paths in the input string.

    A path that cannot be resolved (a missing asset group, asset or seek key,
    or no server that serves it) only warns and is left as written; an
    unknown genome raises locally.

    Args:
        rgc: The Refgenie instance to resolve paths against.
        input: The input string.
        remote: If True, paths are resolved to remote URLs instead of local paths.
        server_urls: Optional list of server URLs to query for remote resolution.
            If not provided, uses all subscribed servers.

    Returns:
        str: The populated input string.
    """

    def resolve(match: re.Match) -> str:
        registry_path = match.group()
        components = rgc.parse_asset_registry_path(registry_path)
        try:
            if remote:
                return rgc.servers.seek_components(components, server_urls=server_urls)
            return rgc.asset.seek_components(asset_registry_path_components=components)
        except (MissingAssetGroupError, MissingAssetError, MissingSeekKeyError, ValueError) as e:
            logger.warning(f"'{registry_path}' refgenie registry path not populated: {e}")
            return registry_path

    return _REGISTRY_PATH.sub(resolve, input)


def populate_registry_paths(
    rgc: "Refgenie",
    input: str | list[str, dict[str, Populateable]],
    remote: bool = False,
    server_urls: list[str] | None = None,
):
    """
    Populate refgenie registry paths in the input, depending on the remote.

    Args:
        rgc: The Refgenie instance to resolve paths against.
        input: The input to populate.
        remote: If True, paths are resolved to remote URLs instead of local paths.
        server_urls: Optional list of server URLs to query for remote resolution.

    Returns:
        str | list[str, dict[str, Populateable]]: The populated input.
    """
    if isinstance(input, str):
        return _populate_string(rgc, input, remote, server_urls)
    elif isinstance(input, list):
        return [populate_registry_paths(rgc, v, remote, server_urls) for v in input]
    elif isinstance(input, dict):
        return {k: populate_registry_paths(rgc, v, remote, server_urls) for k, v in input.items()}
    elif input is None or isinstance(input, (bool, int, float)):
        # Scalars from a parsed YAML block (a looper pipeline interface, say)
        # cannot hold a registry path; pass them through untouched.
        return input
    else:
        raise TypeError(
            f"Invalid input type. Must be str, list, dict, or a scalar, not {type(input)}"
        )
