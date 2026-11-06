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


def _populate_string(
    rgc: "Refgenie",
    input: str,
    remote: bool = False,
    server_urls: list[str] | None = None,
) -> str | None:
    """
    Populate refgenie registry paths in the input string.

    Args:
        rgc: The Refgenie instance to resolve paths against.
        input: The input string.
        remote: If True, paths are resolved to remote URLs instead of local paths.
        server_urls: Optional list of server URLs to query for remote resolution.
            If not provided, uses all subscribed servers.

    Returns:
        str: The populated input string.
    """
    regex_pattern = re.compile(r"refgenie://([A-Za-z0-9_\-/\.\:]+)?")
    for registry_match in re.finditer(regex_pattern, input):
        registry_path = registry_match.group()
        if not (parsed_registry_path := rgc.parse_asset_registry_path(registry_path)):
            logger.info(f"Can't convert non-conforming refgenie registry path: {registry_path}")
            return input
        if remote:
            if parsed_registry_path.asset is None:
                logger.error("Asset name is required to run remote populate")
                return input
            resolved_path = rgc.servers.seek_components(
                parsed_registry_path, server_urls=server_urls
            )
        else:
            # A missing group, asset or seek key on a known genome only warns
            # (below) and leaves the text as written; an unknown genome raises.
            try:
                resolved_path = rgc.asset.seek_components(
                    asset_registry_path_components=parsed_registry_path,
                )
            except (MissingAssetGroupError, MissingAssetError, MissingSeekKeyError):
                resolved_path = None
        if resolved_path is None:
            logger.warning(f"'{registry_path}' refgenie registry path not populated.")
            continue
        input = re.sub(registry_path, resolved_path, input)
    return input


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
