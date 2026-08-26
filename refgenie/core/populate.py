"""Resolution of ``refgenie://`` registry paths embedded in arbitrary input.

The module-level functions here are the implementation;
:class:`RegistryPathPopulateMixin` pins them as methods on
:class:`refgenie.core.facade.Refgenie`. :mod:`refgenie.populator` is the thin
looper-facing wrapper on top, and imports from here -- not the other way round.
"""

import re
from typing import TYPE_CHECKING

from refgenie.logger import logger

if TYPE_CHECKING:  # pragma: no cover - annotation only; importing would cycle
    from refgenie.core.facade import Refgenie

Populateable = str | list[str, dict[str, str]]


def populate_refgenie_registry_paths_in_string(
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
            resolved_path = rgc.asset._seek_remote_by_components(
                asset_registry_path_components=parsed_registry_path,
                server_urls=server_urls,
            )
        else:
            resolved_path = rgc.asset._seek_by_components(
                asset_registry_path_components=parsed_registry_path,
            )
        if resolved_path is None:
            logger.warning(f"'{registry_path}' refgenie registry path not populated.")
            continue
        input = re.sub(registry_path, resolved_path, input)
    return input


def populate_refgenie_registry_paths(
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
        return populate_refgenie_registry_paths_in_string(rgc, input, remote, server_urls)
    elif isinstance(input, list):
        return [populate_refgenie_registry_paths(rgc, v, remote, server_urls) for v in input]
    elif isinstance(input, dict):
        return {
            k: populate_refgenie_registry_paths(rgc, v, remote, server_urls)
            for k, v in input.items()
        }
    else:
        raise TypeError(f"Invalide input type. Must be str, list, or dict, not {type(input)}")


class RegistryPathPopulateMixin:
    """
    Resolution of ``refgenie://`` registry paths embedded in arbitrary input.

    Mixed into :class:`refgenie.core.facade.Refgenie`. Each method is a
    delegation to the module-level functions above. Two methods share a name
    with the function they call; a class body is not an enclosing scope, so
    those calls resolve to the module-level function, not to the method.
    """

    def populate_refgenie_registry_paths_in_string(
        self,
        input: str,
        remote: bool = False,
        server_urls: list[str] | None = None,
    ) -> str | None:
        """Populate refgenie registry paths in the input string."""
        return populate_refgenie_registry_paths_in_string(self, input, remote, server_urls)

    def populate(
        self,
        input: str | list[str, dict[str, Populateable]],
    ):
        """Populate refgenie registry paths in the input with local paths."""
        return populate_refgenie_registry_paths(self, input)

    def populater(
        self,
        input: str | list[str, dict[str, Populateable]],
        server_urls: list[str] | None = None,
    ):
        """Populate refgenie registry paths in the input with remote paths."""
        return populate_refgenie_registry_paths(self, input, remote=True, server_urls=server_urls)

    def populate_refgenie_registry_paths(
        self,
        input: str | list[str, dict[str, Populateable]],
        remote: bool = False,
        server_urls: list[str] | None = None,
    ):
        """Populate refgenie registry paths in the input, depending on the remote."""
        return populate_refgenie_registry_paths(self, input, remote, server_urls)
