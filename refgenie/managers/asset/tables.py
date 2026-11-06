"""
The asset listing table, shared by the local listing and the remote one.

``rgc.asset.table`` (this node's catalog) and ``rgc.servers.assets_table`` (a
subscribed server's catalog) shape their data the same way, so one function
renders both.
"""

from rich.table import Table

from refgenie.logger import logger
from refgenie.models import AssetRegistryPathComponents
from refgenie.utils.tables import SECTION, build_table


def asset_table(
    asset_data: dict[str, list[str]],
    source: str,
    aliases_data: dict[str, str] | None = None,
    include_seek_keys: bool = False,
) -> Table:
    """
    Create a rich Table from asset data.

    Args:
        asset_data: Genome digest -> registry path strings (``group:asset`` or
            ``group.seek_key:asset``).
        source: Where the data came from ("local" or a server URL), for the title.
        aliases_data: Genome digest -> comma-separated alias names.
        include_seek_keys: Whether to add a seek key column.

    Returns:
        Table: A Rich table.
    """
    aliases_data = aliases_data or {}
    title = f"Refgenie assets. Source: {source}"
    if not asset_data:
        logger.warning("No assets found")
        return build_table(title, [], [])

    columns = ["Aliases", "Genome digest", "Asset group", "Asset"]
    if include_seek_keys:
        columns.append("Seek key")

    rows = []
    for genome_digest, asset_strings in asset_data.items():
        # One section per genome, so a long listing stays readable.
        if rows:
            rows.append(SECTION)
        aliases_str = aliases_data.get(genome_digest, "")
        for asset_string in asset_strings:
            components = AssetRegistryPathComponents.parse_registry_path(asset_string)
            row = [aliases_str, genome_digest, components.asset_group, components.asset]
            if include_seek_keys:
                row.append(components.seek_key)
            rows.append(row)
    return build_table(title, columns, rows)
