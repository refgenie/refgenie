"""Data channel router for the refgenie server.

A **data channel** is a static bundle — an ``index.yaml`` plus recipe and
asset-class YAML files — that a refgenie client syncs to build assets. The
canonical channel is published by **refgenie-registry** to GitHub Pages in its
native layout (``asset_classes/<name>.yaml``, ``recipes/<name>/recipe.yaml``).

This server re-serves configured channels as a **transparent, structure-
preserving mirror**: a client sees the SAME ``index.yaml`` and the SAME relative
paths whether it points at the upstream Pages URL or at this server. There is no
filename rewriting and no flattening — the channel served here is structurally
identical to the one on Pages (the index is re-serialized, so bytes may differ).

    GET /v4/data_channel/                       -> list configured channels
    GET /v4/data_channel/<channel>/index.yaml   -> that channel's index (verbatim)
    GET /v4/data_channel/<channel>/<path>       -> 302 to the upstream file

``<path>`` is the channel-relative path taken straight from the channel's own
``index.yaml`` (e.g. ``recipes/bwa_index/recipe.yaml``); the channel name lives
in the URL path segment, so multiple channels coexist without mangling filenames.
"""

import os
import yaml
from fastapi import APIRouter, HTTPException, Response
from refgenie.server.schemas import DataChannelConfig, DataChannelPublic
from urllib.parse import urljoin
from pathlib import Path

router = APIRouter()


def load_data_channels_from_yaml(yaml_path: Path) -> DataChannelConfig:
    """Load data channels configuration from a YAML file and validate with pydantic."""
    with open(yaml_path, "r") as f:
        data_channels = yaml.safe_load(f)
    if not data_channels:
        raise ValueError(f"Data channels configuration file {yaml_path} is empty.")
    return DataChannelConfig(**data_channels)


# Determine config path (env var or default)
DATA_CHANNELS_CONFIG_PATH: Path = Path(
    os.environ.get(
        "DATA_CHANNELS_CONFIG_PATH",
        str(Path(__file__).parent.parent / "data_channels.yaml"),
    )
)
DATA_CHANNELS: DataChannelConfig = load_data_channels_from_yaml(DATA_CHANNELS_CONFIG_PATH)


@router.get("/")
async def list_data_channels() -> dict[str, DataChannelPublic]:
    """List all configured data channels (name, protocol, upstream index_url)."""
    return DATA_CHANNELS.model_dump()


@router.get("/{channel_name}/index.yaml")
@router.head("/{channel_name}/index.yaml")
async def get_channel_index(channel_name: str):
    """Serve a channel's ``index.yaml`` with the same structure as upstream
    (re-serialized, so bytes may differ).

    The client fetches this, then requests each listed relative path back through
    this server, which redirects to the upstream file.
    """
    if (channel := DATA_CHANNELS.root.get(channel_name)) is None:
        raise HTTPException(status_code=404, detail=f"Data channel not found: {channel_name}")

    await channel.fetch_index_yaml()
    index = channel.index_data.model_dump(mode="json")
    return Response(
        content=yaml.dump(index, default_flow_style=False, sort_keys=False),
        media_type="application/x-yaml",
    )


@router.get("/{channel_name}/{file_path:path}")
async def redirect_channel_file(channel_name: str, file_path: str):
    """Redirect a channel-relative path to the upstream file.

    ``file_path`` is resolved against the channel's upstream ``index_url`` with
    ``urljoin`` (which replaces the trailing ``index.yaml``), so
    ``recipes/bwa_index/recipe.yaml`` -> ``<upstream>/recipes/bwa_index/recipe.yaml``.
    """
    if (channel := DATA_CHANNELS.root.get(channel_name)) is None:
        raise HTTPException(status_code=404, detail=f"Data channel not found: {channel_name}")

    target = urljoin(channel.index_url, file_path)
    return Response(
        status_code=302,
        headers={"Location": target},
        media_type="application/x-yaml",
    )
