"""Data channel router: re-serves configured channels as a mirror of upstream.

A data channel is an ``index.yaml`` plus recipe and asset-class YAML files.
Clients see the same index and relative paths as upstream:

    GET /v4/data_channel/                       -> list configured channels
    GET /v4/data_channel/<channel>/index.yaml   -> that channel's index
    GET /v4/data_channel/<channel>/<path>       -> 302 to the upstream file
"""

import os
from pathlib import Path
from urllib.parse import quote, unquote, urljoin, urlsplit

import httpx
import yaml
from fastapi import APIRouter, HTTPException, Response

from refgenie import IndexFile
from refgenie.const import HTTP_HEADERS
from refgenie.server.schemas import DataChannel, DataChannelConfig

router = APIRouter()


def load_data_channels_from_yaml(yaml_path: Path) -> DataChannelConfig:
    """Load data channels configuration from a YAML file and validate with pydantic."""
    with open(yaml_path, "r") as f:
        data_channels = yaml.safe_load(f)
    if not data_channels:
        raise ValueError(f"Data channels configuration file {yaml_path} is empty.")
    return DataChannelConfig(**data_channels)


DATA_CHANNELS_CONFIG_PATH: Path = Path(
    os.environ.get(
        "DATA_CHANNELS_CONFIG_PATH",
        str(Path(__file__).parent.parent / "data_channels.yaml"),
    )
)
DATA_CHANNELS: DataChannelConfig = load_data_channels_from_yaml(DATA_CHANNELS_CONFIG_PATH)


@router.get("/")
async def list_data_channels() -> dict[str, DataChannel]:
    """List all configured data channels (name, protocol, upstream index_url)."""
    return DATA_CHANNELS.model_dump()


@router.get("/{channel_name}/index.yaml")
@router.head("/{channel_name}/index.yaml")
async def get_channel_index(channel_name: str):
    """Serve a channel's ``index.yaml`` (re-serialized, so bytes may differ from upstream)."""
    if (channel := DATA_CHANNELS.root.get(channel_name)) is None:
        raise HTTPException(status_code=404, detail=f"Data channel not found: {channel_name}")

    index = (await _fetch_index(channel.index_url)).model_dump(mode="json")
    return Response(
        content=yaml.dump(index, default_flow_style=False, sort_keys=False),
        media_type="application/x-yaml",
    )


async def _fetch_index(index_url: str) -> IndexFile:
    """Fetch and parse a channel's upstream index.yaml; 502 if either fails."""
    try:
        async with httpx.AsyncClient() as client:
            response = await client.get(index_url, headers=HTTP_HEADERS, timeout=10.0)
            response.raise_for_status()
            return IndexFile(**yaml.safe_load(response.text))
    except (httpx.RequestError, httpx.HTTPStatusError):
        raise HTTPException(status_code=502, detail=f"Failed to fetch index from {index_url}")
    except yaml.YAMLError as e:
        raise HTTPException(status_code=502, detail=f"Invalid YAML from {index_url}") from e


def resolve_channel_file_url(index_url: str, file_path: str) -> str | None:
    """Resolve a channel-relative path against the upstream index_url.

    Returns the upstream URL, or None unless ``file_path`` is a plain relative
    path that stays under the index's directory on the same origin (no
    absolute URLs, ``//host``, traversal, query, fragment or control
    characters). Whether the path is listed in ``index.yaml`` is not checked.
    """
    if not file_path:
        return None
    if file_path.startswith("/") or file_path.startswith("\\"):
        return None
    if any(c in file_path for c in ("\\", "?", "#")):
        return None
    if any(ord(c) < 0x20 or ord(c) == 0x7F for c in file_path):
        return None
    if ":" in file_path.split("/", 1)[0]:
        return None
    if any(seg in (".", "..") for seg in unquote(file_path).split("/")):
        return None

    base = urlsplit(index_url)
    # Re-encode: the route receives a decoded path, and a raw space or
    # non-latin-1 character must not go into the Location header as-is.
    target = urljoin(index_url, quote(file_path, safe="/"))
    t = urlsplit(target)

    if t.scheme != base.scheme or t.hostname != base.hostname or t.port != base.port:
        return None

    base_dir = base.path.rsplit("/", 1)[0] + "/"
    if not t.path.startswith(base_dir) or t.path == base_dir:
        return None

    if t.query or t.fragment:
        return None

    return target


@router.get("/{channel_name}/{file_path:path}")
async def redirect_channel_file(channel_name: str, file_path: str):
    """302 to the upstream file for a channel-relative path, or 404 (see ``resolve_channel_file_url``)."""
    if (channel := DATA_CHANNELS.root.get(channel_name)) is None:
        raise HTTPException(status_code=404, detail=f"Data channel not found: {channel_name}")

    if (target := resolve_channel_file_url(channel.index_url, file_path)) is None:
        raise HTTPException(status_code=404, detail="File not found in data channel")

    return Response(
        status_code=302,
        headers={"Location": target},
    )
