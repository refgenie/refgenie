"""Remote catalog browsing for local mode (``/v1/remote/*``).

Thin wrappers over ``rgc.servers``; do not add a separate HTTP client here.
A ``server_url`` must be a subscription (``server_not_subscribed``
otherwise): these GET routes are readable by bridge origins, so an arbitrary
URL would let a public page aim this machine at any host. Errors return a
generic message and log the cause.
"""

import logging

from fastapi import APIRouter, Depends, HTTPException, Query

from refgenie.core import Refgenie
from refgenie.models import GenomeDigest
from refgenie.server.dependencies import get_refgenie
from refgenie.server.errors import ErrorCode
from refgenie.server.schemas import (
    LocalDataChannelsResponse,
    RemoteAsset,
    RemoteGenome,
    RemoteServersResponse,
)

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/remote", tags=["Remote"])


def _server_urls(rgc: Refgenie, server_url: str | None) -> list[str] | None:
    """One subscribed server, or None meaning every subscription. 404 if not subscribed."""
    if not server_url:
        return None
    subscribed = rgc.servers.find_subscription(server_url)
    if subscribed is None:
        raise HTTPException(
            status_code=404,
            detail={
                "code": str(ErrorCode.SERVER_NOT_SUBSCRIBED),
                "message": "That server is not one of this refgenie's subscriptions.",
            },
        )
    return [subscribed]


def _remote_unavailable() -> HTTPException:
    """A generic 502; the real cause is logged, not returned."""
    return HTTPException(
        status_code=502,
        detail={
            "code": str(ErrorCode.REMOTE_UNAVAILABLE),
            "message": "Could not query the remote server(s). See the refgenie log for details.",
        },
    )


@router.get("/servers", response_model=RemoteServersResponse)
def list_servers(rgc: Refgenie = Depends(get_refgenie)):
    """The subscribed servers, each checked by fetching its ``/openapi.json``.

    An unreachable server reports ``error: "unreachable"``; details go to the log.
    """
    servers = []
    for url in rgc.servers.subscriptions():
        error = None
        try:
            rgc.servers.client(url).openapi_spec
        except Exception as exc:
            error = "unreachable"
            logger.warning(f"Subscribed server {url} is unreachable: {exc}")
        servers.append({"url": url, "subscribed": True, "reachable": error is None, "error": error})
    return {"servers": servers}


@router.get("/data_channels", response_model=LocalDataChannelsResponse)
def list_data_channels(rgc: Refgenie = Depends(get_refgenie)):
    """The data channels this instance syncs recipes and asset classes from (not probed)."""
    return {
        "channels": [
            {
                "name": channel.name,
                "type": channel.type.value,
                "index_address": channel.index_address,
                "description": channel.description,
                "credentials_set": channel.encrypted_credentials is not None,
                "trusted": False,
            }
            for channel in rgc.sources.list_channels()
        ]
    }


@router.get("/genomes", response_model=list[RemoteGenome])
def list_server_genomes(
    server_url: str | None = Query(
        None, description="Query only this server. Defaults to every subscription."
    ),
    rgc: Refgenie = Depends(get_refgenie),
):
    """Genomes available on the subscribed remote servers."""
    server_urls = _server_urls(rgc, server_url)
    try:
        return rgc.servers.list_genomes(server_urls=server_urls)
    except Exception:
        logger.exception("Failed to list remote genomes")
        raise _remote_unavailable()


@router.get("/assets", response_model=list[RemoteAsset])
def list_server_assets(
    genome_digest: str = Query(..., description="The genome to list assets for."),
    server_url: str | None = Query(
        None, description="Query only this server. Defaults to every subscription."
    ),
    rgc: Refgenie = Depends(get_refgenie),
):
    """Assets available for one genome on the subscribed remote servers."""
    server_urls = _server_urls(rgc, server_url)
    try:
        return rgc.servers.list_assets_for_genome(
            GenomeDigest(genome_digest), server_urls=server_urls
        )
    except Exception:
        logger.exception(f"Failed to list remote assets for {genome_digest}")
        raise _remote_unavailable()
