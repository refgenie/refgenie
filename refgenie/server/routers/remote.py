"""Remote catalog browsing for local mode (``/v1/remote/*``).

These endpoints are a thin web skin over the *one* remote client refgenie has:
``rgc.servers.list_genomes`` / ``list_assets_for_genome``, which go
through ``managers/sources/client.py`` (operationId-driven off the remote
server's ``/openapi.json``). Do not add a private httpx layer here with its
own pagination, error mapping or cache: a second remote client is exactly the
divergence the app factory exists to remove, and the SPA caches with
react-query.

There is no default server. An instance with no subscriptions lists nothing;
the web layer does not invent a server the user never configured.

A ``server_url`` query parameter must name one of the subscriptions, or the
request is refused with ``server_not_subscribed`` before any client is built.
These are GET routes with no action-header guard, and an allowlisted bridge
origin can read them, so an arbitrary URL would let a public page aim this
machine at any host (intranet included) and read the answer. The subscription
list is the allowlist. Failures return a generic message; the cause stays in
the log.

The handlers are sync ``def``, so FastAPI runs them in its threadpool; blocking
HTTP inside them is fine and costs no event-loop time.
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
    """One subscribed server, or None meaning "every subscription".

    A server_url that is not a subscription is refused here, before any
    client is built: these are GET routes with no action-header guard, so
    an arbitrary URL would let any allowlisted page aim this machine at
    any host.
    """
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
    """The generic 502: the real cause is logged, never returned to the caller."""
    return HTTPException(
        status_code=502,
        detail={
            "code": str(ErrorCode.REMOTE_UNAVAILABLE),
            "message": "Could not query the remote server(s). See the refgenie log for details.",
        },
    )


@router.get("/servers", response_model=RemoteServersResponse)
def list_servers(rgc: Refgenie = Depends(get_refgenie)):
    """The servers this instance subscribes to, each with its reachability.

    Reachability is one ``/openapi.json`` fetch per server -- the same document
    the client resolves its operations from, so a server that answers here is a
    server the puller can actually use.

    An unreachable server reports the fixed ``error`` string ``"unreachable"``:
    allowlisted bridge origins can read this route, and the raw exception text
    can carry resolver and host details. The full error goes to the log.
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
    """The data channels this instance syncs recipes and asset classes from.

    Lives next to ``/servers`` because both are the upstream sources the user
    configured; a channel is not probed here, ``sync`` reports reachability.
    """
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
