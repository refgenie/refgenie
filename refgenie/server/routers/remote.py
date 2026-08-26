"""Remote catalog browsing for local mode (``/v1/remote/*``).

These endpoints are a thin web skin over the *one* remote client refgenie has:
``rgc.asset.list_remote_genomes`` / ``list_remote_assets_for_genome``, which go
through ``managers/sources/client.py`` (operationId-driven off the remote
server's ``/openapi.json``). The dash's private httpx layer -- its own
pagination, its own error mapping and a TTL cache -- is deleted rather than
ported: a second remote client is exactly the divergence the app factory exists
to remove, and the SPA caches with react-query.

There is no default server. An instance with no subscriptions lists nothing;
the web layer does not invent a server the user never configured.

The handlers are sync ``def``, so FastAPI runs them in its threadpool; blocking
HTTP inside them is fine and costs no event-loop time.
"""

import logging

from fastapi import APIRouter, Depends, HTTPException, Query

from refgenie.core import Refgenie
from refgenie.server.dependencies import get_refgenie
from refgenie.server.schemas import RemoteAsset, RemoteGenome, RemoteServersResponse

logger = logging.getLogger(__name__)

router = APIRouter(prefix="/remote", tags=["Remote"])


def _server_urls(server_url: str | None) -> list[str] | None:
    """One explicit server, or None meaning "every subscription"."""
    return [server_url] if server_url else None


@router.get("/servers", response_model=RemoteServersResponse)
def list_remote_servers(rgc: Refgenie = Depends(get_refgenie)):
    """The servers this instance subscribes to, each with its reachability.

    Reachability is one ``/openapi.json`` fetch per server -- the same document
    the client resolves its operations from, so a server that answers here is a
    server the puller can actually use.
    """
    servers = []
    for url in rgc.configuration.get_server_subscriptions():
        error = None
        try:
            rgc.sources.get_server_client(url).openapi_spec
        except Exception as exc:
            error = str(exc)
            logger.warning(f"Subscribed server {url} is unreachable: {exc}")
        servers.append(
            {"url": url, "subscribed": True, "reachable": error is None, "error": error}
        )
    return {"servers": servers}


@router.get("/genomes", response_model=list[RemoteGenome])
def list_remote_genomes(
    server_url: str | None = Query(
        None, description="Query only this server. Defaults to every subscription."
    ),
    rgc: Refgenie = Depends(get_refgenie),
):
    """Genomes available on the subscribed remote servers."""
    try:
        return rgc.asset.list_remote_genomes(server_urls=_server_urls(server_url))
    except Exception as exc:
        logger.error(f"Failed to list remote genomes: {exc}")
        raise HTTPException(status_code=502, detail=f"Remote server error: {exc}")


@router.get("/assets", response_model=list[RemoteAsset])
def list_remote_assets(
    genome_digest: str = Query(..., description="The genome to list assets for."),
    server_url: str | None = Query(
        None, description="Query only this server. Defaults to every subscription."
    ),
    rgc: Refgenie = Depends(get_refgenie),
):
    """Assets available for one genome on the subscribed remote servers."""
    try:
        return rgc.asset.list_remote_assets_for_genome(
            genome_digest, server_urls=_server_urls(server_url)
        )
    except Exception as exc:
        logger.error(f"Failed to list remote assets for {genome_digest}: {exc}")
        raise HTTPException(status_code=502, detail=f"Remote server error: {exc}")
