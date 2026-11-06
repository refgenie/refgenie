"""``/ping``: the localhost bridge's presence probe.

The public SPA probes ``http://localhost:<port>/ping`` to find a running
``refgenie dash``. Separate from ``/service-info`` so the bridge contract can
change on its own (``bridge_version``); features are in ``capabilities``.
Answers in both modes, with ``Cache-Control: no-store``. CORS decides which
origins may read it.
"""

import logging
import os
import uuid
from pathlib import Path

from fastapi import APIRouter, Depends, Request, Response

from refgenie.const import API_VERSION
from refgenie.core import Refgenie
from refgenie.server.const import ALL_VERSIONS
from refgenie.server.dependencies import get_refgenie
from refgenie.server.schemas import PingBridgeInfo, PingCounts, PingResponse

__all__ = ["BRIDGE_VERSION", "get_instance_id", "router"]

logger = logging.getLogger(__name__)

router = APIRouter(tags=["Default"])

#: Bumped only on a breaking change to the /ping document. A field whose
#: meaning changes gets a new name.
BRIDGE_VERSION = 1


def get_instance_id(home_path: "Path | None" = None) -> str:
    """A random UUID4 created once and stored at ``REFGENIE_HOME_PATH/instance_id``.

    Lets a page recognize an instance; not a credential. The env var is read
    at call time so the id follows the current home path.
    """
    if home_path is None:
        home_path = Path(os.environ.get("REFGENIE_HOME_PATH", Path.home() / ".refgenie"))
    marker = home_path / "instance_id"
    try:
        existing = marker.read_text().strip()
        if existing:
            return existing
    except OSError:
        pass
    new_id = str(uuid.uuid4())
    try:
        home_path.mkdir(parents=True, exist_ok=True)
        marker.write_text(new_id + "\n")
    except OSError:
        # Not fatal: the id just will not survive a restart.
        logger.warning(f"Could not persist instance_id at {marker}")
    return new_id


def _counts(rgc: Refgenie) -> "PingCounts | None":
    """Genome and asset counts, or None on any failure."""
    try:
        counts = rgc.database.summary()
        return PingCounts(genomes=counts["genomes"], assets=counts["assets"])
    except Exception:
        return None


@router.get("/ping", response_model=PingResponse, response_model_exclude_none=True)
def ping(
    request: Request, response: Response, rgc: Refgenie = Depends(get_refgenie)
) -> PingResponse:
    """The bridge handshake document."""
    response.headers["Cache-Control"] = "no-store"
    mode: str = request.app.state.mode
    capabilities: dict = request.app.state.capabilities
    security = getattr(request.app.state, "local_security", None)
    bridge_mode = security.bridge_mode if security is not None else "off"
    expose_paths = security.bridge_expose_paths if security is not None else False
    if expose_paths:
        try:
            instance_label = str(rgc.genome_folder)
        except Exception:
            instance_label = "local refgenie"
    else:
        # No filesystem paths: a home path would leak the username.
        instance_label = "local refgenie" if mode == "local" else "refgenie server"
    return PingResponse(
        bridge_version=BRIDGE_VERSION,
        mode=mode,
        refgenie_version=ALL_VERSIONS["version"],
        api_version=API_VERSION,
        instance_id=get_instance_id(),
        instance_label=instance_label,
        bridge_mode=bridge_mode,
        action_header="X-Refgenie-Action",
        capabilities=capabilities,
        bridge=PingBridgeInfo(actions_cross_origin=bridge_mode == "full"),
        counts=_counts(rgc),
    )
