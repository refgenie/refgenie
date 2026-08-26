"""``/ping``: the localhost bridge's cross-origin presence probe.

The public SPA (https://refgenie.org) probes ``http://localhost:<port>/ping``
to discover a locally running ``refgenie dash``. This endpoint is deliberately:

* **separate from /service-info** -- that is a GA4GH-shaped discovery document
  with a stable published meaning; this is a browser-handshake contract that
  must be free to evolve with the bridge (D1 in the bridge plan);
* **unversioned and unprefixed** -- it carries its own integer
  ``bridge_version``, bumped only on breaking shape changes; feature
  availability is expressed exclusively through ``capabilities``;
* **present in both modes** -- the same-origin local SPA uses it too, and a
  server-mode answer (``mode: "server"``) tells a misconfigured probe what it
  actually reached;
* **cheap and uncacheable** (``Cache-Control: no-store``) -- detection has a
  ~1.5 s client-side budget.

Whether a cross-origin caller may *read* the answer is entirely the CORS
layer's business (``refgenie/server/local/security.py``); under
``bridge_mode == "off"`` this endpoint still answers for same-origin callers
but no public origin gets a CORS grant.
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

#: Bumped ONLY on a breaking change to the /ping contract's shape or field
#: meanings. Field names are never repurposed: a field whose meaning changes
#: gets a new name and the old one is deleted.
BRIDGE_VERSION = 1


def get_instance_id(home_path: "Path | None" = None) -> str:
    """A random UUID4 generated once and persisted at
    ``REFGENIE_HOME_PATH / "instance_id"``.

    It is **not a secret and not a credential** -- it exists so the page can
    remember "this is the same local refgenie I connected to before". It must
    never be accepted as authorization for anything.

    The env var is re-read at call time (not the import-time constant) so the
    id follows the home path a test or a caller sets.
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
        # An unwritable home dir must not break the probe; the id just will
        # not be stable across restarts.
        logger.warning(f"Could not persist instance_id at {marker}")
    return new_id


def _counts(rgc: Refgenie) -> "PingCounts | None":
    """Two cheap ``SELECT count(*)`` queries; ``None`` (field omitted) on any
    failure rather than a slow or erroring probe."""
    try:
        from sqlmodel import Session, func, select

        from refgenie.db.tables import Asset, Genome

        with Session(rgc.database_engine) as session:
            genomes = session.exec(select(func.count()).select_from(Genome)).one()
            assets = session.exec(select(func.count()).select_from(Asset)).one()
        return PingCounts(genomes=genomes, assets=assets)
    except Exception:
        return None


@router.get("/ping", response_model=PingResponse, response_model_exclude_none=True)
def ping(
    request: Request, response: Response, rgc: Refgenie = Depends(get_refgenie)
) -> PingResponse:
    """The bridge handshake document. See the module docstring."""
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
        # A generic label, no filesystem paths: a home-directory path would
        # leak the OS username to any allowlisted origin (threat T6).
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
