"""Collection views declared by the deployment.

A collection is a named, branded slice of the catalog defined by one or more
stores (for example, the Vertebrate Genomes Project assemblies in the ``vgp``
store). ``REFGENIE_COLLECTIONS`` holds a JSON array of :class:`CollectionConfig`
objects; see ``deployment/collections.example.json``. The server validates
them against the stores it loaded and serves the survivors to the web UI in
``/service-info`` under ``refgenie.collections``.
"""

import json
import os
from collections.abc import Iterable
from typing import Literal

from pydantic import BaseModel, Field, HttpUrl, ValidationError

from refgenie.logger import logger

__all__ = [
    "CollectionConfig",
    "CollectionLink",
    "GenomeLinkTemplate",
    "load_collections_from_env",
    "validate_collections",
]

CollectionColumn = Literal["assembly_accession", "haplotype", "assembly_level", "assembly_source"]


class CollectionLink(BaseModel):
    """An external link shown on the collection page."""

    label: str
    url: HttpUrl


class GenomeLinkTemplate(BaseModel):
    """A per-genome external link.

    ``template`` placeholders: ``{assembly_accession}``, ``{species_name}``,
    ``{species_underscored}`` and ``{digest}``. The UI hides the link for a
    genome missing any placeholder it uses.
    """

    label: str
    template: str


class CollectionConfig(BaseModel):
    """One collection view."""

    id: str = Field(pattern=r"^[a-z0-9][a-z0-9-]*$")
    name: str
    #: Used for the placeholder logo badge.
    short_name: str | None = None
    blurb: str
    #: Small-print disclaimer.
    note: str | None = None
    accent_color: str = Field(pattern=r"^#[0-9a-fA-F]{6}$")
    #: ``None`` shows a placeholder badge.
    logo_url: str | None = None
    stores: list[str] = Field(min_length=1)
    links: list[CollectionLink] = []
    genome_links: list[GenomeLinkTemplate] = []
    columns: list[CollectionColumn] = []


def validate_collections(declared: list, store_names: Iterable[str]) -> list[dict]:
    """Validate declared collections against the loaded stores.

    An entry failing validation is logged and skipped; so is one naming a
    store not in ``store_names``, and a duplicate id (the first wins).
    Returns JSON-ready dicts.
    """
    known = set(store_names)
    seen: set[str] = set()
    result: list[dict] = []
    for idx, entry in enumerate(declared):
        try:
            collection = CollectionConfig.model_validate(entry)
        except ValidationError as e:
            logger.error(f"REFGENIE_COLLECTIONS entry {idx} is invalid; skipping it: {e}")
            continue
        missing = [s for s in collection.stores if s not in known]
        if missing:
            logger.warning(
                f"Collection '{collection.id}' names unknown store(s) {', '.join(missing)}; "
                "skipping it."
            )
            continue
        if collection.id in seen:
            logger.warning(f"Duplicate collection id '{collection.id}'; keeping the first.")
            continue
        seen.add(collection.id)
        result.append(collection.model_dump(mode="json"))
    return result


def load_collections_from_env(store_names: Iterable[str]) -> list[dict]:
    """Read and validate ``REFGENIE_COLLECTIONS``; unset or invalid gives ``[]``."""
    raw = os.environ.get("REFGENIE_COLLECTIONS")
    if not raw:
        return []
    try:
        declared = json.loads(raw)
    except json.JSONDecodeError as e:
        logger.error(f"REFGENIE_COLLECTIONS is not valid JSON; ignoring it: {e}")
        return []
    if not isinstance(declared, list):
        logger.error(
            "REFGENIE_COLLECTIONS must be a JSON array of collection objects; ignoring it."
        )
        return []
    return validate_collections(declared, store_names)
