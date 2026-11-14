"""Collection views declared by the deployment.

A collection is a named, branded slice of the catalog defined by one or more
stores (for example, the Vertebrate Genomes Project assemblies in the ``vgp``
store). ``REFGENIE_COLLECTIONS`` holds a JSON array of :class:`CollectionConfig`
objects; see ``deployment/collections.example.json``. The server parses it
once at startup (a bad value stops startup) and serves the collections whose
stores it serves to the web UI in ``/service-info`` under
``refgenie.collections``.
"""

import os
from typing import Literal

from pydantic import BaseModel, Field, HttpUrl, TypeAdapter

__all__ = [
    "CollectionConfig",
    "CollectionLink",
    "GenomeLinkTemplate",
    "load_collections_from_env",
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


def load_collections_from_env() -> list[CollectionConfig]:
    """Parse ``REFGENIE_COLLECTIONS``; unset gives ``[]``.

    Raises:
        pydantic.ValidationError: If the value is not a JSON array of valid collections.
    """
    raw = os.environ.get("REFGENIE_COLLECTIONS")
    return TypeAdapter(list[CollectionConfig]).validate_json(raw) if raw else []
