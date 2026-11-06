"""
API response/request models for the refgenie server, covering each API
surface: the v4 JSON API, the GA4GH DRS endpoints, and data-channel config.
"""

import datetime
from enum import Enum
from typing import Annotated, Any

from urllib.parse import urlsplit

from pydantic import AfterValidator, BaseModel, ConfigDict, Field, RootModel, model_validator
from sqlmodel import SQLModel

from refgenie.db.tables import (
    AssetClassPublic,
    AssetPublic,
    GenomePublic,
    Recipe,
    RecipePublic,
    SeekKeyType,
    Species,
)
from refgenie.models import GenomeListing

# v4 JSON API models


class SpeciesStatistics(BaseModel):
    """Statistics for a single species."""

    genomes: int
    asset_classes: int
    assets: int


class DatabaseSummaryResponse(BaseModel):
    """Response model for database summary endpoint."""

    genomes: int
    asset_groups: int
    assets: int


class AliasResponse(BaseModel):
    """Response for resolving a single alias to its genome data."""

    alias: str
    digest: str
    source: str
    collection: dict[str, Any]
    fhr: dict[str, Any] | None = None


class SpeciesResponse(BaseModel):
    """A genome's species, keyed like a ``species.json`` entry so a client can store it as one."""

    taxon_id: int
    scientific_name: str
    common_name: str | None = None
    #: The species' other names (e.g. "brown rat").
    synonyms: list[str] = Field(default_factory=list)

    @classmethod
    def of(cls, species: Species | None) -> "SpeciesResponse | None":
        if species is None:
            return None
        return cls(
            taxon_id=species.taxon_id,
            scientific_name=species.scientific_name,
            common_name=species.common_name,
            synonyms=species.synonyms,
        )


class GenomeResponse(BaseModel):
    """One genome in a listing. ``common_name`` comes from the species list."""

    digest: str
    aliases: list[str] = Field(default_factory=list)
    description: str | None = None
    asset_count: int = 0
    species_name: str | None = None
    # ``assembly_level`` is on the detail response only.
    common_name: str | None = None
    synonyms: list[str] | None = None
    taxon_id: int | None = None
    assembly_source: str | None = None
    assembly_accession: str | None = None
    species: SpeciesResponse | None = None

    @classmethod
    def of(cls, listing: GenomeListing) -> "GenomeResponse":
        genome = listing.genome
        return cls(
            digest=genome.digest,
            description=genome.description,
            aliases=listing.aliases,
            asset_count=listing.asset_count,
            species_name=genome.species_name,
            common_name=listing.common_name,
            synonyms=genome.synonyms,
            taxon_id=genome.taxon_id,
            assembly_source=genome.assembly_source,
            assembly_accession=genome.assembly_accession,
            species=SpeciesResponse.of(listing.species),
        )


class GenomeDetailResponse(GenomePublic):
    """Detail response for one genome: the core columns plus FHR metadata, taxonomy URI and species."""

    common_name: str | None = None

    #: identifiers.org taxonomy URI derived from ``taxon_id``.
    taxon_uri: str | None = None
    #: The full FHR sidecar (``store.get_fhr_metadata``); may be None.
    fhr: dict[str, Any] | None = None
    species: SpeciesResponse | None = None


class SeekKeyResponse(BaseModel):
    """One seek key of an asset.

    ``value`` is a path relative to the asset directory, or the literal value
    for ``string``/``json`` types. ``from_attributes`` is needed because this
    is filled from ORM rows nested in ``AssetResponse``; the parent's config
    does not reach nested models.
    """

    model_config = ConfigDict(from_attributes=True)

    name: str
    value: str
    description: str | None = None
    type: SeekKeyType
    size: int | None = None


class AssetClassSeekKeyResponse(BaseModel):
    """One seek key an asset class declares: a file pattern, or a typed value."""

    model_config = ConfigDict(from_attributes=True)

    name: str
    value: str | None = None
    description: str | None = None
    type: SeekKeyType


class AssetClassResponse(AssetClassPublic):
    """One asset class with its seek keys, enough for a pulling client to register it."""

    seek_keys: list[AssetClassSeekKeyResponse] = []


class RecipeResponse(RecipePublic):
    """A recipe, with the name of the asset class it builds next to its id."""

    output_asset_class_name: str | None = None

    @classmethod
    def of(cls, recipe: Recipe) -> "RecipeResponse":
        return cls(
            **RecipePublic.model_validate(recipe, from_attributes=True).model_dump(),
            output_asset_class_name=recipe.output_asset_class.name
            if recipe.output_asset_class
            else None,
        )


class AssetNameResponse(BaseModel):
    """One name of an asset, with the provenance of the build behind it.

    Provenance is per name: two builds with identical content share one
    ``Asset`` but each has its own name row. All-null means no build recorded.
    ``from_attributes`` as on ``SeekKeyResponse``.
    """

    model_config = ConfigDict(from_attributes=True)

    name: str
    is_default: bool = False
    build_digest: str | None = None
    build_level1: dict[str, str] | None = None
    build_digest_scheme: str | None = None
    build_timestamp: datetime.datetime | None = None
    refgenie_version: str | None = None
    inputs: dict[str, Any] | None = None
    docker_image: str | None = None
    docker_image_digest: str | None = None
    recipe_id: int | None = None


class AssetResponse(AssetPublic):
    """Public API model for an asset.

    Adds fields read from ``Asset`` ORM properties, including the resolved
    ``serving_modes`` (override, then asset class, then default). Kept off
    ``AssetPublic`` so ``serving_modes`` does not become a DB column.
    """

    serving_modes: list[str]
    asset_class_name: str | None = None
    asset_group_name: str | None = None
    genome_digest: str | None = None
    # Every name this content has in its group, each with its build provenance.
    names: list[AssetNameResponse]
    # Whether any of this asset's names is its group's default.
    is_default: bool | None = None
    # Eager-loaded by ``get_asset``.
    seek_keys: list[SeekKeyResponse] | None = None


# Remote catalog browsing (local mode, GET /v1/remote/*)


class RemoteGenome(BaseModel):
    """A genome on a subscribed remote server (as ``ServerManager.list_genomes`` returns it)."""

    server_url: str
    genome_digest: str
    aliases: list[str] = Field(default_factory=list)
    description: str | None = None


class RemoteServer(BaseModel):
    """A subscribed remote server and whether it currently answers."""

    url: str
    subscribed: bool = True
    reachable: bool = True
    error: str | None = None


class RemoteServersResponse(BaseModel):
    """``GET /v1/remote/servers``."""

    servers: list[RemoteServer] = Field(default_factory=list)


class LocalDataChannel(BaseModel):
    """One data channel this local instance syncs recipes and asset classes from."""

    name: str
    type: str
    index_address: str
    description: str | None = None
    credentials_set: bool = False
    #: Always False: nothing verifies who publishes an index.
    trusted: bool = False


class LocalDataChannelsResponse(BaseModel):
    """``GET /v1/remote/data_channels``."""

    channels: list[LocalDataChannel] = Field(default_factory=list)


class RemoteAsset(BaseModel):
    """An asset on a subscribed remote server (as ``ServerManager.list_assets_for_genome`` returns it)."""

    server_url: str
    genome_digest: str
    asset_group_name: str
    asset_name: str
    asset_digest: str | None = None
    archive_digest: str | None = None
    archive_size: int | None = None


# v4 archives listing (GET /archives)


class ArchiveRecord(SQLModel):
    """One downloadable archive.

    ``digest`` is the asset digest, which ``/archives/{digest}/download``
    accepts. The tarball's checksum is ``tarball_digest``.
    """

    digest: str
    asset_digest: str
    tarball_digest: str | None
    size: int | None
    directory_contents: list[str] | None
    build_commands: list[str] | None
    download_count: int


# DRS Models according to GA4GH DRS specification


class Checksum(BaseModel):
    checksum: str = Field(..., description="The hex-string encoded checksum for the data")
    type: str = Field(..., description="The digest method used to create the checksum")


class AccessURL(BaseModel):
    url: str = Field(
        ...,
        description="A fully resolvable URL that can be used to fetch the actual object bytes",
    )
    headers: list[str] | None = Field(
        None, description="An optional list of headers to include in the HTTP request"
    )


class AccessMethod(BaseModel):
    type: str = Field(..., description="Type of the access method")
    access_url: AccessURL | None = Field(
        None, description="An AccessURL that can be used to fetch the object bytes"
    )
    access_id: str | None = Field(
        None, description="An arbitrary string to be passed to the /access method"
    )
    region: str | None = Field(None, description="Name of the region in the cloud service provider")


class ContentsObject(BaseModel):
    name: str = Field(..., description="A name declared by the bundle author")
    id: str | None = Field(None, description="A DRS identifier of a DrsObject")
    drs_uri: list[str] | None = Field(None, description="A list of full DRS identifier URI paths")
    contents: list["ContentsObject"] | None = Field(None, description="Contents for nested bundles")


class DRSModel(BaseModel):
    id: str = Field(..., description="An identifier unique to this DrsObject")
    name: str | None = Field(None, description="A string that can be used to name a DrsObject")
    self_uri: str = Field(..., description="A drs:// hostname-based URI")
    size: int = Field(..., description="For blobs, the blob size in bytes")
    created_time: datetime.datetime = Field(
        ..., description="Timestamp of content creation in RFC3339"
    )
    updated_time: datetime.datetime | None = Field(
        None, description="Timestamp of content update in RFC3339"
    )
    version: str | None = Field(None, description="A string representing a version")
    mime_type: str | None = Field(
        None, description="A string providing the mime-type of the DrsObject"
    )
    checksums: list[Checksum] = Field(..., description="The checksums of the DrsObject")
    access_methods: list[AccessMethod] = Field(..., description="The list of access methods")
    contents: list[ContentsObject] | None = Field(
        None, description="If set, this DrsObject is a bundle"
    )
    description: str | None = Field(None, description="A human readable description")
    aliases: list[str] | None = Field(None, description="A list of strings for finding metadata")


class ServiceType(BaseModel):
    group: str = Field(..., description="Namespace in reverse domain name format")
    artifact: str = Field(..., description="Name of the API/service")
    version: str = Field(..., description="Version of the API/service")


class Organization(BaseModel):
    name: str = Field(..., description="Name of the organization")
    url: str = Field(..., description="URL of the organization")


class ServiceInfo(BaseModel):
    id: str = Field(..., description="Unique ID of this service")
    name: str = Field(..., description="Name of this service")
    type: ServiceType = Field(..., description="Type of this service")
    description: str | None = Field(None, description="Description of this service")
    organization: Organization = Field(..., description="Organization providing this service")
    contactUrl: str | None = Field(None, description="URL for contacting the developers")
    documentationUrl: str | None = Field(None, description="URL for documentation")
    createdAt: datetime.datetime | None = Field(None, description="Timestamp of service creation")
    updatedAt: datetime.datetime | None = Field(None, description="Timestamp of service update")
    environment: str | None = Field(None, description="Environment the service is running in")
    version: str = Field(..., description="Version of this service")


# For recursive model
ContentsObject.model_rebuild()


# Data-channel config models


def validate_data_channel_name(name: str) -> str:
    if "__" in name:
        raise ValueError("Name cannot contain '__'")
    return name


DataChannelName = Annotated[
    str,
    AfterValidator(validate_data_channel_name),
]


class DataChannelProtocol(str, Enum):
    """Supported data-channel protocols."""

    HTTP = "http"
    HTTPS = "https"


class DataChannel(BaseModel):
    """A data channel's public config: name, protocol and ``index.yaml`` URL."""

    name: DataChannelName
    protocol: DataChannelProtocol
    index_url: str

    @model_validator(mode="after")
    def _validate_index_url(self) -> "DataChannel":
        """``index_url`` must be an http(s) URL with a host, matching ``protocol``.

        The redirect guard resolves against it, so a bad URL fails at config load.
        """
        parts = urlsplit(self.index_url)
        if parts.scheme not in ("http", "https") or not parts.hostname:
            raise ValueError(
                f"index_url must be an http(s) URL with a host, got: {self.index_url!r}"
            )
        if parts.scheme != self.protocol.value:
            raise ValueError(
                f"index_url scheme {parts.scheme!r} does not match protocol {self.protocol.value!r}"
            )
        return self


class DataChannelConfig(RootModel):
    """Data channels by name, each mirrored at ``/v4/data_channel/<channel>/``."""

    root: dict[DataChannelName, DataChannel]


# ---------------------------------------------------------------------------
# /ping -- the localhost bridge's presence probe
# ---------------------------------------------------------------------------


class PingCounts(BaseModel):
    """Cheap headline numbers for the connect UI ("12 genomes, 47 assets")."""

    genomes: int
    assets: int


class PingBridgeInfo(BaseModel):
    """Bridge-only flags, kept apart from the shared ``capabilities`` keys."""

    #: True only under ``bridge_mode == "full"``: bridge origins may pull.
    actions_cross_origin: bool = False


class PingResponse(BaseModel):
    """The ``/ping`` document: the localhost bridge's cross-origin handshake."""

    #: Constant; a check against port collisions, not authentication.
    service: str = "refgenie"
    #: Bumped only on a breaking change to this document. Features are in ``capabilities``.
    bridge_version: int
    mode: str
    refgenie_version: str
    api_version: str
    #: A random UUID stored under REFGENIE_HOME_PATH, so a page can recognize
    #: an instance. Not a credential.
    instance_id: str
    instance_label: str
    bridge_mode: str
    #: The header state-changing requests must carry (presence, any value).
    action_header: str
    #: The same keys /service-info emits. A missing key reads as false.
    capabilities: dict[str, bool]
    bridge: PingBridgeInfo
    #: None when counting is slow or the database is unreachable.
    counts: PingCounts | None = None
