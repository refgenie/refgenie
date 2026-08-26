"""
API response/request models for the refgenie server, covering each API
surface: the v4 JSON API, the GA4GH DRS endpoints, and data-channel config.
"""

import datetime
from enum import Enum
from typing import Annotated, Any

import httpx
import yaml
from fastapi import HTTPException
from pydantic import AfterValidator, BaseModel, ConfigDict, Field, RootModel
from sqlmodel import SQLModel

from refgenie import IndexFile
from refgenie.const import HTTP_HEADERS
from refgenie.db.tables import AssetPublic, GenomePublic, SeekKeyType

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


class GenomeResponse(BaseModel):
    """
    Response model for listing genomes.
    """

    digest: str
    aliases: list[str] = Field(default_factory=list)
    description: str | None = None
    asset_count: int = 0
    species_name: str | None = None
    # Queryable/faceted core columns. ``assembly_level`` is intentionally
    # omitted from the list response -- it is a detail-page facet, not a scan
    # column.
    common_name: str | None = None
    taxon_id: int | None = None
    assembly_source: str | None = None
    assembly_accession: str | None = None


class GenomeDetailResponse(GenomePublic):
    """Detail response for a single genome.

    Extends the ``GenomePublic`` core columns (digest, description,
    species_name, common_name, taxon_id, assembly_*, remote_url) with the
    read-only FHR sidecar blob and the derived taxonomy URI, so the UI detail
    page gets the whole long tail in one request instead of a second call to
    ``GET /aliases/{name}``.
    """

    #: Derived from ``taxon_id`` (identifiers.org taxonomy URI) so the UI does
    #: not recompute it.
    taxon_uri: str | None = None
    #: The full FHR sidecar (``store.get_fhr_metadata``); may be None.
    fhr: dict[str, Any] | None = None


class SeekKeyResponse(BaseModel):
    """One seek key of an asset, as the asset detail page consumes it.

    ``value`` is the path relative to the asset directory for path-typed seek
    keys, and the literal value for ``string`` / ``json`` types.
    """

    name: str
    value: str
    description: str | None = None
    type: SeekKeyType
    size: int | None = None


class AssetNameResponse(BaseModel):
    """One name of an asset, with the provenance of the build behind it.

    Provenance is per-name, not per-content: two builds can produce identical
    content and share one ``Asset``, and each gets its own name row. All-null
    provenance means the name records no build.

    ``from_attributes`` is required here and not on ``SeekKeyResponse``:
    ``Asset.names`` is a plain property, so it validates through pydantic's
    generic path rather than SQLModel's relationship handling, which supplies
    the setting itself.
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
    """
    Public API response model for an asset.

    Extends AssetPublic with the *resolved* serving_modes (override -> asset_class
    -> default) and the asset class name. These are read directly from the
    ``Asset`` ORM properties of the same name via ``from_attributes``, so the
    server exposes the resolved serving mode rather than only the raw override.

    This is intentionally a non-table subclass: adding ``serving_modes`` as a
    mapped field on ``AssetPublic`` would collide with the ``Asset.serving_modes``
    ORM property and create a spurious DB column.
    """

    serving_modes: list[str] | None = None
    asset_class_name: str | None = None
    # Group name and genome digest, read from the asset's group via the
    # ``Asset.asset_group_name`` / ``Asset.genome_digest`` ORM properties, so
    # clients can render "group / asset" labels and genome links without a
    # per-asset fetch of the asset group.
    asset_group_name: str | None = None
    genome_digest: str | None = None
    # Every name this content has in the asset's group, each with the build
    # behind it. Read from the ``Asset.names`` ORM property (which returns the
    # ``asset_names`` rows), so a client can adopt any of them on pull *and*
    # carry over the build provenance. See AssetName.
    names: list[AssetNameResponse] | None = None
    # Whether this asset is its group's default, read from the
    # ``Asset.is_default`` ORM property (any of its ``asset_names`` rows is
    # flagged). The manage UI's set-default control reflects this.
    is_default: bool | None = None
    # The asset's seek keys (the individual files it exposes). The web UI's
    # asset page lists these, so they must arrive with the asset rather than in
    # a second request; ``get_asset`` eager-loads the relationship.
    seek_keys: list[SeekKeyResponse] | None = None


# Remote catalog browsing (local mode, GET /v1/remote/*)


class RemoteGenome(BaseModel):
    """A genome available on a subscribed remote server.

    Mirrors the dict ``AssetPuller.list_remote_genomes`` already returns; the
    web layer types it for OpenAPI rather than reshaping it.
    """

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


class RemoteAsset(BaseModel):
    """An asset available for a genome on a subscribed remote server.

    Mirrors ``AssetPuller.list_remote_assets_for_genome``.
    """

    server_url: str
    genome_digest: str
    asset_group_name: str
    asset_name: str
    asset_digest: str | None = None
    archive_digest: str | None = None
    archive_size: int | None = None


# v4 archives listing (GET /archives)


class ArchiveRecord(SQLModel):
    """One downloadable archive, as the UI's asset tables consume it.

    ``digest`` is the ASSET digest: it is what ``/archives/{digest}/download``
    accepts, and clients (refgenie-ui Genome/Asset pages) plug this field
    straight into that route. The tarball's own checksum is ``tarball_digest``.
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
    """
    Enum for supported protocols.

    Attributes:
        HTTP: Hypertext Transfer Protocol
        HTTPS: Hypertext Transfer Protocol Secure
    """

    HTTP = "http"
    HTTPS = "https"


class DataChannelPublic(BaseModel):
    """
    DataChannelPublic model to represent a public data channel.

    Attributes:
        protocol: The protocol used (http or https)
        index_url: The URL to the index.yaml file
    """

    name: DataChannelName
    protocol: DataChannelProtocol
    index_url: str


class DataChannel(DataChannelPublic):
    """
    DataChannel model to represent a data channel.

    Attributes:
        name: The name of the data channel, used as a prefix for files. It can't include __.
        protocol: The protocol used (http or https)
        index_url: The URL to the index.yaml file
        index_data: Parsed YAML content from the index.yaml file
    """

    protocol: DataChannelProtocol
    index_url: str
    index_data: IndexFile = None

    async def fetch_index_yaml(self):
        """
        Fetch and parse the channel's index.yaml into ``self.index_data``.

        Raises:
            HTTPException(502) on fetch failure or invalid YAML.
        """
        try:
            async with httpx.AsyncClient() as client:
                response = await client.get(self.index_url, headers=HTTP_HEADERS, timeout=10.0)
                response.raise_for_status()
                self.index_data = IndexFile(**yaml.safe_load(response.text))
        except (httpx.RequestError, httpx.HTTPStatusError):
            raise HTTPException(
                status_code=502, detail=f"Failed to fetch index from {self.index_url}"
            )
        except yaml.YAMLError as e:
            raise HTTPException(
                status_code=502, detail=f"Invalid YAML from {self.index_url}"
            ) from e


class DataChannelConfig(RootModel):
    """
    DataChannelConfig model to represent the mapping of data channels and is a container for methods associated with data channels.

    Channels are re-served per-name at /v4/data_channel/<channel>/ as a
    structure-preserving mirror of their upstream index (see the data_channel
    router). There is no cross-channel aggregation/rewriting.
    """

    root: dict[DataChannelName, DataChannel]


# ---------------------------------------------------------------------------
# /ping -- the localhost bridge's presence probe
# ---------------------------------------------------------------------------


class PingCounts(BaseModel):
    """Cheap headline numbers for the connect UI ("12 genomes, 47 assets")."""

    genomes: int
    assets: int


class PingBridgeInfo(BaseModel):
    """Bridge-only flags, nested so they never mix into ``capabilities``
    (which is the same key set /service-info emits, shared with the SPA)."""

    #: True only under ``bridge_mode == "full"``: whether an allowlisted
    #: public origin may POST /v1/actions/pull.
    actions_cross_origin: bool = False


class PingResponse(BaseModel):
    """The ``/ping`` document: the cross-origin handshake of the localhost
    bridge (see refgenie/server/routers/ping.py). ``/service-info`` is the
    same-origin bootstrap; this is deliberately separate, unversioned (it
    carries its own ``bridge_version``), present in both modes, cheap, and
    served with ``Cache-Control: no-store``.
    """

    #: Constant. Together with ``bridge_version`` this is how a page decides
    #: it is really talking to refgenie -- a sanity check against accidental
    #: port collisions, NOT authentication: any local process can serve this
    #: document, and nothing lets a web page authenticate a loopback peer.
    service: str = "refgenie"
    #: Bumped only on a breaking change to this document's shape or field
    #: meanings. Feature availability is expressed exclusively through
    #: ``capabilities`` -- never gate a UI feature on a version number.
    bridge_version: int
    mode: str
    refgenie_version: str
    api_version: str
    #: A random UUID persisted under REFGENIE_HOME_PATH. Not a secret and not
    #: a credential -- it exists so the page can remember "this is the same
    #: local refgenie I connected to before". Never accepted as authorization.
    instance_id: str
    instance_label: str
    bridge_mode: str
    #: The header state-changing requests must carry (presence, any value).
    action_header: str
    #: The same capability key set /service-info emits -- one shared
    #: vocabulary, no bridge-specific renames. A missing key reads as false.
    capabilities: dict[str, bool]
    bridge: PingBridgeInfo
    #: Omitted (None) when counting is not cheap or the database is not
    #: reachable -- detection has a ~1.5s client budget.
    counts: PingCounts | None = None
