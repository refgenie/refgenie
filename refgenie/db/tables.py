"""SQLModel table definitions for the refgenie metadata database.

These models have filesystem side effects, and the rule governing them is:
**the catalog commits first, and the filesystem follows.** ``db/events.py``
registers mapper listeners on ``Asset``, ``Alias`` and ``StagedAsset`` that
resolve which files a delete makes obsolete -- which can only be done while the
row is still there -- and queue them via ``db/cleanup.py``. A session-level
``after_commit`` listener runs the queue; ``after_rollback`` discards it. So
``session.delete(asset); session.commit()`` still removes that asset's content
directory, and deleting a ``Genome`` still cascades into its assets and their
files, but a transaction that rolls back leaves every byte in place alongside
the rows it restored.

No code path may destroy managed data before the catalog agrees it is gone.
Where a manager does filesystem work outside a listener, the same order applies:
commit, then clean up, and make the cleanup safe to re-run.

Every listener and exactly what it touches on disk is tabulated in
``docs/development.md`` under "Filesystem side effects of the ORM"; read that
before writing code that inserts or deletes these rows.

Cascade and route-ownership rationale lives in ``docs/design-notes.md``.
"""

import os
import re
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Optional

from pydantic import model_validator
from sqlalchemy import (
    BigInteger,
    ForeignKeyConstraint,
    Index,
    UniqueConstraint,
    func,
    text,
)
from sqlmodel import JSON, Column, DateTime, Field, Relationship, SQLModel
from yaml import dump as ydump

SA_RELATIONSHIP_KWARGS_BASE = dict(lazy="select")
SA_RELATIONSHIP_KWARGS_ORPHAN = SA_RELATIONSHIP_KWARGS_BASE | dict(cascade="all, delete-orphan")


class SeekKeyType(Enum):
    """
    The seek key type enumeration.
    """

    # Path-based types
    file = "file"
    directory = "directory"
    prefix = "prefix"
    # Non-path types
    string = "string"
    json = "json"


PATH_SEEK_KEY_TYPES = {SeekKeyType.file, SeekKeyType.directory, SeekKeyType.prefix}


class ServingMode(str, Enum):
    """The serving mode for an asset class or asset."""

    file = "file"
    archive = "archive"
    none = "none"


VALID_SERVING_MODES = frozenset(mode.value for mode in ServingMode)


def _check_serving_modes(modes: list[str]) -> None:
    """Raise ValueError if any of ``modes`` is not a valid :class:`ServingMode`."""
    invalid = set(modes) - VALID_SERVING_MODES
    if invalid:
        raise ValueError(
            f"Invalid serving modes: {invalid}. Valid modes: {sorted(VALID_SERVING_MODES)}"
        )


def is_path_type(seek_key_type: SeekKeyType) -> bool:
    """Return True if the seek key type represents a filesystem path."""
    return seek_key_type in PATH_SEEK_KEY_TYPES


class RemoteType(Enum):
    """
    The remote type enumeration.
    """

    s3 = "s3"
    http = "http"
    https = "https"


class StoreType(str, Enum):
    """How a :class:`Store` is opened.

    ``remote`` stores are ``rgstore.json`` manifests reached over http(s)/S3;
    ``on_disk`` stores are local ``.refget_store`` directories.
    """

    remote = "remote"
    on_disk = "on_disk"


class GenomePublic(SQLModel):
    """
    A reference genome model.
    """

    digest: str = Field(default=None, primary_key=True, max_length=32, min_length=32)
    description: str | None
    #: Scientific (species) name, e.g. "Homo sapiens".
    species_name: str | None = Field(default=None, index=True)
    #: Common name, e.g. "human". Promoted from the FHR sidecar so it is
    #: queryable/faceted alongside the scientific name.
    common_name: str | None = None
    #: NCBI taxonomy id, e.g. 9606. Indexed so genomes can be grouped by taxon.
    taxon_id: int | None = Field(default=None, index=True)
    #: Assembly provider (NCBI / UCSC / Ensembl / 1000G). Indexed: a top facet.
    assembly_source: str | None = Field(default=None, index=True)
    #: Assembly accession (GCA_/GCF_...). Indexed so it is directly searchable.
    assembly_accession: str | None = Field(default=None, index=True)
    #: Assembly level (chromosome / scaffold / contig). A cheap second facet.
    assembly_level: str | None = None
    remote_url: str | None = Field(default=None)


class Genome(GenomePublic, table=True):
    """
    A reference genome model.
    """

    #: Name of the :class:`Store` that authoritatively owns this genome's
    #: sequence bytes and collection metadata. Nullable only transiently during
    #: migration/ingestion; the federated router uses it to dispatch reads. Kept
    #: off :class:`GenomePublic` so it stays an internal routing column, not a
    #: wire field.
    store_name: str | None = Field(default=None, foreign_key="store.name", index=True)
    store: Optional["Store"] = Relationship(
        back_populates="genomes", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    aliases: list["Alias"] = Relationship(
        back_populates="genome", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN
    )
    asset_groups: list["AssetGroup"] = Relationship(
        back_populates="genome", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )

    def __str__(self):
        parts = [self.digest]
        if self.species_name:
            parts.append(self.species_name)
        if self.description:
            parts.append(self.description)
        return f"{parts[0]} ({', '.join(parts[1:])})" if len(parts) > 1 else parts[0]


class AliasPublic(SQLModel):
    """
    A reference genome alias model.
    """

    name: str = Field(primary_key=True)
    genome_digest: str = Field(default=None, foreign_key="genome.digest")


class Alias(AliasPublic, table=True):
    """
    A reference genome alias model.
    """

    #: Owning store of the genome this alias points at. Denormalized from
    #: ``genome.store_name`` so collisions and qualified ``store::alias`` lookups
    #: are answerable without re-reading the stores. Kept off ``AliasPublic`` so
    #: it stays an internal routing column, not a wire field.
    store_name: str | None = Field(default=None, foreign_key="store.name", index=True)
    genome: Optional["Genome"] = Relationship(
        back_populates="aliases", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )

    def __str__(self):
        return self.name


class AssetGroupPublic(SQLModel):
    """
    A reference asset group model.
    """

    id: int | None = Field(default=None, primary_key=True)
    name: str
    description: str | None
    genome_digest: str = Field(default=None, foreign_key="genome.digest")
    asset_class_id: int = Field(default=None, foreign_key="assetclass.id")

    def __repr__(self):
        return f"{self.genome_digest}/{self.name}"


class AssetGroup(AssetGroupPublic, table=True):
    """
    A reference asset group model.
    """

    genome_digest: str = Field(default=None, foreign_key="genome.digest", index=True)
    genome: "Genome" = Relationship(
        back_populates="asset_groups",
        sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE,
    )
    assets: list["Asset"] = Relationship(
        back_populates="asset_group",
        sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN,
    )
    asset_class: Optional["AssetClass"] = Relationship(back_populates="asset_groups")
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )

    # A genome cannot have two asset groups with the same name. AssetName's
    # scoped uniqueness (asset_group_id, name) is only as strong as this.
    __table_args__ = (UniqueConstraint("genome_digest", "name", name="unique_genome_group_name"),)

    def __str__(self):
        return f"{self.genome.digest}/{self.name}"


class AssetLink(SQLModel, table=True):
    """
    A table to indicate the self-referential many-to-many relationship between Assets
    that encaptulates the parent-child relationship.
    """

    id: int | None = Field(default=None, primary_key=True)
    parent_digest: str | None = Field(default=None, foreign_key="asset.digest")
    child_digest: str | None = Field(default=None, foreign_key="asset.digest")


class AssetPublic(SQLModel):
    """
    A reference asset model.
    """

    digest: str | None = Field(
        max_length=64, min_length=64, primary_key=True, unique=True, nullable=True
    )
    name: str
    description: str | None
    recipe_id: int | None = Field(
        default=None,
        foreign_key="recipe.id",
        description="The recipe used to build the asset",
    )
    asset_group_id: int | None = Field(default=None, foreign_key="assetgroup.id")
    size: int | None = Field(sa_column=Column(BigInteger()))
    serving_modes_override: list[str] | None = Field(
        sa_column=Column(JSON, nullable=True),
        default=None,
    )
    colocate: list[dict[str, str]] | None = Field(
        sa_column=Column(JSON, nullable=True),
        default=None,
        description="Colocation metadata: parent files to symlink into this asset after pull",
    )

    @model_validator(mode="after")
    def validate_serving_modes_override(self) -> "AssetPublic":
        if self.serving_modes_override is not None:
            if not self.serving_modes_override:
                raise ValueError("serving_modes_override must not be empty when set")
            _check_serving_modes(self.serving_modes_override)
        return self


class Asset(AssetPublic, table=True):
    """
    A reference asset model.
    """

    asset_group_id: int | None = Field(default=None, foreign_key="assetgroup.id", index=True)
    asset_group: "AssetGroup" = Relationship(
        back_populates="assets", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    staged: list["StagedAsset"] = Relationship(
        back_populates="asset", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN
    )
    seek_keys: list["SeekKey"] = Relationship(
        back_populates="asset", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN
    )
    asset_names: list["AssetName"] = Relationship(
        back_populates="asset", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    path: str | None = Field(default=None)
    parents: list["Asset"] = Relationship(
        back_populates="children",
        link_model=AssetLink,
        sa_relationship_kwargs=dict(
            primaryjoin="Asset.digest==AssetLink.child_digest",
            secondaryjoin="Asset.digest==AssetLink.parent_digest",
            lazy="joined",
        ),
    )
    children: list["Asset"] = Relationship(
        back_populates="parents",
        link_model=AssetLink,
        sa_relationship_kwargs=dict(
            primaryjoin="Asset.digest==AssetLink.parent_digest",
            secondaryjoin="Asset.digest==AssetLink.child_digest",
            lazy="joined",
        ),
    )

    @property
    def serving_modes(self) -> list[str]:
        """
        Return the resolved serving modes for this asset.

        Uses the asset's own serving_modes override if set, otherwise falls back
        to the AssetClass serving_modes, otherwise defaults to ["archive"].
        """
        if self.serving_modes_override:
            return self.serving_modes_override
        if (
            self.asset_group
            and self.asset_group.asset_class
            and self.asset_group.asset_class.serving_modes
        ):
            return self.asset_group.asset_class.serving_modes
        return ["archive"]

    @property
    def asset_class_name(self) -> str | None:
        """Return the name of the asset's asset class, if resolvable."""
        if self.asset_group and self.asset_group.asset_class:
            return self.asset_group.asset_class.name
        return None

    @property
    def asset_group_name(self) -> str | None:
        """Return the name of the asset's group, if resolvable."""
        if self.asset_group:
            return self.asset_group.name
        return None

    @property
    def genome_digest(self) -> str | None:
        """Return the genome digest of the asset's group, if resolvable."""
        if self.asset_group:
            return self.asset_group.genome_digest
        return None

    @property
    def names(self) -> list["AssetName"]:
        """Every name row this content has in its group.

        Content is addressed by digest; names live in the ``assetname`` table,
        each carrying the provenance of the build that produced it under that
        name. ``AssetResponse`` serializes these whole -- not just the name
        strings -- so a client can adopt any name on pull *and* record the
        build behind it. Use ``an.name`` for the bare string.
        """
        return list(self.asset_names)

    @property
    def is_default(self) -> bool:
        """Whether this asset is its group's default (any of its names is).

        The flag lives on ``assetname`` rows; ``AssetResponse`` serializes this
        so the web UI's set-default control can show the current default
        without a per-group lookup.
        """
        return any(an.is_default for an in self.asset_names)

    @property
    def serves_files(self) -> bool:
        return "file" in self.serving_modes

    @property
    def serves_archive(self) -> bool:
        return "archive" in self.serving_modes

    @property
    def is_metadata_only(self) -> bool:
        return self.serving_modes == ["none"]

    @property
    def seek_keys_dict(self) -> dict[str, Path]:
        """
        Path-based seek keys only, as a dictionary of name -> relative path.
        Used in the build process for asset inputs; recipe templates access
        these as `values.assets[<name>].seek_keys_dict[<key>]`.

        Non-path seek keys (string, json) are excluded -- use seek() to access them.

        Returns:
            A dictionary of path-based seek keys with their corresponding relative paths.
        """
        if self.path is None:
            raise ValueError("Incomplete asset: path is not set")
        return {
            seek_key.name: Path(self.path) / seek_key.value
            for seek_key in self.seek_keys
            if is_path_type(seek_key.type)
        }

    @property
    def registry_path(self) -> str:
        """
        The registry path of the asset.

        Returns:
            str: The registry path.
        """
        return f"{self.asset_group.genome.digest}/{self.asset_group.name}:{self.name}"

    def __str__(self) -> str:
        return f"{self.registry_path}"


# ---------------------------------------------------------------------------
# Asset naming: three concepts, one content digest
#
# 1. Asset.name       -- A denormalized "canonical" name stamped on the Asset
#                        row at insert time.  Removal/update endpoints accept
#                        only this name.  It is a publication label, not a
#                        lookup axis.
#
# 2. AssetName.name   -- Rows in the assetname table.  An asset can have many
#                        names (aliases) within its group; this is the
#                        authoritative name axis.  E.g. "0.7.17" and "0.7.19"
#                        may be peer AssetName rows pointing at the same
#                        content digest.
#
# 3. AssetName.is_default -- A boolean flag on exactly one AssetName per
#                        group.  It controls what bare `group`-without-`:name`
#                        lookups resolve to.  This is explicitly NOT the same
#                        as the canonical name -- "default" can be a third
#                        peer alongside "0.7.17" and "0.7.19".
#
# So an asset can be *named* "0.7.17", *canonically named* "0.7.19" (the
# Asset.name value), with "default" as a separate AssetName flagged
# is_default -- and removal accepts only "0.7.19".
# ---------------------------------------------------------------------------


class AssetName(SQLModel, table=True):
    """
    A human-readable name for an asset's content, in a given asset group.

    Content is addressed by digest (``Asset.digest``); this table is the name
    axis (see the naming overview above). Unlike ``Alias``, asset names are
    unique per group, not globally -- every group has a ``default``.
    """

    id: int | None = Field(default=None, primary_key=True)
    name: str
    # Denormalized onto the row (derivable via asset_digest -> group) precisely
    # so the scoped uniqueness below is enforceable in the schema, not in prose.
    asset_group_id: int = Field(
        foreign_key="assetgroup.id", index=True, sa_column_kwargs={"nullable": False}
    )
    asset_digest: str = Field(foreign_key="asset.digest", index=True)
    # Group-default-ness lives here as a flag rather than as a pointer on
    # AssetGroup, so the same-group invariant holds by construction and a
    # deleted asset cannot leave a dangling default pointer.
    is_default: bool = Field(default=False)

    # --- Build provenance ---------------------------------------------------
    # Provenance describes *a build*, not the content it produced, so it lives
    # on the name row rather than as seek keys on the asset row. Two builds can
    # yield byte-identical content and share one Asset; each still gets its own
    # name row. Hanging provenance off the content silently dropped the second
    # build's, which is the bug these columns fix.
    #
    # All-NULL means the name was not created by a build: an incomplete-asset
    # placeholder, a pull from a server reporting none, or a legacy row whose
    # build could not be reconstructed from the catalog (see migration
    # 7c2f4b9e1a83 for which rows those are and why).
    build_digest: str | None = Field(default=None, index=True)
    #: Level-1 per-attribute digests, so build_digest is re-derivable from the row.
    build_level1: dict[str, str] | None = Field(
        sa_column=Column(JSON, nullable=True), default=None
    )
    #: Which digest algorithm build_digest/build_level1 were computed under.
    #: See refgenie.utils.build.BUILD_DIGEST_SCHEME.
    build_digest_scheme: str | None = Field(default=None)
    build_timestamp: datetime | None = Field(default=None)
    refgenie_version: str | None = Field(default=None)
    #: The level-2 record: every attribute of the build, inherent and
    #: non-inherent, in its raw form.
    inputs: dict[str, Any] | None = Field(sa_column=Column(JSON, nullable=True), default=None)
    docker_image: str | None = Field(default=None)
    docker_image_digest: str | None = Field(default=None)
    #: The recipe this build ran -- the asset is shared across builds, the recipe is not.
    recipe_id: int | None = Field(default=None, foreign_key="recipe.id")

    asset: "Asset" = Relationship(
        back_populates="asset_names", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )

    __table_args__ = (
        UniqueConstraint("asset_group_id", "name", name="unique_group_asset_name"),
        # At most one default name per group. Partial unique index -- SQLite
        # supports the WHERE clause; Postgres does too.
        Index(
            "unique_default_per_group",
            "asset_group_id",
            unique=True,
            sqlite_where=text("is_default"),
            postgresql_where=text("is_default"),
        ),
        # Not unique: two names can legitimately record the same build_digest
        # -- a rebuild with the same recipe/genome/inputs, or two rows whose
        # differences (a param's type, which slot an input filled) are not
        # inherent under the current scheme. Every insert still records full
        # provenance; a repeat is one build with two names, not a collision.
        Index(
            "unique_group_build_digest",
            "asset_group_id",
            "build_digest",
            unique=False,
        ),
    )

    def __str__(self) -> str:
        return self.name


class SeekKey(SQLModel, table=True):
    """
    A reference seek key model.
    """

    id: int | None = Field(default=None, primary_key=True)
    name: str
    value: str
    description: str | None
    type: SeekKeyType
    asset_digest: str | None = Field(default=None, foreign_key="asset.digest")
    asset: "Asset" = Relationship(
        back_populates="seek_keys", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    size: int | None = Field(sa_column=Column(BigInteger()))
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )


class AssetClassSeekKey(SQLModel, table=True):
    """
    A reference asset class seek key model.
    """

    id: int | None = Field(default=None, primary_key=True)
    name: str
    value: str | None = Field(default=None)
    description: str | None
    type: SeekKeyType
    asset_class_id: int = Field(default=None, foreign_key="assetclass.id")
    asset_class: Optional["AssetClass"] = Relationship(
        back_populates="seek_keys", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )

    def match_file(self, directory_path: Path) -> Path:
        """
        Match the file based on the seek key value and type.

        Args:
            directory_path: The directory path.
        """
        if not is_path_type(self.type):
            raise TypeError(
                f"Seek key '{self.name}' has type '{self.type.value}' which is not path-based. "
                f"Cannot match files for non-path seek keys."
            )
        error_string = f"files found that match '{self.value}' in {directory_path}"
        if not directory_path.exists():
            raise FileNotFoundError(f"Directory not found: {directory_path}")
        if self.type == SeekKeyType.file:
            # Ignore files starting with a dot
            matched_files = [
                f
                for f in directory_path.glob(re.sub(r"{[^}]+}", "*", self.value))
                if not f.name.startswith(".")
            ]
            if not matched_files:
                raise FileNotFoundError(f"No {error_string}")
            if len(matched_files) > 1:
                raise ValueError(f"Multiple {error_string}")
            return matched_files[0]

        elif self.type == SeekKeyType.prefix:
            # match common prefix of the files in the directory, ignore dotfiles
            files = [file for file in directory_path.iterdir() if not file.name.startswith(".")]
            if not files:
                raise FileNotFoundError(f"No {error_string}")
            common_prefix = os.path.commonprefix(list(set(file.name for file in files)))
            return directory_path / common_prefix.rstrip(".")

        elif self.type == SeekKeyType.directory:
            matched_directory_path = directory_path / self.value
            if not matched_directory_path.exists():
                raise FileNotFoundError(f"No {error_string}")
            return matched_directory_path


class AssetClassLink(SQLModel, table=True):
    """
    A table to indicate the self-referential many-to-many relationship between AssetClasses
    that encaptulates the parent-child relationship.
    """

    id: int | None = Field(default=None, primary_key=True)
    parent_id: int | None = Field(default=None, foreign_key="assetclass.id")
    child_id: int | None = Field(default=None, foreign_key="assetclass.id")


class AssetClassPublic(SQLModel):
    """
    A reference asset class model.
    """

    id: int | None = Field(default=None, primary_key=True)
    name: str
    version: str
    description: str | None
    serving_modes: list[str] = Field(
        sa_column=Column(JSON, nullable=False),
        default=["file"],
    )

    @model_validator(mode="after")
    def validate_serving_modes(self) -> "AssetClassPublic":
        if not self.serving_modes:
            raise ValueError("serving_modes must not be empty")
        _check_serving_modes(self.serving_modes)
        return self


class AssetClass(AssetClassPublic, table=True):
    """
    A reference asset class model.
    """

    asset_groups: list["AssetGroup"] = Relationship(
        back_populates="asset_class", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    recipes: list["Recipe"] = Relationship(
        back_populates="output_asset_class",
        sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE,
    )
    recipes_input: list["RecipeAssetClassesInputs"] = Relationship(
        back_populates="asset_class",
        sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN,
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    seek_keys: list["AssetClassSeekKey"] = Relationship(
        back_populates="asset_class",
        sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN,
    )
    children: list["AssetClass"] = Relationship(
        back_populates="parents",
        link_model=AssetClassLink,
        sa_relationship_kwargs=dict(
            primaryjoin="AssetClass.id==AssetClassLink.parent_id",
            secondaryjoin="AssetClass.id==AssetClassLink.child_id",
            lazy="joined",
        ),
    )
    parents: list["AssetClass"] = Relationship(
        back_populates="children",
        link_model=AssetClassLink,
        sa_relationship_kwargs=dict(
            primaryjoin="AssetClass.id==AssetClassLink.child_id",
            secondaryjoin="AssetClass.id==AssetClassLink.parent_id",
            lazy="joined",
        ),
    )

    __table_args__ = (UniqueConstraint("name", "version", name="unique_name_and_version"),)

    @property
    def serves_files(self) -> bool:
        """True if this asset class serves individual files."""
        return "file" in self.serving_modes

    @property
    def serves_archive(self) -> bool:
        """True if this asset class serves archive tarballs."""
        return "archive" in self.serving_modes

    @property
    def is_metadata_only(self) -> bool:
        """True if this asset class only serves metadata (no files or archives)."""
        return self.serving_modes == ["none"]

    def to_yaml(self) -> str:
        """
        Convert the model to a YAML string that could be used to create a new instance of the model.

        Returns:
            str: The YAML string.
        """
        asset_class_dict = self.model_dump(exclude={"created_at", "updated_at", "id"})
        asset_class_dict["seek_keys"] = {}
        for asset_class_seek_key in self.seek_keys:
            seek_key_entry = {
                "type": asset_class_seek_key.type.value,
                "description": asset_class_seek_key.description,
            }
            if asset_class_seek_key.value is not None:
                seek_key_entry["value"] = asset_class_seek_key.value
            asset_class_dict["seek_keys"][asset_class_seek_key.name] = seek_key_entry
        return ydump(asset_class_dict, default_flow_style=False)

    def __str__(self) -> str:
        return f"{self.name} (v{self.version})"


InputEntities = dict[str, dict[str, Any]]


class RecipePublic(SQLModel):
    """
    A reference recipe model.
    """

    id: int | None = Field(default=None, primary_key=True)
    name: str
    version: str
    description: str | None
    output_asset_class_id: int = Field(index=True, default=None, foreign_key="assetclass.id")
    command_templates: list[str] = Field(sa_column=Column(JSON, nullable=False))
    input_params: InputEntities | None = Field(sa_column=Column(JSON, nullable=True))
    input_files: InputEntities | None = Field(sa_column=Column(JSON, nullable=True))
    input_assets: InputEntities | None = Field(sa_column=Column(JSON, nullable=True))
    docker_image: str | None
    custom_seek_keys: dict[str, str] | None = Field(sa_column=Column(JSON), default_factory=dict)
    default_asset: str

    #: Ordered, gitignore-style glob list naming which of the files this recipe
    #: produces are inherent to the built asset's identity, and therefore which
    #: ones its content digest is computed over. Last match wins; a ``!`` prefix
    #: marks a match incidental. NULL means the recipe declares nothing, which
    #: includes every file. See ``refgenie.utils.build.get_dir_digest``.
    inherent: list[str] | None = Field(sa_column=Column(JSON, nullable=True), default=None)


class Recipe(RecipePublic, table=True):
    """
    A reference recipe model.
    """

    output_asset_class: "AssetClass" = Relationship(
        back_populates="recipes", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    input_asset_classes: Optional[list["RecipeAssetClassesInputs"]] = Relationship(
        back_populates="recipe", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )

    __table_args__ = (UniqueConstraint("name", "version"),)

    def to_yaml(self) -> str:
        """
        Convert the model to a YAML string that could be used to create a new instance of the model.

        Returns:
            str: The YAML string.
        """
        return ydump(
            self.model_dump(exclude={"created_at", "updated_at", "id"}),
            default_flow_style=False,
        )

    def __str__(self) -> str:
        return f"{self.name} (v{self.version})"


class RecipeAssetClassesInputs(SQLModel, table=True):
    """
    A table to indicate the many-to-many relationship between recipes and asset inputs
    """

    id: int | None = Field(default=None, primary_key=True)
    recipe_id: int = Field(index=True, default=None, foreign_key="recipe.id")
    asset_class_id: int = Field(index=True, default=None, foreign_key="assetclass.id")
    name: str
    default: str
    recipe: "Recipe" = Relationship(back_populates="input_asset_classes")
    asset_class: "AssetClass" = Relationship(back_populates="recipes_input")
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )


class DataChannelType(str, Enum):
    """
    The data channel type enumeration.
    """

    ftp = "ftp"
    http = "http"
    https = "https"
    local = "local"


class DataChannel(SQLModel, table=True):
    """
    A reference data channel model.
    """

    id: int | None = Field(default=None, primary_key=True)
    # double underscores are prohibited in the name
    name: str = Field(unique=True, regex=r"^[^_]+$")
    description: str | None
    type: DataChannelType
    index_address: str  # Full path/URL to the index.yaml file
    # Store encrypted credentials as a string
    encrypted_credentials: str | None = Field(sa_column=Column(JSON), default=None)
    configuration_id: int = Field(default=None, foreign_key="configuration.id")
    configuration: Optional["Configuration"] = Relationship(
        back_populates="data_channels",
        sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE,
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )

    def __str__(self) -> str:
        return f"{self.name} ({self.type.value})"


class RemoteAssetLink(SQLModel, table=True):
    """Links a Remote to the StagedAssets it serves.

    Records are created during staging with pushed=False (intent to push).
    refgenie push sets pushed=True after successful upload.
    """

    remote_id: int = Field(foreign_key="remote.id", primary_key=True)
    asset_digest: str = Field(primary_key=True)
    mode: str = Field(primary_key=True)
    pushed: bool = Field(default=False)

    __table_args__ = (
        ForeignKeyConstraint(
            ["asset_digest", "mode"],
            ["stagedasset.asset_digest", "stagedasset.mode"],
        ),
    )


class Remote(SQLModel, table=True):
    """
    A reference remotes model.
    """

    id: int | None = Field(default=None, primary_key=True)
    prefix: str
    type: RemoteType = Field()
    description: str | None
    push_command: str | None = Field(default=None)
    configuration_id: int = Field(default=None, foreign_key="configuration.id")
    configuration: Optional["Configuration"] = Relationship(
        back_populates="remotes", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    staged_assets: list["StagedAsset"] = Relationship(
        back_populates="remotes",
        link_model=RemoteAssetLink,
        sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE,
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )

    def __str__(self) -> str:
        return f"Remote {self.type.value} ({self.prefix})"


class StorePublic(SQLModel):
    """A refget store the server federates over.

    The ``store`` table is the single source of truth for what a refgenie
    service serves. Each row is one physically separate store; content is
    digest-addressed, so identical genomes across stores dedup automatically.
    The only true conflict is two stores mapping the same alias to different
    collection digests, resolved by ``priority`` (lower integer wins).
    """

    name: str = Field(primary_key=True)
    #: Store root URL (remote) or local ``.refget_store`` path (on_disk).
    url: str
    type: StoreType = Field(default=StoreType.remote)
    #: Lower integer = higher priority. The highest-priority store wins any
    #: alias tie and is the default mount for service-info.
    priority: int = Field(default=100)
    enabled: bool = Field(default=True)
    description: str | None = None


class Store(StorePublic, table=True):
    """A refget store the server federates over."""

    genomes: list["Genome"] = Relationship(
        back_populates="store", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )

    def __str__(self) -> str:
        return f"Store {self.name} ({self.type.value}, priority={self.priority})"


class ConfigurationPublic(SQLModel):
    """
    A reference configuration model.
    """

    version: int = Field(unique=True)
    servers: list[str] = Field(sa_column=Column(JSON), default_factory=list)
    genome_folder: str
    genome_stage_folder: str | None


class Configuration(ConfigurationPublic, table=True):
    """
    A reference configuration model.
    """

    id: int | None = Field(default=None, primary_key=True)
    remotes: list["Remote"] = Relationship(
        back_populates="configuration",
        sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN,
    )
    data_channels: list["DataChannel"] = Relationship(
        back_populates="configuration",
        sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_ORPHAN,
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )


class StagedAssetPublic(SQLModel):
    asset_digest: str = Field(default=None, foreign_key="asset.digest", index=True)
    mode: str = Field(max_length=10)  # "file" or "archive"
    directory_contents: list[str] = Field(sa_column=Column(JSON))
    build_commands: list[str] = Field(default=None, sa_column=Column(JSON))
    download_count: int = Field(default=0)
    tarball_digest: str | None = Field(default=None, max_length=64)
    tarball_size: int | None = Field(default=None, sa_column=Column(BigInteger()))


class StagedAsset(StagedAssetPublic, table=True):
    """
    A staged asset record. Each asset can have up to two StagedAsset records:
    one for mode="file" and one for mode="archive".

    File-mode: a directory symlink in genome_stage_folder pointing to genome_folder.
    Archive-mode: a real .tgz tarball in genome_stage_folder.
    """

    asset: "Asset" = Relationship(
        back_populates="staged", sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE
    )
    remotes: list["Remote"] = Relationship(
        back_populates="staged_assets",
        link_model=RemoteAssetLink,
        sa_relationship_kwargs=SA_RELATIONSHIP_KWARGS_BASE,
    )
    updated_at: datetime | None = Field(
        sa_column=Column(DateTime(), onupdate=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )
    created_at: datetime | None = Field(
        sa_column=Column(DateTime(), default=func.now()),
        default_factory=lambda: datetime.now(timezone.utc),
    )

    # Logical identity is (asset_digest, mode), enforced as a unique
    # constraint; SQLModel needs a single PK field, so id is the physical PK.
    __table_args__ = (UniqueConstraint("asset_digest", "mode", name="unique_asset_mode"),)
    id: int | None = Field(default=None, primary_key=True)

    def __str__(self):
        return f"StagedAsset({self.asset_digest}, mode={self.mode})"


class AlembicVersion(SQLModel, table=True):
    """
    Table to manage the alembic version:

    """

    __tablename__ = "alembic_version"

    version_num: str = Field(primary_key=True)
