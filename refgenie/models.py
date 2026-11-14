import re
from collections.abc import Iterable
from dataclasses import dataclass, field
from functools import partial
from pathlib import Path
from typing import Annotated, Any, Generic, TypeVar

from pydantic import BaseModel
from pydantic.functional_validators import BeforeValidator
from pydantic_core import core_schema

from refgenie.const import (
    ASSET_NAME_CHAR_WHITELIST,
    ASSET_REGISTRY_PATH_COMPONENT_CHAR_BLACKLIST,
)
from refgenie.db.tables import Asset, Genome, Recipe, Species

T = TypeVar("T")


def validate_string(
    v: str,
    character_whitelist: set[str] | None = None,
    character_blacklist: set[str] | None = None,
) -> str | None:
    """
    Check a string against an optional character whitelist and blacklist.

    Args:
        v: The input string to validate.
        character_whitelist: The set of whitelisted characters.
        character_blacklist: The set of blacklisted characters.

    Returns:
        The validated string.

    Raises:
        AssertionError: If the input string contains any of the blacklisted characters or
            if it contains any characters not in the whitelisted set.
    """
    if character_whitelist is not None:
        assert all(char in character_whitelist for char in v), (
            f"Invalid string: {v}. Strings can only contain the following characters: "
            f"{character_whitelist}"
        )
    if character_blacklist is not None:
        assert all(char not in character_blacklist for char in v), (
            f"Invalid string: {v}. Strings cannot contain the following characters: "
            f"{character_blacklist}"
        )
    return v


validate_asset_registry_path_component = partial(
    validate_string,
    character_blacklist=ASSET_REGISTRY_PATH_COMPONENT_CHAR_BLACKLIST,
    character_whitelist=None,
)
validate_asset_name = partial(
    validate_string,
    character_whitelist=ASSET_NAME_CHAR_WHITELIST,
    character_blacklist=None,
)


AssetRegistryPathComponent = Annotated[str, BeforeValidator(validate_asset_registry_path_component)]
AssetNameStr = Annotated[
    str,
    BeforeValidator(validate_asset_registry_path_component),
    BeforeValidator(validate_asset_name),
]


_GENOME_DIGEST_RE = re.compile(r"[A-Za-z0-9_-]{32}")


class GenomeDigest(str):
    """The sequence-collection (level 0) digest that identifies a genome.

    32 characters of URL-safe base64 (GA4GH sha512t24u); a malformed value
    raises ValueError. It validates a known digest and is never used to tell a
    digest from an alias.
    """

    def __new__(cls, value: str) -> "GenomeDigest":
        if not isinstance(value, str) or not _GENOME_DIGEST_RE.fullmatch(value):
            raise ValueError(f"'{value}' is not a genome digest (32 URL-safe base64 characters)")
        return super().__new__(cls, value)

    @classmethod
    def __get_pydantic_core_schema__(cls, source: Any, handler: Any) -> core_schema.CoreSchema:
        return core_schema.no_info_after_validator_function(cls, core_schema.str_schema())


class GenomeAlias(str):
    """A human-readable genome name, e.g. ``hg38``.

    Construct with ``GenomeAlias(s)``. An empty name, or one containing ``/``
    (which would split a registry path), raises ValueError. A ``store::name``
    qualified alias is allowed.
    """

    @staticmethod
    def is_valid(value: object) -> bool:
        """Whether ``value`` can be a genome alias: a non-empty str with no '/'."""
        return isinstance(value, str) and bool(value) and "/" not in value

    def __new__(cls, value: str) -> "GenomeAlias":
        if not cls.is_valid(value):
            raise ValueError(f"'{value}' is not a genome alias (non-empty, no '/')")
        return super().__new__(cls, value)

    @classmethod
    def __get_pydantic_core_schema__(cls, source: Any, handler: Any) -> core_schema.CoreSchema:
        return core_schema.no_info_after_validator_function(cls, core_schema.str_schema())


def display_alias(names: Iterable[str]) -> GenomeAlias | None:
    """The shortest plain alias of ``names`` (ties alphabetical), the name to
    show a person: no iGenomes path ('/') and no store-qualified 'store::name'."""
    plain = [n for n in names if GenomeAlias.is_valid(n) and "::" not in n]
    return GenomeAlias(min(plain, key=lambda n: (len(n), n))) if plain else None


def first_usable_alias(names: Iterable[str]) -> GenomeAlias | None:
    """The first of ``names`` that is a valid alias here, skipping iGenomes paths
    such as 'Homo_sapiens/NCBI/GRCh38'; None if there is none."""
    return next((GenomeAlias(n) for n in names if GenomeAlias.is_valid(n)), None)


class BuildParams(BaseModel):
    """
    A model for build input parameters which include:

    - optional dict mapping asset input to value, like `{'fasta': 'hg38/fasta:default'}"`
    - optional dict mapping parameter name to value, like `{'cores': 4}`
    - optional dict mapping file name to file path, like `{'fasta': '/path/to/fasta.fa'}`
    """

    assets: dict[str, str] | None = None
    params: dict[str, Any] | None = None
    files: dict[str, Path] | None = None

    def populate_with_defaults_from_recipe(self, recipe: Recipe) -> None:
        """
        Populates the build parameters with default values from the recipe.

        ```
        input_params:
            mersize:
                default: "30"
                description: The mer size.
            minocc:
                default: "2"
                description: The minimum occurrence number for the mers to index.
        ```

        Args:
            recipe: The recipe to get the default values from.
        """
        # implemented only for params for now
        if self.params or recipe.input_params is None:
            return
        self.params = {k: v["default"] for k, v in recipe.input_params.items() if "default" in v}


class AssetRegistryPathComponents(BaseModel):
    """
    A model for asset registry path components.

    The asset registry path is a string that represents a reference asset in the
    refgenie registry. It is composed of the following components:

    - genome: the genome name, optional
    - asset_group: the asset group name
    - seek_key: the seek key name, optional
    - asset: the asset name, optional
    """

    protocol: AssetRegistryPathComponent | None = None
    genome: AssetRegistryPathComponent | None = None
    asset_group: AssetRegistryPathComponent
    seek_key: AssetRegistryPathComponent | None = None
    asset: AssetNameStr | None = None

    @classmethod
    def parse_registry_path(cls, asset_registry_path: str) -> "AssetRegistryPathComponents":
        """
        Parses asset registry path it into its components.
        """
        components = {}

        if "/" not in asset_registry_path:
            asset_registry_path = f"/{asset_registry_path}"

        pattern = re.compile(
            r"(?:(?P<protocol>[0-9a-zA-Z._-]+)(?:::|:\/\/))?(?P<genome>\S+)?\/(?P<asset_group>[^:]+)(?:\:(?P<asset>\S+))?$"
        )

        match = pattern.match(asset_registry_path)

        if match:
            if match.group("protocol"):
                components["protocol"] = match.group("protocol")
            components["genome"] = match.group("genome")
            components["asset"] = match.group("asset")
            asset_group = match.group("asset_group")
            if "." in asset_group:
                components["asset_group"] = asset_group.split(".")[0]
                components["seek_key"] = asset_group.split(".")[1]
            else:
                components["asset_group"] = asset_group
        return cls(**components)

    def to_registry_path(self) -> str:
        """
        Converts the asset registry path components into a string.
        """
        components = []
        if self.protocol:
            components.append(f"{self.protocol}::")
        if self.genome:
            components.append(f"{self.genome}/")
        components.append(f"{self.asset_group}")
        if self.seek_key:
            components.append(f".{self.seek_key}")
        if self.asset:
            components.append(f":{self.asset}")
        return "".join(components)


class BuildCommandValues(BaseModel):
    """
    A model for build command values which include:

    - asset_group_name: a `str` representing the asset group name
    - genome_digest: a `str` representing the digest uniquely identifying the genome
    - output_folder: a `pathlib.Path` object indicating the output folder
    - assets: a `dict[str, Asset]` mapping asset group names to resolved asset names: specified by the user, or defaults
    - params: a `dict[str, Any]` mapping of parameter names to values: specified by the user, or defaults
    - files: a `dict[str, Path]` mapping of file names to file paths: specified by the user, or defaults
    - custom_seek_keys: a `dict[str, Any]` mapping of custom seek key values, may include software versions, etc.
    """

    asset_group_name: str
    genome_digest: str
    assets: dict[str, Asset] | None = None
    params: dict[str, Any] | None = None
    files: dict[str, Path] | None = None
    output_folder: Path | None = None
    genome_folder: Path
    custom_seek_keys: dict[str, Any] | None = None
    refget_store_path: str | None = None


class InputEntity(BaseModel):
    """
    A model for input entities, like files, parameters, and assets.
    """

    description: str
    default: str | int | None = None

    def __repr__(self) -> str:
        return f"({self.description}) default={self.default}"

    def __str__(self) -> str:
        return self.__repr__()


class DataChannelSyncReport(BaseModel):
    """What one ``SourceManager.sync_channel`` run did, item by item.

    Shared by the CLI and the local web action.
    """

    channel: str
    asset_classes_added: int = 0
    asset_classes_skipped: int = 0
    asset_classes_failed: int = 0
    recipes_added: int = 0
    recipes_skipped: int = 0
    recipes_failed: int = 0
    #: One human-readable line per failed item, in the order they were tried.
    errors: list[str] = []

    @property
    def ok(self) -> bool:
        return not self.asset_classes_failed and not self.recipes_failed


@dataclass(frozen=True)
class Page(Generic[T]):
    """One page of a list read, as a manager returns it.

    Services turn it into their own wire envelope. A dataclass, not a Pydantic
    model, so table rows pass through unvalidated.
    """

    items: list[T]
    offset: int
    limit: int
    total: int


@dataclass(frozen=True)
class GenomeListing:
    """A genome as the listings show it: the row, its names, its asset count, its species."""

    genome: Genome
    #: From the mode-selected alias backend, read for the whole page at once.
    aliases: list[str]
    asset_count: int
    #: The genome's species from the species list; None when the list has no
    #: row for its taxon (or it has no taxon).
    species: Species | None = None

    @property
    def common_name(self) -> str | None:
        """The species list's common name for the genome's taxon."""
        return self.species.common_name if self.species else None


@dataclass(frozen=True)
class GenomeDetail(GenomeListing):
    """A genome listing plus its FHR record from the store that owns it."""

    fhr: dict | None = None
    #: Base URL (ending in ``/``) of the remote refget store holding the
    #: genome's sequences; None for on-disk or unknown stores.
    store_url: str | None = None


@dataclass(frozen=True)
class AliasEntry:
    """One ``(namespace, alias)`` pair from a refget store's alias tables."""

    namespace: str
    alias: str


@dataclass(frozen=True)
class SequenceRow:
    """One sequence of a genome's collection, as the genome page lists it."""

    name: str
    length: int
    #: Bare GA4GH digest, without the ``SQ.`` prefix level 2 carries.
    sha512t24u: str
    md5: str | None = None
    alphabet: str | None = None
    description: str | None = None
    aliases: list[AliasEntry] = field(default_factory=list)


@dataclass(frozen=True)
class CollectionSequences:
    """A genome's sequence collection, row by row, from the store that owns it."""

    digest: str
    store: str
    #: Level-1 attribute digests (``names``, ``lengths``, ``sequences``, ...).
    level1: dict[str, str]
    sequences: list[SequenceRow]
    collection_aliases: list[AliasEntry]
