import re
from functools import partial
from pathlib import Path
from typing import Annotated, Any

from pydantic import BaseModel
from pydantic.functional_validators import BeforeValidator

from refgenie.const import (
    ASSET_NAME_CHAR_WHITELIST,
    ASSET_REGISTRY_PATH_COMPONENT_CHAR_BLACKLIST,
)
from refgenie.db.tables import Asset, Recipe


def validate_string(
    v: str,
    character_whitelist: set[str] | None = None,
    character_blacklist: set[str] | None = None,
) -> str | None:
    """
    Validates the input string by checking if it contains any of the blacklisted
    characters and if it contains any characters not in the whitelisted set.

    The behavior is dictated by whether the character_whitelist and character_blacklist
    are provided, otherwise the corresponding checks are skipped.

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
