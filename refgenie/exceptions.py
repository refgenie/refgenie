"""Package exception types"""

import abc

__all__ = [
    "AssetClassExistsError",
    "AssetExistsError",
    "ConfigError",
    "DownloadJsonError",
    "FederatedAliasError",
    "InvalidSeqColError",
    "MissingAliasError",
    "MissingAssetClassError",
    "MissingAssetError",
    "MissingAssetGroupError",
    "MissingBuildInputError",
    "MissingConfigDataError",
    "MissingGenomeError",
    "MissingRecipeError",
    "MissingRemoteError",
    "MissingSeekKeyError",
    "MissingStagedAssetError",
    "MissingStoreError",
    "NoArchiveError",
    "PullFailedError",
    "PullSkipped",
    "RecipeExistsError",
    "RefgenieError",
    "RemoteDigestMismatchError",
    "ServerCannotServe",
    "StoreExistsError",
    "UnboundEnvironmentVariablesError",
]

class RefgenieError(Exception):
    """Base exception type for this package"""

    __metaclass__ = abc.ABCMeta


class InvalidSeqColError(RefgenieError):
    """Object was not validated successfully according to schema."""

    def __init__(self, message, errors):
        super().__init__(message)
        self.message = message
        self.errors = errors

    def __str__(self):
        return f"InvalidSeqColError ({self.message}): {self.errors}"


class MissingStagedAssetError(RefgenieError):
    """Error type for request of an unavailable staged asset."""

    def __init__(self, context: str):
        super(MissingStagedAssetError, self).__init__(f"Staged asset not found: {context}")


class DownloadJsonError(RefgenieError):
    """Non-OK response from a JSON download attempt"""

    def __init__(self, resp):
        super(DownloadJsonError, self).__init__(
            "No response provided" if resp is None else f"Reason: {resp}"
        )
        self.response = resp


class ConfigError(RefgenieError):
    """Error in refgenie configuration."""

    pass


class MissingAliasError(RefgenieError):
    """Error type for request of an unavailable alias."""

    def __init__(self, alias):
        super(MissingAliasError, self).__init__(f"Genome alias '{alias}' not found.")


class FederatedAliasError(RefgenieError):
    """A local alias would shadow a name a federated store already claims."""

    def __init__(self, alias, genome_digest):
        super(FederatedAliasError, self).__init__(
            f"Alias '{alias}' already names federated genome {genome_digest}. "
            f"Choose another name, or remove the store that supplies it."
        )


class MissingGenomeError(RefgenieError):
    """Error type for request of unknown genome/assembly."""

    def __init__(self, genome):
        super(MissingGenomeError, self).__init__(f"Genome '{genome}' not found.")


class MissingAssetGroupError(RefgenieError):
    """Error type for request of an unavailable asset group."""

    def __init__(self, genome, asset_group):
        super(MissingAssetGroupError, self).__init__(
            f"Asset group '{genome}/{asset_group}' not found."
        )


class MissingAssetError(RefgenieError):
    """Error type for request of an unavailable asset."""

    def __init__(
        self,
        genome: str | None = None,
        asset_group: str | None = None,
        asset: str | None = None,
        digest: str | None = None,
    ):
        # either digest or genome/asset_group/asset must be provided
        if digest is None and (genome is None or asset_group is None or asset is None):
            raise ValueError(
                "Either 'digest' or 'genome', 'asset_group', and 'asset' must be provided."
            )
        super(MissingAssetError, self).__init__(
            f"Asset '{genome}/{asset_group}:{asset}' not found."
            if digest is None
            else f"Asset identified with digest '{digest}' not found."
        )


class MissingSeekKeyError(RefgenieError):
    """Error type for request of an unavailable asset seek key."""

    def __init__(self, genome, asset_group, asset, seek_key):
        super(MissingSeekKeyError, self).__init__(
            f"Seek key '{genome}/{asset_group}.{seek_key}:{asset}' not found."
        )


class MissingRecipeError(RefgenieError):
    """Error type for request of an unavailable recipe."""

    def __init__(self, recipe):
        super(MissingRecipeError, self).__init__(f"Recipe '{recipe}' not found.")


class RecipeExistsError(RefgenieError):
    """Recipe with this name and version already exists."""

    def __init__(self, recipe_name, recipe_version):
        super(RecipeExistsError, self).__init__(
            f"Recipe '{recipe_name}' version '{recipe_version}' already exists."
        )


class MissingRemoteError(RefgenieError):
    """Error type for request of an unavailable remote."""

    def __init__(self, remote_type):
        super(MissingRemoteError, self).__init__(f"Remote '{remote_type.value}' not found.")


class MissingStoreError(RefgenieError):
    """Error type for request of an unavailable store."""

    def __init__(self, name):
        super(MissingStoreError, self).__init__(f"Store '{name}' not found.")


class StoreExistsError(RefgenieError):
    """Error type raised when adding a store whose name already exists."""

    def __init__(self, name):
        super(StoreExistsError, self).__init__(f"Store '{name}' already exists.")


class MissingAssetClassError(RefgenieError):
    """Error type for request of an unavailable asset class."""

    def __init__(self, asset_class):
        super(MissingAssetClassError, self).__init__(f"Asset class '{asset_class}' not found.")


class AssetClassExistsError(RefgenieError):
    """Asset class with this name and version already exists."""

    def __init__(self, asset_class_name, asset_class_version):
        super(AssetClassExistsError, self).__init__(
            f"Asset class '{asset_class_name}' version '{asset_class_version}' already exists."
        )


class MissingConfigDataError(RefgenieError):
    """Missing required configuration instance items"""

    pass


class UnboundEnvironmentVariablesError(RefgenieError):
    """Use of environment variable that isn't bound to a value."""

    pass


class RemoteDigestMismatchError(RefgenieError):
    """Remote digest of the parent asset does not match its local counterpart"""

    def __init__(self, asset, local_digest, remote_digest):
        msg = (
            "This asset is built from parent asset '{}', but for this parent, the remote does not "
            "match your local asset (local: {}; remote: {}). Refgenie will not pull this asset "
            "because the remote version was not built from the same parent asset you have locally.".format(
                asset, local_digest, remote_digest
            )
        )
        super(RemoteDigestMismatchError, self).__init__(msg)


class PullFailedError(RefgenieError):
    """Pull operation failed."""

    pass


class AssetExistsError(PullFailedError):
    """Asset already exists (use force to overwrite)."""

    pass


class NoArchiveError(PullFailedError):
    """No archive found on server."""

    pass


class ServerCannotServe(RefgenieError):
    """A particular server cannot serve the requested asset.

    Raised inside a pull attempt so the PullTransaction rolls back what that
    attempt created; the server loop catches it, logs, and tries the next
    server. ``final_error`` carries the exception pull() should raise if no
    server can serve the asset.
    """

    def __init__(self, message: str, final_error: Exception | None = None):
        super().__init__(message)
        self.final_error = final_error


class PullSkipped(RefgenieError):
    """The pull was deliberately skipped (user declined, or force=False).

    Raised inside a pull attempt so the PullTransaction rolls back what that
    attempt created; pull() catches it and returns None without trying any
    further servers.
    """

    pass


class MissingBuildInputError(RefgenieError):
    """Missing required build input (files, params, assets)."""

    pass
