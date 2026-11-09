"""Package exception types"""

__all__ = [
    "AssetClassExistsError",
    "AssetExistsError",
    "BuildFailedError",
    "CatalogConflictError",
    "ConfigError",
    "ContainerEngineMissingError",
    "CustomSeekKeyError",
    "DatabaseMigrationError",
    "FederatedAliasError",
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
    "NoSubscriptionsError",
    "ParentPullError",
    "PullFailedError",
    "PullSkipped",
    "RecipeExistsError",
    "RefgenieError",
    "RemoteDigestMismatchError",
    "FastaDigestMismatchError",
    "ServerCannotServe",
    "SpeciesFileError",
    "StoreExistsError",
    "UnknownStoreError",
    "UnservedStoreError",
]


class RefgenieError(Exception):
    """Base exception type for this package"""


class MissingStagedAssetError(RefgenieError):
    """Error type for request of an unavailable staged asset."""

    def __init__(self, context: str):
        super(MissingStagedAssetError, self).__init__(f"Staged asset not found: {context}")


class ConfigError(RefgenieError):
    """Error in refgenie configuration."""

    pass


class ContainerEngineMissingError(RefgenieError):
    """A recipe runs in a container image, and this host has nothing to run it with."""

    def __init__(self, image: str):
        super().__init__(
            f"This recipe runs in a container image ({image}); install Docker or Apptainer."
        )
        self.image = image


class CustomSeekKeyError(RefgenieError):
    """A recipe's custom seek key command did not produce a usable value.

    The web build form shows the message as is, so it is written for the person
    building, with the failed command last for the recipe's maintainer.
    """

    def __init__(
        self,
        seek_key: str,
        command: str,
        reason: str,
        image: str | None = None,
        engine: str | None = None,
    ):
        if image is None:
            where = "on this computer"
        elif engine is None:
            where = f"inside the container image {image}"
        else:
            where = f"inside the container image {image}, with {str(engine).capitalize()}"
        super().__init__(
            f"Refgenie could not work out a name for this asset. The recipe names it after "
            f"the {seek_key} of the tool it runs, which it looks up by running a command "
            f"{where}. {reason} "
            f"You can type your own name instead (the Asset name box, or --asset on the "
            f"command line) and the build will go ahead. "
            f"Whoever maintains this recipe needs the command that failed: {command}"
        )
        self.seek_key = seek_key
        self.command = command
        self.image = image
        self.engine = engine


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


class UnservedStoreError(MissingGenomeError):
    """A genome lives in a store this refgenie does not serve.

    A :class:`MissingGenomeError`, so callers that answer 404 for a missing
    genome answer 404 here too.
    """

    def __init__(self, genome, store_name: str | None):
        self.genome = genome
        self.store_name = store_name
        if store_name is None:
            message = (
                f"Genome '{genome}' has no store, and this server has no store of its own "
                "to read it from."
            )
        else:
            message = (
                f"Genome '{genome}' is in store '{store_name}', which this server does not serve."
            )
        RefgenieError.__init__(self, message)


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

    def __init__(self, ref):
        super(MissingRemoteError, self).__init__(f"Remote '{ref}' not found.")


class MissingStoreError(RefgenieError):
    """Error type for request of an unavailable store."""

    def __init__(self, name):
        super(MissingStoreError, self).__init__(f"Store '{name}' not found.")


class StoreExistsError(RefgenieError):
    """Error type raised when adding a store whose name already exists."""

    def __init__(self, name):
        super(StoreExistsError, self).__init__(f"Store '{name}' already exists.")


class UnknownStoreError(RefgenieError, ValueError):
    """A genome filter named a store this server did not load."""

    def __init__(self, name: str, known: list[str]):
        listed = ", ".join(known) if known else "none"
        super().__init__(f"Unknown store '{name}'. Known stores: {listed}.")


class SpeciesFileError(RefgenieError):
    """A store's species.json could not be read, or its top level is not the documented shape."""


class CatalogConflictError(RefgenieError):
    """A write clashed with a row another writer committed first."""


class DatabaseMigrationError(RefgenieError):
    """The database schema could not be brought to the version this code needs."""

    def __init__(self, current, target, cause):
        super(DatabaseMigrationError, self).__init__(
            f"Could not migrate the database from schema '{current}' to '{target}': {cause}"
        )


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


class FastaDigestMismatchError(RefgenieError, ValueError):
    """A FASTA file's sequence collection digest is not the one expected."""

    def __init__(self, fasta_name, expected_digest, actual_digest):
        self.expected_digest = expected_digest
        self.actual_digest = actual_digest
        super().__init__(f"{fasta_name} has digest {actual_digest}, expected {expected_digest}")


class PullFailedError(RefgenieError):
    """Pull operation failed.

    ``reasons`` maps each server that was tried to why it could not serve the
    asset. It is filled in when the puller has run out of servers, and empty
    for a failure that is not about any one server.
    """

    def __init__(self, message: str = "", reasons: dict[str, str] | None = None):
        super().__init__(message)
        self.reasons: dict[str, str] = dict(reasons or {})


class ParentPullError(PullFailedError):
    """An asset needs a parent that is not here, and the parent could not be pulled.

    The message names the parent (genome, asset group, asset name, digest),
    the asset that needs it, and why the parent's pull failed.
    """


class NoSubscriptionsError(PullFailedError):
    """There is no server to pull from: nothing is subscribed."""

    def __init__(
        self,
        message: str = (
            "No server subscriptions found. Subscribe to a server first "
            "(refgenie subscribe -s <url>)."
        ),
    ):
        super().__init__(message)


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
    """The pull was skipped (user declined, or force=False).

    Raised inside a pull attempt so the PullTransaction rolls back what that
    attempt created; pull() catches it and returns None without trying any
    further servers.
    """

    pass


class MissingBuildInputError(RefgenieError):
    """Missing required build input (files, params, assets)."""

    pass


class BuildFailedError(RefgenieError):
    """A recipe ran but did not produce its asset."""

    pass
