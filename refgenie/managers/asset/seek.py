"""
AssetSeekMixin - seek, seek-key, and group-default resolution for AssetManager.

Owns the read side of asset addressing: resolving a registry path to a local
alias-tree path or content path, remote file-level URLs, seek-key lookups, and
the per-group default asset name.

This is a mixin, not a collaborator object, so every method keeps its name on
``AssetManager``; it is not usable standalone.
"""

from pathlib import Path
from typing import Any, TYPE_CHECKING

from sqlalchemy import update as sa_update
from sqlalchemy.orm import selectinload
from sqlmodel import Session, select

from refgenie.db.tables import (
    Asset,
    AssetClassSeekKey,
    AssetName,
    SeekKey,
    SeekKeyType,
    is_path_type,
)
from refgenie.exceptions import (
    MissingAliasError,
    MissingAssetError,
    MissingSeekKeyError,
)
from refgenie.logger import logger
from refgenie.managers.asset.queries import seek_key_by_name_stmt
from refgenie.managers.queries import one_or_raise
from refgenie.managers.sources.api_ids import API_ID_ALIAS_DIGEST
from refgenie.models import AssetRegistryPathComponents

if TYPE_CHECKING:
    from refgenie.managers.alias import AliasManager
    from refgenie.managers.genome import GenomeManager
    from refgenie.managers.sources.manager import SourceManager


def _payload_is_path_type(seek_key_type: "str | SeekKeyType | None") -> bool:
    """Return True if a server payload's seek-key ``type`` is a path type.

    The value arrives from JSON as the enum's string value (e.g. ``"file"``),
    but may also be a :class:`SeekKeyType`. Both are normalized here.
    """
    if seek_key_type is None:
        return False
    try:
        return is_path_type(SeekKeyType(seek_key_type))
    except ValueError:
        return False


class AssetSeekMixin:
    """
    Seek/seek-key/default resolver mixed into ``AssetManager``.

    Relies on state and methods provided by ``AssetManager`` (attributes below,
    plus catalog lookups like ``get``/``get_group`` and
    ``_resolve_genome_digest``), resolved through the MRO at runtime.
    """

    # Provided by AssetManager.__init__ / ResourceManager.
    _genome_folder: Path
    _alias_folder: Path
    _alias_manager: "AliasManager"
    _genome_manager: "GenomeManager"
    _source_manager: "SourceManager"

    # === Seek Operations ===

    def seek(
        self,
        genome_name: str,
        asset_group_name: str,
        asset_name: str | None = None,
        seek_key_name: str | None = None,
        force_exists: bool = False,
        abs_path: bool = False,
    ) -> str:
        """
        Get the path to an asset file, or the value of a non-path seek key.

        By default the returned path is under the human-readable alias tree
        (``alias/<alias>/<group>/<name>/...``), named for the alias the caller
        supplied. Pass ``abs_path=True`` for the digest-addressed content path
        under ``data/`` instead.

        Args:
            genome_name: The name of the genome (an alias).
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            seek_key_name: The name of the seek key.
            force_exists: Whether to raise an error if the path doesn't exist.
            abs_path: Return the content (``data/``) path rather than the alias path.

        Returns:
            str: The path to the asset file (for path-based seek keys), or the
                value directly (for non-path seek keys like string or json).
        """
        return self._seek_by_components(
            asset_registry_path_components=AssetRegistryPathComponents(
                genome=genome_name,
                asset_group=asset_group_name,
                asset=asset_name,
                seek_key=seek_key_name,
            ),
            force_exists=force_exists,
            abs_path=abs_path,
        )

    def _seek_by_components(
        self,
        asset_registry_path_components: AssetRegistryPathComponents,
        force_exists: bool = False,
        abs_path: bool = False,
    ) -> str:
        """
        Get the path to an asset file from registry path components, or the value
        directly for non-path seek keys.

        The default result is the alias-tree path built from exactly the alias
        and asset name the caller supplied, with the seek value's genome-digest
        substring rewritten to the alias (matching the alias tree's per-file
        rename). ``abs_path`` returns the digest-addressed content path instead.

        Args:
            asset_registry_path_components: The asset registry path components.
            force_exists: Whether to raise an error if the path doesn't exist.
            abs_path: Return the content (``data/``) path rather than the alias path.

        Returns:
            str: The path to the asset file (for path-based seek keys), or the
                value directly (for non-path seek keys like string or json).
        """
        logger.debug(f"Seeking {asset_registry_path_components}")
        if (
            asset_registry_path_components.genome is None
            or asset_registry_path_components.asset_group is None
        ):
            raise ValueError("Genome name and asset name are required")
        genome_token = asset_registry_path_components.genome
        # The token may be an alias or a genome digest. Aliases win; a bare
        # digest for a known genome is accepted too. ``supplied_alias`` records
        # which form the caller used, so the alias-tree path below can key off
        # the alias the caller named (and only when they named one).
        try:
            genome_digest = self._alias_manager.resolve(genome_token)
            supplied_alias = genome_token
        except MissingAliasError:
            if not self._genome_manager.exists(genome_token):
                raise
            genome_digest = genome_token
            supplied_alias = None
        asset_group_name = asset_registry_path_components.asset_group
        asset_name = asset_registry_path_components.asset or self.get_default(
            asset_group_name, genome_digest=genome_digest
        )
        seek_key = self.get_seek_key(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            seek_key_name=asset_registry_path_components.seek_key
            or self.get_default_seek_key(asset_group_name, asset_name, genome_digest=genome_digest),
        )

        if not is_path_type(seek_key.type):
            # Non-path types: return the value directly as a string
            return seek_key.value

        if seek_key.asset.path is None:
            raise ValueError(
                f"Seek key {seek_key.name} has no path defined for asset {seek_key.asset}"
            )

        # Choose the alias tree key. When the caller supplied a digest there is
        # no tree keyed by it; fall back to a local alias if the genome has one,
        # otherwise return the content (``data/``) path.
        alias_for_tree = supplied_alias
        if alias_for_tree is None and not abs_path:
            local_aliases = self._alias_manager.get_for_genome(genome_digest)
            alias_for_tree = local_aliases[0] if local_aliases else None

        if abs_path or alias_for_tree is None:
            seek_key_path = (self._genome_folder / seek_key.asset.path / seek_key.value).absolute()
        else:
            # Alias-tree path, named for the alias the caller supplied. The alias
            # tree rewrites the genome digest to the alias in filenames, so the
            # seek value is rewritten to match.
            rewritten_value = seek_key.value.replace(genome_digest, alias_for_tree)
            seek_key_path = (
                self._alias_folder / alias_for_tree / asset_group_name / asset_name / rewritten_value
            ).absolute()

        if force_exists and not seek_key_path.exists():
            raise FileNotFoundError(f"Seek key path not found: {seek_key_path}")
        return str(seek_key_path)

    def seek_remote(
        self,
        genome_name: str,
        asset_group_name: str,
        asset_name: str | None = None,
        seek_key: str | None = None,
        server_urls: list[str] | None = None,
    ) -> str:
        """
        Seek a remote path to an asset via the file-level endpoint.

        For assets with file serving mode, returns a direct URL to the file.
        For non-path seek keys, returns the value directly.

        Args:
            genome_name: The name of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            seek_key: The seek key.
            server_urls: Optional list of server URLs to query. If not provided,
                         uses all subscribed servers.

        Returns:
            str: The remote URL to the asset file, or the value directly for non-path seek keys.
        """
        return self._seek_remote_by_components(
            asset_registry_path_components=AssetRegistryPathComponents(
                genome=genome_name,
                asset_group=asset_group_name,
                asset=asset_name,
                seek_key=seek_key,
            ),
            server_urls=server_urls,
        )

    def _seek_remote_by_components(
        self,
        asset_registry_path_components: AssetRegistryPathComponents,
        server_urls: list[str] | None = None,
    ) -> str:
        """
        Get the remote path to an asset from registry path components.

        For non-path seek keys, the value is returned directly (no remote URL needed).
        For path-based seek keys, constructs a URL to the file-level endpoint.

        Args:
            asset_registry_path_components: The asset registry path components.
            server_urls: Optional list of server URLs to query. If not provided,
                         uses all subscribed servers.

        Returns:
            str: The remote URL to the asset file, or the value directly for non-path seek keys.

        Raises:
            ValueError: If the asset does not support file-level access.
        """
        # Check if this is a non-path seek key -- if so, return value directly
        try:
            seek_key = self.get_seek_key(
                genome_digest=self._alias_manager.resolve(asset_registry_path_components.genome),
                asset_group_name=asset_registry_path_components.asset_group,
                asset_name=(
                    asset_registry_path_components.asset
                    or self.get_default(
                        asset_registry_path_components.asset_group,
                        genome_name=asset_registry_path_components.genome,
                    )
                ),
                seek_key_name=asset_registry_path_components.seek_key
                or self.get_default_seek_key(
                    asset_registry_path_components.asset_group,
                    asset_registry_path_components.asset
                    or self.get_default(
                        asset_registry_path_components.asset_group,
                        genome_name=asset_registry_path_components.genome,
                    ),
                    genome_name=asset_registry_path_components.genome,
                ),
            )
            if not is_path_type(seek_key.type):
                return seek_key.value
        except Exception:
            pass  # Fall through to remote path resolution for path-based keys

        server_urls = server_urls or self._source_manager.get_subscriptions()
        if asset_registry_path_components.genome is None:
            raise ValueError(f"Genome name is required: {asset_registry_path_components=}")

        genome_token = asset_registry_path_components.genome
        for server_url in server_urls:
            client = self._source_manager.get_server_client(server_url)

            # Resolve the genome to a digest for THIS server. A local alias wins;
            # otherwise ask the server to resolve it as an alias (read-only -- we
            # write nothing locally); failing that, treat the token as a digest and
            # let the asset-group query below reject an unknown one.
            try:
                genome_digest = self._alias_manager.resolve(genome_token)
            except MissingAliasError:
                genome_digest = None
                try:
                    data = client.get(
                        operation_id=API_ID_ALIAS_DIGEST,
                        url_format_params={"name": genome_token},
                    )
                    if isinstance(data, dict):
                        genome_digest = data.get("digest")
                except Exception:
                    genome_digest = None
                if genome_digest is None:
                    genome_digest = genome_token

            # Look up asset group on server
            try:
                asset_groups = client.get_asset_groups(
                    params={
                        "genome_digest": genome_digest,
                        "asset_group_name": asset_registry_path_components.asset_group,
                    }
                )
            except Exception:
                continue
            if not asset_groups:
                continue

            # Fetch the group's assets from the server and choose one from the
            # response, rather than consulting the local DB (which need not have
            # this genome).
            asset_group_id = asset_groups[0]["id"]
            try:
                if asset_registry_path_components.asset is not None:
                    assets = client.get_assets(
                        params={
                            "asset_group_id": asset_group_id,
                            "name": asset_registry_path_components.asset,
                        }
                    )
                else:
                    assets = client.get_assets(params={"asset_group_id": asset_group_id})
            except Exception:
                continue
            if not assets:
                continue

            if asset_registry_path_components.asset is not None:
                asset_metadata = assets[0]
            else:
                # Default asset: the one flagged is_default, else the sole asset.
                default_assets = [a for a in assets if a.get("is_default") is True]
                if len(default_assets) == 1:
                    asset_metadata = default_assets[0]
                elif len(assets) == 1:
                    asset_metadata = assets[0]
                else:
                    continue  # ambiguous: no default and more than one asset

            serving_modes = asset_metadata.get("serving_modes", ["archive"])
            if "file" not in serving_modes:
                raise ValueError(
                    f"Asset '{asset_registry_path_components}' does not support file-level access "
                    f"(serving_modes={serving_modes}). Pull it locally with 'refgenie pull'."
                )

            # Choose the seek key from the server payload, not the local DB.
            seek_keys = asset_metadata.get("seek_keys") or []
            if asset_registry_path_components.seek_key is not None:
                chosen = next(
                    (
                        sk
                        for sk in seek_keys
                        if sk.get("name") == asset_registry_path_components.seek_key
                    ),
                    None,
                )
            else:
                chosen = self._default_seek_key_from_payload(
                    seek_keys, asset_registry_path_components.asset_group
                )
            if chosen is None:
                continue

            # Non-path seek keys return their value directly, no URL needed.
            if not _payload_is_path_type(chosen.get("type")):
                return chosen["value"]

            asset_digest = asset_metadata["digest"]
            file_path = chosen["value"]
            return f"{server_url}/v4/assets/{asset_digest}/files/{file_path}"

        raise ValueError("Failed to resolve remote seek path from any subscribed server.")

    @staticmethod
    def _default_seek_key_from_payload(
        seek_keys: list[dict], asset_group_name: str
    ) -> dict | None:
        """
        Pick the default seek key from a server asset payload's ``seek_keys``.

        Mirrors :meth:`get_default_seek_key`: prefer path-typed keys; if exactly
        one, use it; otherwise the key whose name equals the asset group name;
        else the first path key. With no path keys, fall back to the first seek
        key. Returns ``None`` when there are no seek keys at all.
        """
        if not seek_keys:
            return None
        path_keys = [sk for sk in seek_keys if _payload_is_path_type(sk.get("type"))]
        if not path_keys:
            return seek_keys[0]
        if len(path_keys) == 1:
            return path_keys[0]
        for sk in path_keys:
            if sk.get("name") == asset_group_name:
                return sk
        return path_keys[0]

    # === SeekKey/Default Operations ===

    def get_seek_key(
        self,
        asset_group_name: str,
        asset_name: str,
        seek_key_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
    ) -> SeekKey:
        """
        Get a seek key by its name.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            seek_key_name: The name of the seek key.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.

        Returns:
            SeekKey: The seek key.
        """
        genome_digest = self._resolve_genome_digest(genome_digest=genome_digest, genome_name=genome_name)
        statement = seek_key_by_name_stmt(
            genome_digest,
            asset_group_name,
            asset_name,
            seek_key_name,
            options=[selectinload(SeekKey.asset)],
        )
        with self._database_session as session:
            return one_or_raise(
                session,
                statement,
                MissingSeekKeyError(
                    genome=genome_digest,
                    asset_group=asset_group_name,
                    asset=asset_name,
                    seek_key=seek_key_name,
                ),
                unique=True,
            )

    def seek_key_exists(
        self,
        asset_group_name: str,
        asset_name: str,
        seek_key_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
    ) -> bool:
        """
        Check if a seek key exists.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            seek_key_name: The name of the seek key.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.

        Returns:
            bool: Whether the seek key exists.
        """
        genome_digest = self._resolve_genome_digest(genome_digest=genome_digest, genome_name=genome_name)
        statement = seek_key_by_name_stmt(
            genome_digest, asset_group_name, asset_name, seek_key_name
        )
        with self._database_session as session:
            result = session.exec(statement)
            return bool(result.first())

    def get_default_seek_key(
        self,
        asset_group_name: str,
        asset_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
    ) -> str:
        """
        Get the default seek key for an asset. Prefers path-based seek keys,
        since the default is used for symlinks, build inputs, etc.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.

        Returns:
            str: The default seek key name.
        """
        genome_digest = self._resolve_genome_digest(genome_digest=genome_digest, genome_name=genome_name)
        asset = self.get(
            genome_digest=genome_digest,
            genome_name=genome_name,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )
        # Filter to path-based seek keys for default selection
        path_seek_keys = [sk for sk in asset.seek_keys if is_path_type(sk.type)]
        if not path_seek_keys:
            # Fallback: if no path-based keys, return first seek key
            return asset.seek_keys[0].name
        if len(path_seek_keys) == 1:
            return path_seek_keys[0].name
        for seek_key in path_seek_keys:
            if seek_key.name == asset_group_name:
                return seek_key.name
        return path_seek_keys[0].name

    def list_seek_keys(
        self,
        genome_name: str,
        asset_group_name: str,
        asset_name: str | None = None,
    ) -> list[str]:
        """
        List the seek key names available for a given asset.

        If ``asset_name`` is not provided, the default asset for the asset
        group is used (matching the defaulting behavior of :meth:`seek`).

        Args:
            genome_name: The name of the genome (or alias).
            asset_group_name: The name of the asset group.
            asset_name: Optional name of the asset. Defaults to the asset
                group's default asset.

        Returns:
            list[str]: Seek key names defined on the resolved asset.
        """
        asset_name = asset_name or self.get_default(asset_group_name, genome_name=genome_name)
        asset = self.get(
            genome_name=genome_name,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )
        return [sk.name for sk in asset.seek_keys]

    def set_default(
        self,
        asset_group_name: str,
        asset_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
    ):
        """
        Set the default asset for an asset group.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.
        """
        genome_digest = self._resolve_genome_digest(genome_digest=genome_digest, genome_name=genome_name)
        asset_group = self.get_group(genome_digest=genome_digest, asset_group_name=asset_group_name)
        with self._database_session as session:
            # Look up the exact AssetName row for the name the caller passed,
            # scoped to this group. A name from another group is simply not found
            # here -- that is the same-group guard, by construction. The name the
            # caller chose is what gets flagged, so a non-canonical name is not
            # silently canonicalized.
            target = one_or_raise(
                session,
                select(AssetName).where(
                    AssetName.asset_group_id == asset_group.id,
                    AssetName.name == asset_name,
                ),
                MissingAssetError(
                    genome=genome_digest, asset_group=asset_group_name, asset=asset_name
                ),
            )
            # Clear the current default first, then set the new one. The partial
            # unique index forbids two defaults; a single transaction makes the
            # intermediate zero-default state invisible.
            session.exec(
                sa_update(AssetName)
                .where(
                    AssetName.asset_group_id == asset_group.id,
                    AssetName.is_default == True,  # noqa: E712 - SQL boolean, not Python identity
                )
                .values(is_default=False)
            )
            target.is_default = True
            session.add(target)
            session.commit()
        logger.info(f"Set default asset: '{genome_digest}/{asset_group_name}:{asset_name}'")

    def get_default(
        self,
        asset_group_name: str,
        *,
        genome_name: str | None = None,
        genome_digest: str | None = None,
    ) -> str | None:
        """
        Get the default asset for an asset group.

        Args:
            asset_group_name: The name of the asset group.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.

        Returns:
            str | None: The name of the default asset, or None if no default is set.
        """
        return self._read_default(
            asset_group_name=asset_group_name,
            genome_name=genome_name,
            genome_digest=genome_digest,
        )

    def _read_default(
        self,
        asset_group_name: str,
        *,
        genome_digest: str | None = None,
        genome_name: str | None = None,
    ):
        """
        Read default asset for an asset group from database.

        Args:
            asset_group_name: The name of the asset group.
            genome_digest: The digest of the genome.
            genome_name: The name of the genome.

        Returns:
            str: The name of the default asset, or None.
        """
        genome_digest = self._resolve_genome_digest(genome_digest=genome_digest, genome_name=genome_name)
        asset_group = self.get_group(genome_digest=genome_digest, asset_group_name=asset_group_name)
        with self._database_session as session:
            row = session.exec(
                select(AssetName).where(
                    AssetName.asset_group_id == asset_group.id,
                    AssetName.is_default == True,  # noqa: E712 - SQL boolean, not Python identity
                )
            ).one_or_none()
            return row.name if row else None

    @staticmethod
    def _bind_seek_key(
        asset: Asset,
        asset_class_seek_key: AssetClassSeekKey,
        directory_path: Path,
        session: Session,
        custom_seek_key_value: str | None = None,
    ) -> SeekKey:
        """
        Add a seek key to an asset.

        Args:
            asset: The asset.
            asset_class_seek_key: The seek key to add.
            directory_path: The directory path of the asset.
            session: The database session.
            custom_seek_key_value: For non-path seek key types, the resolved value
                to persist. Ignored for path-based seek key types.

        Returns:
            SeekKey: The added seek key.
        """
        if is_path_type(asset_class_seek_key.type):
            matched_file = asset_class_seek_key.match_file(directory_path=directory_path)
            value = (
                asset_class_seek_key.value
                if asset_class_seek_key.type == SeekKeyType.directory
                else matched_file.name
            )
        else:
            if custom_seek_key_value is None:
                raise ValueError(
                    f"Non-path seek key '{asset_class_seek_key.name}' requires a value "
                    f"from recipe custom_seek_keys"
                )
            value = custom_seek_key_value

        asset_seek_key = SeekKey(
            name=asset_class_seek_key.name,
            value=value,
            description=asset_class_seek_key.description,
            type=asset_class_seek_key.type,
            asset=asset,
        )
        session.add(asset_seek_key)
        return asset_seek_key

    @classmethod
    def _bind_asset_class_seek_keys(
        cls,
        asset: Asset,
        asset_class: Any,
        directory_path: Path,
        custom_seek_keys: dict[str, "str | tuple[str, SeekKeyType]"] | None,
        session: Session,
    ) -> None:
        """
        Bind every seek key of an asset class against a directory, plus any extras.

        Verifies that all files declared as seek keys by the asset class are present
        in the directory, persists them, and then persists custom seek keys that the
        asset class does not declare.

        Args:
            asset: The asset to bind the seek keys to.
            asset_class: The asset class declaring the seek keys.
            directory_path: The absolute directory path to resolve seek keys against.
            custom_seek_keys: Custom seek key values, keyed by seek key name.
            session: The database session.
        """
        for asset_class_seek_key in asset_class.seek_keys:
            raw = (custom_seek_keys or {}).get(asset_class_seek_key.name)
            csv = raw[0] if isinstance(raw, tuple) else raw
            cls._bind_seek_key(
                asset=asset,
                asset_class_seek_key=asset_class_seek_key,
                directory_path=directory_path,
                session=session,
                custom_seek_key_value=csv,
            )
        class_sk_names = {sk.name for sk in asset_class.seek_keys}
        cls._persist_extra_seek_keys(asset, class_sk_names, custom_seek_keys, session)

    @staticmethod
    def _persist_extra_seek_keys(
        asset: Asset,
        asset_class_seek_key_names: set[str],
        custom_seek_keys: dict[str, "str | tuple[str, SeekKeyType]"] | None,
        session: Session,
    ) -> None:
        """Persist custom seek keys not declared in the AssetClass."""
        for name, raw_value in (custom_seek_keys or {}).items():
            if name not in asset_class_seek_key_names:
                if isinstance(raw_value, tuple):
                    sk_value, sk_type = raw_value
                else:
                    sk_value, sk_type = raw_value, SeekKeyType.string
                session.add(
                    SeekKey(
                        name=name,
                        value=sk_value,
                        description=None,
                        type=sk_type,
                        asset=asset,
                    )
                )
