"""
SeekKeyManager - reading an asset's seek keys, plus the helpers that bind them.

Reached as ``rgc.asset.seek_key``. Its asset lookups go straight to
``queries.asset_by_name_stmt``, so it holds no reference to ``AssetManager``.

The module functions write seek-key rows into a caller's open session (the
content write path) or pick a default from a server payload (remote seek).
None of them use manager state.
"""

from pathlib import Path
from typing import Any, TYPE_CHECKING

from sqlalchemy.engine import Engine
from sqlalchemy.orm import selectinload
from sqlmodel import Session

from refgenie.db.tables import (
    Asset,
    AssetClassSeekKey,
    SeekKey,
    SeekKeyType,
    is_path_type,
)
from refgenie.exceptions import MissingAssetError, MissingSeekKeyError
from refgenie.managers.asset.queries import asset_by_name_stmt, seek_key_by_name_stmt
from refgenie.managers.base import ResourceManager
from refgenie.managers.queries import one_or_raise
from refgenie.models import GenomeDigest

if TYPE_CHECKING:
    from refgenie.managers.asset.group import AssetGroupManager


class SeekKeyManager(ResourceManager):
    """Manager for reading the seek keys of an asset."""

    def __init__(self, database_engine: Engine, groups: "AssetGroupManager"):
        """
        Initialize the SeekKeyManager.

        Args:
            database_engine: The database engine.
            groups: The AssetGroupManager, for a group's default asset.
        """
        super().__init__(database_engine)
        self._groups = groups

    def _asset(self, genome_digest: GenomeDigest, asset_group_name: str, asset_name: str) -> Asset:
        """The named asset, with its seek keys loaded."""
        statement = asset_by_name_stmt(
            genome_digest,
            asset_group_name,
            asset_name,
            options=[selectinload(Asset.seek_keys)],
        )
        with self._database_session as session:
            return one_or_raise(
                session,
                statement,
                MissingAssetError(
                    genome=genome_digest, asset_group=asset_group_name, asset=asset_name
                ),
                unique=True,
            )

    def get(
        self,
        asset_group_name: str,
        asset_name: str,
        seek_key_name: str,
        *,
        genome_digest: GenomeDigest,
    ) -> SeekKey:
        """
        Get a seek key by its name.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            seek_key_name: The name of the seek key.
            genome_digest: The genome digest.

        Returns:
            SeekKey: The seek key.
        """
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

    def exists(
        self,
        asset_group_name: str,
        asset_name: str,
        seek_key_name: str,
        *,
        genome_digest: GenomeDigest,
    ) -> bool:
        """
        Check if a seek key exists.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            seek_key_name: The name of the seek key.
            genome_digest: The genome digest.

        Returns:
            bool: Whether the seek key exists. False for an unknown genome.
        """
        statement = seek_key_by_name_stmt(
            genome_digest, asset_group_name, asset_name, seek_key_name
        )
        with self._database_session as session:
            result = session.exec(statement)
            return bool(result.first())

    def get_default(
        self,
        asset_group_name: str,
        asset_name: str,
        *,
        genome_digest: GenomeDigest,
    ) -> str:
        """
        Get the default seek key for an asset. Prefers path-based seek keys,
        since the default is used for symlinks, build inputs, etc.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
            genome_digest: The genome digest.

        Returns:
            str: The default seek key name.
        """
        asset = self._asset(genome_digest, asset_group_name, asset_name)
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

    def resolve(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str | None = None,
        seek_key_name: str | None = None,
    ) -> tuple[str, SeekKey]:
        """
        The seek key a registry path means, filling in whatever it leaves out.

        A missing asset name becomes the group's default asset, and a missing
        seek key name the asset's default seek key. Local seek and remote seek
        (for a local non-path value) both start here.

        Args:
            genome_digest: The genome digest.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset. Defaults to the group's default.
            seek_key_name: The name of the seek key. Defaults to the asset's default.

        Returns:
            (asset name, seek key). The name is returned too because the alias
            tree is keyed by the name asked for, which may not be the asset's
            canonical name.

        Raises:
            MissingAssetGroupError, MissingAssetError, MissingSeekKeyError: If
                the group, asset or seek key is not here.
        """
        asset_name = asset_name or self._groups.get_default(
            asset_group_name, genome_digest=genome_digest
        )
        seek_key = self.get(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            seek_key_name=seek_key_name
            or self.get_default(asset_group_name, asset_name, genome_digest=genome_digest),
        )
        return asset_name, seek_key

    def list_all(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str | None = None,
    ) -> list[str]:
        """
        List the seek key names available for a given asset.

        If ``asset_name`` is not provided, the default asset for the asset
        group is used (matching the defaulting behavior of ``AssetManager.seek``).

        Args:
            genome_digest: The genome digest.
            asset_group_name: The name of the asset group.
            asset_name: Optional name of the asset. Defaults to the asset
                group's default asset.

        Returns:
            list[str]: Seek key names defined on the resolved asset.
        """
        asset_name = asset_name or self._groups.get_default(
            asset_group_name, genome_digest=genome_digest
        )
        asset = self._asset(genome_digest, asset_group_name, asset_name)
        return [sk.name for sk in asset.seek_keys]


def payload_is_path_type(seek_key_type: "str | SeekKeyType | None") -> bool:
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


def default_seek_key_from_payload(seek_keys: list[dict], asset_group_name: str) -> dict | None:
    """
    Pick the default seek key from a server asset payload's ``seek_keys``.

    Mirrors :meth:`SeekKeyManager.get_default`: prefer path-typed keys; if
    exactly one, use it; otherwise the key whose name equals the asset group
    name; else the first path key. With no path keys, fall back to the first
    seek key. Returns ``None`` when there are no seek keys at all.
    """
    if not seek_keys:
        return None
    path_keys = [sk for sk in seek_keys if payload_is_path_type(sk.get("type"))]
    if not path_keys:
        return seek_keys[0]
    if len(path_keys) == 1:
        return path_keys[0]
    for sk in path_keys:
        if sk.get("name") == asset_group_name:
            return sk
    return path_keys[0]


def bind_seek_key(
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


def bind_asset_class_seek_keys(
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
        bind_seek_key(
            asset=asset,
            asset_class_seek_key=asset_class_seek_key,
            directory_path=directory_path,
            session=session,
            custom_seek_key_value=csv,
        )
    class_sk_names = {sk.name for sk in asset_class.seek_keys}
    persist_extra_seek_keys(asset, class_sk_names, custom_seek_keys, session)


def persist_extra_seek_keys(
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
