"""
Asset management: :class:`AssetManager` and the pieces it is assembled from.

``manager.py`` holds asset records: lookups, removal, local seek, listing and
rename. Only ``AssetManager`` is exported. Five smaller managers are composed
in as public attributes, none holding a reference back to ``AssetManager``:

- ``rgc.asset.content`` -- ``content.py`` (``AssetContentManager``): the
  content write path (``add``, ``adopt_name``), plus the
  stateless helpers it uses.
- ``rgc.asset.group`` -- ``group.py`` (``AssetGroupManager``): asset groups and
  the default asset of each group.
- ``rgc.asset.seek_key`` -- ``seek_key.py`` (``SeekKeyManager``): seek-key
  lookups and default filling (``resolve``), plus module functions that bind
  seek keys during a write.
- ``rgc.asset.tree`` -- ``alias_tree.py`` (``AliasTree``): rendering and
  purging the per-alias ``alias/`` and ``builds/`` trees on disk.
- ``rgc.asset.links`` -- ``links.py`` (``AssetLinkManager``): parent/child
  links between assets.

``tables.py`` holds the asset listing table, and ``colocation.py`` the
parent-file symlinks placed in child directories by build and pull.
"""

from refgenie.managers.asset.manager import AssetManager

__all__ = ["AssetManager"]
