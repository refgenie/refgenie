"""
Asset management: :class:`AssetManager` and the pieces it is assembled from.

``manager.py`` holds asset records: lookups, removal, local seek, listing and
rename. Only ``AssetManager`` is exported. It inherits nothing but
``ResourceManager``; five smaller managers are composed in as public
attributes, each taking its dependencies in its constructor and none holding a
reference back to ``AssetManager``:

- ``rgc.asset.content`` -- ``content.py`` (``AssetContentManager``): the
  content write path (``add``, ``adopt_name``, ``add_incomplete``), plus the
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

Callers use these attributes directly; ``AssetManager`` does not wrap them.
``queries.py`` holds the shared SELECT builders for ``genome/group:name`` lookups,
``tables.py`` the asset listing table, and ``colocation.py`` the parent-file
symlinks that the build manager and the puller place in child directories.

Pulling is ``rgc.transfer`` (``refgenie.managers.transfer``) and building is
``rgc.build`` (``refgenie.managers.build.BuildManager``). Both depend on
``rgc.asset`` one way.
What a server has, and remote seek, is ``rgc.servers``
(``refgenie.managers.sources.servers``).
"""

from refgenie.managers.asset.manager import AssetManager

__all__ = ["AssetManager"]
