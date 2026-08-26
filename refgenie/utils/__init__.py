"""
Leaf modules shared across refgenie's subpackages.

Import from the submodules directly (e.g. ``from refgenie.utils.build import
get_build_dir``); this package deliberately re-exports nothing, so that
``import refgenie.utils`` stays lightweight. Some of these modules encode
refgenie's on-disk conventions — ``build.py`` owns the asset-identity digest
scheme and the ``builds/`` layout, ``symlinks.py`` the alias-ownership rule,
``staging.py`` the staged-archive path convention. Nothing in ``utils/`` may
import ``refgenie.server``, ``refgenie.db``, or any other refgenie subpackage
above it in the dependency graph.
"""
