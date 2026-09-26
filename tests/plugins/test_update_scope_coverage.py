"""Every public method that commits a local-state change carries ``@update_scope``.

A static check, in the style of ``tests/test_layering.py``: it reads the source,
so it needs nothing importable and costs nothing at runtime. It keeps a future
mutator from forgetting to fire ``post_update``. A method that commits but is
not a change to local assets, aliases or genomes is named below with a reason.
"""

import ast
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[2] / "refgenie"

#: (file, class) whose public committing methods must be scoped.
CHECKED = [
    ("managers/asset/manager.py", "AssetManager"),
    ("managers/asset/content.py", "AssetContentManager"),
    ("managers/asset/group.py", "AssetGroupManager"),
    ("managers/genome.py", "GenomeManager"),
]

#: "Class.method" -> why it commits without firing post_update.
NOT_A_LOCAL_STATE_CHANGE = {
    "AssetContentManager.add_incomplete": (
        "a placeholder row with no files; content.add records asset_added when it completes"
    ),
    "GenomeManager.apply_fhr_columns": "genome metadata only",
    "GenomeManager.set_store_name": "federation bookkeeping, not a local asset",
}


def _commits(func: ast.FunctionDef) -> bool:
    return any(
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and node.func.attr == "commit"
        for node in ast.walk(func)
    )


def _decorator_names(func: ast.FunctionDef) -> set[str]:
    names = set()
    for dec in func.decorator_list:
        if isinstance(dec, ast.Name):
            names.add(dec.id)
        elif isinstance(dec, ast.Attribute):
            names.add(dec.attr)
    return names


def _public_methods(path: Path, class_name: str) -> list[ast.FunctionDef]:
    tree = ast.parse(path.read_text())
    (cls,) = [n for n in tree.body if isinstance(n, ast.ClassDef) and n.name == class_name]
    return [n for n in cls.body if isinstance(n, ast.FunctionDef) and not n.name.startswith("_")]


def test_every_committing_mutator_is_scoped():
    missing = []
    for rel, class_name in CHECKED:
        for func in _public_methods(PACKAGE / rel, class_name):
            name = f"{class_name}.{func.name}"
            if not _commits(func) or name in NOT_A_LOCAL_STATE_CHANGE:
                continue
            if "update_scope" not in _decorator_names(func):
                missing.append(name)
    assert not missing, (
        "these methods commit a change but are not marked @update_scope, so no "
        f"post_update fires for them: {missing}. Mark them and record a Change "
        "after the commit, or name them in NOT_A_LOCAL_STATE_CHANGE with a reason."
    )


def test_the_allowlist_names_real_committing_methods():
    """A stale entry would hide nothing today and everything tomorrow."""
    committing = {
        f"{class_name}.{func.name}"
        for rel, class_name in CHECKED
        for func in _public_methods(PACKAGE / rel, class_name)
        if _commits(func)
    }
    assert set(NOT_A_LOCAL_STATE_CHANGE) <= committing
