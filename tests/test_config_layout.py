"""
`refgenie/config/__init__.py` indexes; it does not define.

The rule for a package `__init__.py`: a docstring, imports that re-export, and
`__all__`. Every definition -- an env-var read, a `mkdir`, the process-wide
`config` object, a warning printed at import -- lives in a named module
(`config/env.py`, `config/settings.py`) so that reading or reloading it is a
deliberate act on a file you can point at, not a side effect of importing the
package.

The check reads the file statically, so it costs nothing at runtime.
"""

import ast
from pathlib import Path

import pytest

CONFIG_INIT = Path(__file__).resolve().parents[1] / "refgenie" / "config" / "__init__.py"


def _offending_statements(tree: ast.Module) -> list[str]:
    """Top-level statements that are neither a docstring, an import, nor `__all__`."""
    offenders = []
    for index, node in enumerate(tree.body):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            continue
        if (
            index == 0
            and isinstance(node, ast.Expr)
            and isinstance(node.value, ast.Constant)
            and isinstance(node.value.value, str)
        ):
            continue  # the module docstring
        if (
            isinstance(node, ast.Assign)
            and len(node.targets) == 1
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id == "__all__"
        ):
            continue
        offenders.append(f"line {node.lineno}: {ast.unparse(node).splitlines()[0]}")
    return offenders


@pytest.mark.unit
def test_config_init_only_reexports():
    tree = ast.parse(CONFIG_INIT.read_text(), filename=str(CONFIG_INIT))
    offenders = _offending_statements(tree)
    assert not offenders, (
        "refgenie/config/__init__.py must hold only a docstring, imports and "
        "__all__; move these into refgenie/config/env.py or settings.py:\n  "
        + "\n  ".join(offenders)
    )


@pytest.mark.unit
def test_config_init_declares_all():
    tree = ast.parse(CONFIG_INIT.read_text(), filename=str(CONFIG_INIT))
    has_all = any(
        isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "__all__" for t in node.targets)
        for node in tree.body
    )
    assert has_all, "refgenie/config/__init__.py must declare an explicit __all__"
