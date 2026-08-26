"""
The package's dependency direction, enforced.

`refgenie` is layered: configuration and leaf helpers at the bottom, the
database and domain models above them, then the managers, then the `Refgenie`
facade, and finally the four entry points (CLI, server, MCP, Snakefile) that
compose it. An import may point down a layer or sideways within one; it may
never point up. That rule is what lets `refgenie.utils` promise (in its own
docstring) that it drags nothing else in, and what keeps `refgenie.core` free
of the import cycle that used to force method-local imports in
`core/populate.py`.

The check reads import statements statically, so it costs nothing at runtime
and does not depend on any module being importable.
"""

import ast
from collections import defaultdict
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1] / "refgenie"

#: Component -> layer. A component is a subpackage of `refgenie` or a module at
#: its root. Higher numbers sit higher in the graph and may import lower ones.
LAYERS = {
    "config": 0,
    # Leaves: no refgenie imports at all, beyond `config` for the logger.
    "logger": 1,
    "const": 1,
    "exceptions": 1,
    "progress": 1,
    "utils": 2,
    "db": 3,
    "models": 4,
    "catalog_transfer": 4,
    "managers": 5,
    "core": 6,
    # Entry points and the root facade that re-exports for them.
    "refgenie": 7,
    "populator": 7,
    "cli": 7,
    "server": 7,
    "mcp": 7,
    "snakefile": 7,
}

#: `populator.py` is the looper integration surface: a dotted path that appears
#: verbatim in users' looper configs. It is a wrapper over the facade and must
#: stay one -- the registry-path implementation belongs in `core/populate.py`.
POPULATOR_MAY_IMPORT = {"core"}


def _is_module(dotted: str) -> bool:
    """Whether a dotted `refgenie...` name is a module or package in the tree."""
    parts = dotted.split(".")
    if parts[0] != "refgenie" or len(parts) == 1:
        return False
    base = PACKAGE.joinpath(*parts[1:])
    return base.is_dir() or base.with_suffix(".py").is_file()


def _component(module: str) -> str | None:
    """The layered component a dotted `refgenie...` module belongs to."""
    parts = module.split(".")
    if parts[0] != "refgenie":
        return None
    if len(parts) == 1:
        return "refgenie"
    head = parts[1]
    if (PACKAGE / head).is_dir() or (PACKAGE / f"{head}.py").is_file():
        return head
    return "refgenie"


def _imports(path: Path, module: str) -> set[str]:
    """Every refgenie component the module at `path` imports, itself excluded.

    `from refgenie import progress` names the root package but reaches a root
    module, so each imported name is tried as a submodule of the target before
    the target itself is charged with the dependency.
    """
    package = module if path.name == "__init__.py" else module.rsplit(".", 1)[0]
    found = set()
    for node in ast.walk(ast.parse(path.read_text())):
        if isinstance(node, ast.Import):
            targets = [alias.name for alias in node.names]
        elif isinstance(node, ast.ImportFrom):
            if node.level:
                base = package.split(".")[: -(node.level - 1) or None]
                target = ".".join(base + ([node.module] if node.module else []))
            else:
                target = node.module or ""
            targets = [
                submodule if _is_module(submodule := f"{target}.{alias.name}") else target
                for alias in node.names
            ]
        else:
            continue
        for target in targets:
            component = _component(target)
            if component is not None:
                found.add(component)
    return found - {_component(module)}


def _graph() -> dict[str, set[str]]:
    edges: dict[str, set[str]] = defaultdict(set)
    for path in sorted(PACKAGE.rglob("*.py")):
        if "migrations" in path.parts or "__pycache__" in path.parts:
            continue
        parts = path.relative_to(PACKAGE.parent).with_suffix("").parts
        module = ".".join(parts[:-1] if path.name == "__init__.py" else parts)
        source = _component(module)
        for target in _imports(path, module):
            edges[source].add(target)
    return edges


def test_every_component_has_a_declared_layer():
    """A new subpackage or root module must be placed in the graph, not left
    to slip past the check unnoticed."""
    present = {
        path.stem if path.is_file() else path.name
        for path in PACKAGE.iterdir()
        if (path.is_dir() and (path / "__init__.py").exists())
        or (path.is_file() and path.suffix == ".py" and path.name != "__init__.py")
    }
    assert present <= set(LAYERS), f"components with no declared layer: {present - set(LAYERS)}"


def test_imports_never_point_up_a_layer():
    violations = [
        f"{source} -> {target}"
        for source, targets in _graph().items()
        for target in targets
        if LAYERS[target] > LAYERS[source]
    ]
    assert not violations, "imports pointing up the dependency graph: " + ", ".join(
        sorted(violations)
    )


def test_utils_imports_nothing_above_itself():
    """`refgenie/utils/__init__.py` promises this in prose; hold it to it."""
    assert _graph()["utils"] <= {"config", "logger", "const", "exceptions", "progress"}


def test_populator_stays_a_thin_wrapper_over_the_facade():
    assert _graph()["populator"] <= POPULATOR_MAY_IMPORT


def test_models_names_one_thing():
    """`models.py` means the domain types, and only `refgenie/models.py` is it.

    Repeated basenames are fine here -- `manager.py`, `const.py`, `queries.py`
    and others are deliberately reused and scoped by their package. `models` is
    the exception because the word is ambiguous: domain types, ORM tables, or
    the JSON a client posts. Wire models go in a `schemas.py`; an internal
    record that is not on the wire stays with the code that owns it. See
    "Module naming" in `docs/design-notes.md`.

    This drifted once: a rename pass settled on `schemas.py`, then the commit
    adding the jobs and actions packages brought back two `models.py` files.
    """
    found = sorted(
        str(path.relative_to(PACKAGE.parent))
        for path in PACKAGE.rglob("models.py")
        if "__pycache__" not in path.parts
    )
    assert found == ["refgenie/models.py"], (
        "models.py must name the domain types and appear once; wire models "
        f"belong in a schemas.py. Found: {found}"
    )
