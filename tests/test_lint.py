"""
Package hygiene: the exception-hierarchy contract, and comment pointers.

Lint is a dedicated CI job (`ruff check .`), not a test.
"""

import inspect
import re
from pathlib import Path

#: Named in a docstring only as deleted history, never as a live pointer.
_DELETED_TEST_MODULES = {"tests/test_dash.py"}


def test_every_exception_is_exported_and_in_the_hierarchy():
    """`except RefgenieError` must catch everything the package raises, and
    `__all__` must name every exception it defines."""
    from refgenie import exceptions

    defined = {
        name
        for name, obj in vars(exceptions).items()
        if inspect.isclass(obj)
        and issubclass(obj, BaseException)
        and obj.__module__ == exceptions.__name__
    }

    assert defined - {"RefgenieError"} <= set(exceptions.__all__) - {"RefgenieError"}
    assert set(exceptions.__all__) == defined
    for name in defined - {"RefgenieError"}:
        assert issubclass(getattr(exceptions, name), exceptions.RefgenieError), name


def test_every_test_module_named_in_a_comment_exists():
    """A comment that points at a renamed or merged test module is a lie that
    ages into the tree; only `_DELETED_TEST_MODULES` may be named as history."""
    root = Path(__file__).resolve().parents[1]
    named = {
        match
        for path in (*root.glob("refgenie/**/*.py"), *root.glob("tests/**/*.py"))
        for match in re.findall(r"tests/test_[a-z_0-9]+\.py", path.read_text())
    }
    missing = sorted(m for m in named - _DELETED_TEST_MODULES if not (root / m).exists())
    assert not missing, f"comments name nonexistent test modules: {missing}"
