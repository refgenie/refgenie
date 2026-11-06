"""Shared test support: plain functions, classes and constants.

Fixtures and hooks live in ``tests/conftest.py``; everything importable lives
here, one module per concern. Top-level imports are limited to the stdlib,
``pytest`` and ``unittest.mock``. Imports of ``refgenie``, ``sqlmodel``,
``fastapi`` and friends go inside the functions that need them: the
``tests/e2e/cli_compat`` suite must stay library-free, and the mock-only
tests must run without the ``server`` extras.
"""
