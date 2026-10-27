"""Refgenie's side of outside workflow tools (looper, Snakemake).

These are bridges, not front doors like ``cli/``, ``server/`` and ``mcp/``.
Each module is imported by a dotted path a user types into the outside tool's
config (``refgenie.integrations.looper.populate``) or by one CLI command
(``refgenie generate snakefile``). Nothing here imports the outside tool
itself, so refgenie carries no dependency on it. This package re-exports
nothing; import from the submodules.
"""
