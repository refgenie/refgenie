"""
Model Context Protocol server exposing read-only refgenie queries to LLM tools.

``stdio.py`` is the ``refgenie-mcp`` console script; the production server
mounts the same tools over Streamable HTTP from :mod:`refgenie.server.main`.
This package must not import ``fastapi`` or ``refgenie.server``
(enforced by ``tests/test_layering.py``).
"""
