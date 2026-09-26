"""
Model Context Protocol server exposing read-only refgenie queries to LLM tools.

``tools.py`` builds the ``MCPServer`` and registers the tools (list and search
genomes; list assets, asset classes and recipes; look up a digest; compare two
genomes). Each is a query against a lazily created :class:`Refgenie`, and
``set_refgenie`` lets a host hand in its own. ``stdio.py`` is the
``refgenie-mcp`` console script for local use over stdio; the production
server mounts the same tool set over Streamable HTTP from
:mod:`refgenie.server.main`.

The ``mcp`` dependency ships only in the ``mcp`` extra, so importing
``tools.py`` without it raises a clear ImportError. This package must stay
free of ``fastapi`` and ``refgenie.server``: it is installable without the
server stack, and ``tests/test_layering.py`` enforces that.
"""
