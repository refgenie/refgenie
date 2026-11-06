"""
The ``refgenie`` command line. ``main.py`` is the entry point; each command
family lives in ``commands/`` with its models and handlers.

This is the only package allowed to turn interactive prompts on (see
``main_cli``); the server and MCP entry points never read stdin.
"""
