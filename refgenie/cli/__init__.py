"""
The ``refgenie`` command line: pydantic-settings models, argparse glue, and
the handlers that call into :class:`refgenie.core.Refgenie`.

``main.py`` is the console-script entry point: it preprocesses argv, builds
the parser from ``parser.py``, and hands the parsed command to the handler
looked up in ``dispatch.py``. Everything below ``main`` is imported lazily so
``--help`` and ``--version`` stay fast. ``framework.py`` holds the parser
machinery (``CliList``, the help renderer); ``messages.py`` the help strings;
``errors.py`` the exit codes and ``fail()``, the one way a handler reports
failure. The commands themselves live in ``commands/``, one module per family,
each holding both models and handlers; ``commands/helpers.py`` is what they
share. ``build_fasta.py`` is a separate console script
(``refgenie-build-fasta``) called from recipe shell templates.

This is the only package allowed to turn interactive prompts on (see
``main_cli``); the server and MCP entry points never read stdin.
"""
