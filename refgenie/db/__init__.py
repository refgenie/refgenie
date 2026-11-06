"""
The SQL catalog, and the only place SQL is written.

``tables.py`` defines the tables. ``events.py`` and ``cleanup.py`` remove the
files a delete makes obsolete, after the commit. ``migrations/`` is the alembic
environment; ``TARGET_ALEMBIC_VERSION`` in :mod:`refgenie.const` must name its
head. ``queries/`` holds every statement the managers run, one module per
entity. Nothing in ``db/`` imports a manager, the server, or the CLI.
"""
