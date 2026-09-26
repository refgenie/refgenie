"""
The SQL catalog: SQLModel tables, their filesystem side effects, and the
alembic migrations that keep the schema current.

``tables.py`` defines every table (``Genome``, ``Alias``, ``Asset``,
``AssetGroup``, ``AssetClass``, ``Recipe``, ``Store``, ``StagedAsset``,
``Configuration``, ...). ``events.py`` registers the SQLAlchemy listeners that
resolve which files a delete makes obsolete, and ``cleanup.py`` is the queue
that removes them only after the transaction commits: the catalog commits
first, and the filesystem follows. ``migrations/`` is the alembic environment;
``migrations/utils.py`` runs it against a database URL, and
``TARGET_ALEMBIC_VERSION`` in :mod:`refgenie.const` must name the chain head.

Engine creation and the session lifecycle live in :mod:`refgenie.core`, not
here, and the ``Refgenie`` root is what calls ``register_events``. Nothing in
``db/`` imports a manager, the server, or the CLI; ``refgenie.utils`` and
``refgenie.logger`` are its only refgenie dependencies.
"""
