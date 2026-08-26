"""Register federated stores on server boot.

The ``store`` registry is the single source of truth for what the server
federates over, but a fresh deployment starts with an empty table. This module
idempotently upserts a declared set of stores at startup so a container comes up
serving the right stores with no manual ``refgenie store add``.

The declaration is a JSON array in the ``REFGENIE_STORES`` environment variable,
each entry ``{"name", "url", "type"?, "priority"?, "description"?, "enabled"?}``.
Genome/alias rows still come from the published catalog import (or an explicit
``refgenie store sync``); this only manages the registry rows.
"""

import json
import os

from refgenie.db.tables import Store, StoreType
from refgenie.logger import logger


def bootstrap_stores_from_env(rg) -> None:
    """Apply the ``REFGENIE_STORES`` declaration to ``rg``'s store registry."""
    raw = os.environ.get("REFGENIE_STORES")
    if not raw:
        return
    try:
        declared = json.loads(raw)
    except json.JSONDecodeError as e:
        logger.error(f"REFGENIE_STORES is not valid JSON; ignoring it: {e}")
        return
    if not isinstance(declared, list):
        logger.error("REFGENIE_STORES must be a JSON array of store objects; ignoring it.")
        return
    apply_store_declaration(rg, declared)


def apply_store_declaration(rg, declared: list[dict]) -> None:
    """Idempotently upsert each declared store into the registry.

    Insert if new; update url/type/priority/enabled/description when they drift.
    Never removes rows not in the declaration -- a store added out of band stays.
    """
    from sqlmodel import Session, select

    with Session(rg.database_engine) as session:
        for entry in declared:
            name = entry.get("name")
            url = entry.get("url")
            if not name or not url:
                logger.error(f"Skipping store declaration with no name/url: {entry!r}")
                continue
            store_type = StoreType(entry.get("type", "remote"))
            priority = int(entry.get("priority", 100))
            enabled = bool(entry.get("enabled", True))
            description = entry.get("description")

            existing = session.exec(select(Store).where(Store.name == name)).one_or_none()
            if existing is None:
                session.add(
                    Store(
                        name=name,
                        url=url,
                        type=store_type,
                        priority=priority,
                        enabled=enabled,
                        description=description,
                    )
                )
                logger.info(f"Registered store '{name}' ({store_type.value}, priority={priority})")
            else:
                existing.url = url
                existing.type = store_type
                existing.priority = priority
                existing.enabled = enabled
                existing.description = description
                logger.info(f"Updated store '{name}' from REFGENIE_STORES")
        session.commit()
