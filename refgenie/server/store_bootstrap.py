"""Register federated stores on server boot.

Upserts the stores declared in ``REFGENIE_STORES`` (a JSON array of
``{"name", "url", "type"?, "priority"?, "description"?, "enabled"?}``) so a
fresh container needs no manual ``refgenie store add``. Only registry rows;
genomes come from the catalog import or ``refgenie store sync``.
"""

import json
import os

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
    rg.store.apply_declaration(declared)
