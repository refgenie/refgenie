"""Looper pre-submit hook: ``refgenie.integrations.looper.populate``.

Name it in a pipeline interface:

    pre_submit:
      python_functions:
      - refgenie.integrations.looper.populate

Templates then read ``{ refgenie[sample.genome].fasta.chrom_sizes }``, and
``refgenie://`` paths anywhere in the pipeline block are resolved. Looper calls
this once per sample; the Refgenie instance and its path view are built once
per process and reused, so each path is looked up at most once per run.

The database comes from ``var_templates.refgenie_db_config`` in the pipeline
interface when that is set and non-empty, else from ``$REFGENIE_DB_CONFIG_PATH``
(the ``Refgenie()`` default).
"""

from functools import cache
from pathlib import Path

from refgenie.core import Refgenie
from refgenie.core.paths import AssetPaths


def populate(namespaces: dict) -> dict:
    """Set ``namespaces["refgenie"]`` and resolve ``refgenie://`` in the pipeline block.

    Returns:
        dict: ``{"pipeline": ...}``, the pipeline block with every
        ``refgenie://`` path resolved, for looper to merge back in. A path to a
        missing asset is logged and left as written; an unknown genome raises.
    """
    pipeline = namespaces["pipeline"]
    db_config = (pipeline.get("var_templates") or {}).get("refgenie_db_config") or ""
    rg, paths = _open(db_config.strip() or None)
    # Set in place: looper merges returned dicts leaf by leaf, which would
    # force a full walk of the lazy view.
    namespaces["refgenie"] = paths
    return {"pipeline": rg.populate(dict(pipeline))}


@cache
def _open(db_config: str | None) -> tuple[Refgenie, AssetPaths]:
    rg = Refgenie(database_config_path=Path(db_config) if db_config else None)
    return rg, rg.paths()
