"""
Move catalog *metadata* from a build node to a server, with S3 as the carrier.

The build node exports a "publish catalog" (a fresh SQLite file holding only
what a server needs to list and serve the pushed assets) and uploads it next to
the asset files. The server imports it at startup and on a schedule; this
import is the only writer of the server's catalog besides download counters.

Export makes two changes:

- Aliases from both backends (the RefgetStore and the SQL ``alias`` table) are
  written as SQL rows, since a server reads only the SQL table.
- Pushed links are re-targeted at one https remote for the public mirror of the
  stage folder, since download redirects only use http(s) remotes.
"""

import logging
import tempfile
from contextlib import contextmanager
from pathlib import Path

from refgenie.db import publish
from refgenie.db.tables import RemoteType

logger = logging.getLogger(__name__)


def _assert_no_dropped_aliases(sql_aliases: list[dict], genome_digests, aliases) -> None:
    """Fail the export if it would publish fewer names than the SQL alias table holds.

    Raises:
        RuntimeError: If any SQL alias row for an exported genome has no
            counterpart in the aliases about to be written.
    """
    sql_names = {r["name"] for r in sql_aliases if r["genome_digest"] in genome_digests}
    dropped = sql_names - {a["name"] for a in aliases}
    if dropped:
        sample = ", ".join(sorted(dropped)[:10])
        raise RuntimeError(
            f"Export would drop {len(dropped)} of {len(sql_names)} alias rows the "
            f"catalog holds for exported genomes (e.g. {sample}). The alias "
            f"manager is not reading the whole alias space; the exported genomes "
            f"would be served without names."
        )


def export_publish_catalog(rg, dest_path: str | Path, https_prefix: str) -> dict[str, int]:
    """
    Write the publishable subset of ``rg``'s catalog to a fresh SQLite file.

    Only pushed (asset_digest, mode) pairs are included, re-linked to one https
    remote at ``https_prefix``. Local paths are stripped and download counters
    start at zero.

    Returns a table -> row-count summary.
    """
    dest_path = Path(dest_path)
    if dest_path.exists():
        dest_path.unlink()
    dest_path.parent.mkdir(parents=True, exist_ok=True)

    pushed_pairs, rows = publish.read_publishable(rg.database_engine)
    genome_digests = {g["digest"] for g in rows["genome"]}

    # The alias manager unions the RefgetStore and SQL aliases.
    aliases = [
        {"name": a.name, "genome_digest": a.genome_digest}
        for a in rg.alias.list_all()
        if a.genome_digest in genome_digests
    ]
    _assert_no_dropped_aliases(rows["alias"], genome_digests, aliases)

    for a in rows["asset"]:
        a["path"] = None  # build-node-local; meaningless (and ugly) in public

    config_rows = sorted(rows["configuration"], key=lambda c: c["id"])
    config_rows = config_rows[-1:]  # newest only; its stage folder is the strip-prefix
    for c in config_rows:
        c["plugin_settings"] = {}  # build-node-local (paths); never published

    staged = [s for s in rows["stagedasset"] if (s["asset_digest"], s["mode"]) in pushed_pairs]
    for s in staged:
        s["download_count"] = 0  # server-owned counter, not build state

    published = {digest for digest, _ in pushed_pairs}
    remote_row = {
        "id": 1,
        "prefix": https_prefix,
        "type": RemoteType.https,
        "name": "Public https mirror of the build stage folder",
        "push_command": None,
        "configuration_id": config_rows[0]["id"] if config_rows else None,
    }

    table_rows: dict[str, list[dict]] = {
        **rows,
        "configuration": config_rows,
        "remote": [remote_row],
        "alias": aliases,
        "assetlink": [
            link
            for link in rows["assetlink"]
            if link["parent_digest"] in published and link["child_digest"] in published
        ],
        "stagedasset": staged,
        "remoteassetlink": [
            {"remote_id": 1, "asset_digest": digest, "mode": mode, "pushed": True}
            for digest, mode in sorted(pushed_pairs)
        ],
    }
    summary = publish.write_catalog(dest_path, table_rows)
    logger.info(f"Exported publish catalog to {dest_path}: {summary}")
    return summary


@contextmanager
def _local_copy(source: str):
    """Yield a local path for ``source``, downloading it first if it is a URL."""
    if not source.startswith(("http://", "https://")):
        yield source
        return
    import httpx

    with tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False) as tmp:
        with httpx.stream("GET", source, follow_redirects=True, timeout=60) as resp:
            resp.raise_for_status()
            for chunk in resp.iter_bytes():
                tmp.write(chunk)
    try:
        yield tmp.name
    finally:
        Path(tmp.name).unlink(missing_ok=True)


def import_publish_catalog(engine, source: str) -> dict[str, int]:
    """
    Upsert a publish-catalog artifact (local path or http(s) URL) into
    ``engine``'s catalog.

    Insert-or-update in one transaction. Nothing is deleted (ORM deletes would
    queue file removals), except that each species taxon in the artifact is
    rewritten. Retracted rows stay until :func:`prune_publish_catalog` runs.
    The artifact's schema stamp must match this code's.
    """
    with _local_copy(source) as local_path:
        summary = publish.import_catalog(engine, local_path)

    logger.info(f"Imported publish catalog from {source}: {summary}")
    return summary


def prune_publish_catalog(engine, source: str, dry_run: bool = False) -> dict[str, list[str]]:
    """
    Delete from ``engine`` the assets that the artifact at ``source`` retracted.

    Never run automatically. See :func:`refgenie.db.publish.prune_catalog`.
    """
    with _local_copy(source) as local_path:
        removed = publish.prune_catalog(engine, local_path, dry_run=dry_run)
    verb = "Would prune" if dry_run else "Pruned"
    logger.info(
        f"{verb} {len(removed['asset'])} asset(s) and {len(removed['assetgroup'])} "
        f"asset group(s) retracted by {source}"
    )
    return removed
