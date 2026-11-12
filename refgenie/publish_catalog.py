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

It also publishes the store list (enabled remote stores only), so a server
learns which stores to read from the catalog. Genomes in a disabled store are
left out. Genomes the build node made itself have no
store; when a *catalog store* is named, they are exported as living in it.
"""

import logging
import tempfile
from contextlib import contextmanager
from pathlib import Path

from refgenie.db import publish
from refgenie.db.tables import RemoteType, Store, StoreType
from refgenie.models import GenomeAlias

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


def _publishable_stores(
    store_rows: list[dict],
    catalog_store_name: str | None,
    catalog_store_url: str | None,
    catalog_store_priority: int,
) -> list[dict]:
    """The store rows to publish: the remote ones, plus the catalog's own store.

    An ``on_disk`` row names a build-node path, so it is never published. The
    catalog store replaces any build-catalog row of the same name.

    Raises:
        ValueError: If only one of the catalog store's name and url is given.
    """
    if bool(catalog_store_name) != bool(catalog_store_url):
        raise ValueError("A catalog store needs both a name and a url")
    stores = [
        s
        for s in store_rows
        if s["type"] == StoreType.remote and s["enabled"] and s["name"] != catalog_store_name
    ]
    skipped = sorted(s["name"] for s in store_rows if s["type"] != StoreType.remote)
    if skipped:
        logger.info(f"Not publishing on-disk store(s): {', '.join(skipped)}")
    if catalog_store_name:
        stores.append(
            Store(
                name=catalog_store_name,
                url=catalog_store_url,
                type=StoreType.remote,
                priority=catalog_store_priority,
                enabled=True,
                description="Genomes built for this catalog",
            ).model_dump()
        )
    return stores


def _omit_disabled_store_genomes(
    rows: dict[str, list[dict]], pushed_pairs: set[tuple[str, str]]
) -> tuple[set[tuple[str, str]], dict[str, int]]:
    """Drop, in place, every genome in a disabled store and everything hanging off it.

    A disabled store is not published, so neither are its genomes, their SQL
    aliases, asset groups, assets, asset names, seek keys, links and staged
    records. Returns the pushed pairs that remain and ``{store: genomes omitted}``.
    """
    disabled = {s["name"] for s in rows["store"] if not s["enabled"]}
    omitted: dict[str, int] = {}
    dropped: set[str] = set()
    for g in rows["genome"]:
        if g["store_name"] in disabled:
            dropped.add(g["digest"])
            omitted[g["store_name"]] = omitted.get(g["store_name"], 0) + 1
    if not dropped:
        return pushed_pairs, omitted
    rows["genome"] = [g for g in rows["genome"] if g["digest"] not in dropped]
    rows["alias"] = [a for a in rows["alias"] if a["genome_digest"] not in dropped]
    group_ids = {g["id"] for g in rows["assetgroup"] if g["genome_digest"] in dropped}
    rows["assetgroup"] = [g for g in rows["assetgroup"] if g["id"] not in group_ids]
    assets = {a["digest"] for a in rows["asset"] if a["asset_group_id"] in group_ids}
    rows["asset"] = [a for a in rows["asset"] if a["digest"] not in assets]
    rows["assetname"] = [n for n in rows["assetname"] if n["asset_group_id"] not in group_ids]
    rows["seekkey"] = [k for k in rows["seekkey"] if k["asset_digest"] not in assets]
    rows["stagedasset"] = [s for s in rows["stagedasset"] if s["asset_digest"] not in assets]
    rows["assetlink"] = [
        link
        for link in rows["assetlink"]
        if link["parent_digest"] not in assets and link["child_digest"] not in assets
    ]
    for name, count in sorted(omitted.items()):
        logger.info(f"Not publishing disabled store '{name}': omitting its {count} genome(s)")
    return {(d, m) for d, m in pushed_pairs if d not in assets}, omitted


def _assert_store_names_published(stores: list[dict], **tables: list[dict]) -> None:
    """Fail the export if any exported row names a store it does not publish.

    A server import of such a catalog would fail on the store foreign key.

    Raises:
        RuntimeError: Naming each unpublished store and the rows pointing at it.
    """
    published = {s["name"] for s in stores}
    missing: dict[str, dict[str, int]] = {}
    for table, rows in tables.items():
        for r in rows:
            name = r.get("store_name")
            if name is not None and name not in published:
                counts = missing.setdefault(name, {})
                counts[table] = counts.get(table, 0) + 1
    if missing:
        detail = "; ".join(
            f"'{name}' ({', '.join(f'{n} {t}' for t, n in sorted(c.items()))})"
            for name, c in sorted(missing.items())
        )
        raise RuntimeError(
            f"Exported rows name stores the catalog does not publish: {detail}. "
            f"Only remote stores are published; register them with "
            f"'refgenie store apply', or name the catalog store."
        )


def export_publish_catalog(
    rg,
    dest_path: str | Path,
    https_prefix: str,
    catalog_store_name: str | None = None,
    catalog_store_url: str | None = None,
    catalog_store_priority: int = 30,
) -> dict[str, int]:
    """
    Write the publishable subset of ``rg``'s catalog to a fresh SQLite file.

    Only pushed (asset_digest, mode) pairs are included, re-linked to one https
    remote at ``https_prefix``. Local paths are stripped and download counters
    start at zero. Enabled remote stores are published; on-disk ones are not.
    A disabled store is left out along with its genomes and everything that
    hangs off them (aliases, asset groups, assets, staged records).

    Args:
        rg: The build node's Refgenie instance.
        dest_path: Where to write the artifact.
        https_prefix: Public https base URL mirroring the stage folder.
        catalog_store_name: The public store holding the genomes this node
            built itself (those with no store). When given, it is published
            and those genomes are exported as living in it.
        catalog_store_url: That store's public root URL. Required with the name.
        catalog_store_priority: That store's priority (lower wins alias ties).

    Returns a table -> row-count summary.

    Raises:
        ValueError: If only one of the catalog store's name and url is given.
        RuntimeError: If an exported genome or alias names a store the
            catalog does not publish.
    """
    pushed_pairs, rows = publish.read_publishable(rg.database_engine)
    pushed_pairs, _ = _omit_disabled_store_genomes(rows, pushed_pairs)
    stores = _publishable_stores(
        rows["store"], catalog_store_name, catalog_store_url, catalog_store_priority
    )
    if catalog_store_name:
        for g in rows["genome"]:
            if g["store_name"] is None:
                g["store_name"] = catalog_store_name
    genome_digests = {g["digest"] for g in rows["genome"]}

    # The alias manager unions the RefgetStore and SQL aliases.
    aliases = [
        {"name": a.name, "genome_digest": a.genome_digest}
        for a in rg.alias.list_all()
        if a.genome_digest in genome_digests
    ]
    _assert_no_dropped_aliases(rows["alias"], genome_digests, aliases)
    _assert_store_names_published(stores, genome=rows["genome"], alias=aliases)

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
        "store": stores,
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
    dest_path = Path(dest_path)
    if dest_path.exists():
        dest_path.unlink()
    dest_path.parent.mkdir(parents=True, exist_ok=True)
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
    queue file removals), except that each species taxon and each carried
    genome's alias set is rewritten; an alias with '/' is skipped. Retracted
    assets, genomes and stores stay until :func:`prune_publish_catalog` runs.
    The artifact's schema stamp must match this code's.
    """
    with _local_copy(source) as local_path:
        summary = publish.import_catalog(engine, local_path, valid_alias=GenomeAlias.is_valid)

    logger.info(f"Imported publish catalog from {source}: {summary}")
    return summary


def prune_publish_catalog(engine, source: str, dry_run: bool = False) -> dict[str, list[str]]:
    """
    Delete from ``engine`` what the artifact at ``source`` no longer has: retracted
    assets, genomes it dropped that have no assets, and empty disabled stores.

    Never run automatically. See :func:`refgenie.db.publish.prune_catalog`.
    """
    with _local_copy(source) as local_path:
        removed = publish.prune_catalog(engine, local_path, dry_run=dry_run)
    verb = "Would prune" if dry_run else "Pruned"
    logger.info(
        f"{verb} {len(removed['asset'])} asset(s), {len(removed['assetgroup'])} asset "
        f"group(s), {len(removed['genome'])} genome(s), {len(removed['alias'])} alias(es) "
        f"and {len(removed['store'])} store(s) no longer in {source}"
    )
    if removed["kept_genome"]:
        logger.warning(
            f"Kept {len(removed['kept_genome'])} genome(s) the catalog no longer carries, "
            f"because they still have assets"
        )
    return removed
