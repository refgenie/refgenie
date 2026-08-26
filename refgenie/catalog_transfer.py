"""
Move catalog *metadata* from a build node to a server, with S3 as the carrier.

The build node exports a filtered "publish catalog" -- a fresh SQLite file in
the current schema holding only what a server needs to list and serve the
published assets -- and uploads it next to the asset files. The server imports
that artifact at startup and on a schedule. Files never travel this path; they
are already on S3. Metadata never travels any other path; the server's SQL
catalog is written only by this import (plus its own download counters).

Two deliberate transformations happen at export, both consequences of how the
server reads its catalog:

- The build node splits its aliases across two backends: the RefgetStore holds
  the names of the genomes it built, and the SQL ``alias`` table holds the names
  of the genomes it federates over. A server reads the SQL table alone, so both
  halves are materialized into rows here. Copying tables alone would leave every
  locally built genome nameless.
- Download redirects only consider remotes of type https/http
  (``server/remote_assets.py``), while the build node pushes through an
  s3-type remote. The export therefore emits a single https remote pointing at
  the public mirror of the stage folder, and re-targets the pushed links at it.
"""

import logging
import tempfile
from pathlib import Path

from sqlalchemy import create_engine, select
from sqlmodel import SQLModel

from refgenie.const import TARGET_ALEMBIC_VERSION
from refgenie.db.tables import (
    AlembicVersion,
    Alias,
    Asset,
    AssetClass,
    AssetClassLink,
    AssetClassSeekKey,
    AssetGroup,
    AssetLink,
    AssetName,
    Configuration,
    Genome,
    Recipe,
    RecipeAssetClassesInputs,
    Remote,
    RemoteAssetLink,
    RemoteType,
    SeekKey,
    StagedAsset,
)

logger = logging.getLogger(__name__)

# FK-safe order for both export inserts and import upserts.
_ORDERED_MODELS = [
    Configuration,
    Remote,
    Genome,
    Alias,
    AssetClass,
    AssetClassSeekKey,
    AssetClassLink,
    Recipe,
    RecipeAssetClassesInputs,
    AssetGroup,
    Asset,
    AssetName,
    AssetLink,
    SeekKey,
    StagedAsset,
    RemoteAssetLink,
]


def _rows(conn, table, where=None):
    """All rows of a table as plain dicts, optionally filtered."""
    stmt = select(table)
    if where is not None:
        stmt = stmt.where(where)
    return [dict(r._mapping) for r in conn.execute(stmt)]


def _assert_no_dropped_aliases(conn, alias_table, genome_digests, aliases) -> None:
    """Fail the export if it would publish fewer names than the catalog holds.

    Reading aliases through a backend that sees only one half of the alias
    space produces a catalog whose genomes are served as bare digests, and it
    does so silently: every table has rows, just not enough of them. So compare
    what is about to be written against the SQL table, which is the half an
    incomplete backend drops.

    Raises:
        RuntimeError: If any SQL alias row for an exported genome has no
            counterpart in the aliases about to be written.
    """
    sql_names = {
        r.name
        for r in conn.execute(select(alias_table.c.name, alias_table.c.genome_digest))
        if r.genome_digest in genome_digests
    }
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

    Only assets with at least one pushed RemoteAssetLink are included, and
    only their pushed (asset_digest, mode) pairs are re-linked -- against a
    single synthesized https remote at ``https_prefix`` -- so every row in the
    artifact describes something a client can actually download. Local build
    paths are stripped; download counters start at zero.

    Returns a table -> row-count summary.
    """
    dest_path = Path(dest_path)
    if dest_path.exists():
        dest_path.unlink()
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    dest_engine = create_engine(f"sqlite:///{dest_path}")
    SQLModel.metadata.create_all(dest_engine)

    t = {m.__tablename__: m.__table__ for m in _ORDERED_MODELS}
    summary: dict[str, int] = {}

    with rg.database_engine.connect() as src, dest_engine.begin() as dest:
        pushed_pairs = {
            (r.asset_digest, r.mode)
            for r in src.execute(
                select(
                    t["remoteassetlink"].c.asset_digest, t["remoteassetlink"].c.mode
                ).where(t["remoteassetlink"].c.pushed == True)  # noqa: E712
            )
        }
        published = {digest for digest, _ in pushed_pairs}

        genomes = _rows(src, t["genome"])
        genome_digests = {g["digest"] for g in genomes}

        # Aliases come from the mode-selected manager rather than the SQL table
        # because the build node keeps half of them in its RefgetStore. The
        # manager unions both halves; the SQL table alone would drop every
        # locally built name.
        aliases = [
            {"name": a.name, "genome_digest": a.genome_digest}
            for a in rg.alias.list_all()
            if a.genome_digest in genome_digests
        ]
        _assert_no_dropped_aliases(src, t["alias"], genome_digests, aliases)

        groups = _rows(src, t["assetgroup"])
        published_group_ids = {
            r.asset_group_id
            for r in src.execute(
                select(t["asset"].c.asset_group_id).where(t["asset"].c.digest.in_(published))
            )
        }
        groups = [g for g in groups if g["id"] in published_group_ids]

        assets = _rows(src, t["asset"], t["asset"].c.digest.in_(published))
        for a in assets:
            a["path"] = None  # build-node-local; meaningless (and ugly) in public

        config_rows = _rows(src, t["configuration"])
        config_rows.sort(key=lambda c: c["id"])
        config_rows = config_rows[-1:]  # newest only; its stage folder is the strip-prefix

        staged = [
            s
            for s in _rows(src, t["stagedasset"])
            if (s["asset_digest"], s["mode"]) in pushed_pairs
        ]
        for s in staged:
            s["download_count"] = 0  # server-owned counter, not build state

        remote_row = {
            "id": 1,
            "prefix": https_prefix,
            "type": RemoteType.https,
            "description": "Public https mirror of the build stage folder",
            "push_command": None,
            "configuration_id": config_rows[0]["id"] if config_rows else None,
        }

        table_rows: dict[str, list[dict]] = {
            "configuration": config_rows,
            "remote": [remote_row],
            "genome": genomes,
            "alias": aliases,
            "assetclass": _rows(src, t["assetclass"]),
            "assetclassseekkey": _rows(src, t["assetclassseekkey"]),
            "assetclasslink": _rows(src, t["assetclasslink"]),
            "recipe": _rows(src, t["recipe"]),
            "recipeassetclassesinputs": _rows(src, t["recipeassetclassesinputs"]),
            "assetgroup": groups,
            "asset": assets,
            "assetname": _rows(src, t["assetname"], t["assetname"].c.asset_digest.in_(published)),
            "assetlink": [
                link
                for link in _rows(src, t["assetlink"])
                if link["parent_digest"] in published and link["child_digest"] in published
            ],
            "seekkey": _rows(src, t["seekkey"], t["seekkey"].c.asset_digest.in_(published)),
            "stagedasset": staged,
            "remoteassetlink": [
                {"remote_id": 1, "asset_digest": digest, "mode": mode, "pushed": True}
                for digest, mode in sorted(pushed_pairs)
            ],
        }

        for model in _ORDERED_MODELS:
            name = model.__tablename__
            rows = table_rows[name]
            if rows:
                dest.execute(model.__table__.insert(), rows)
            summary[name] = len(rows)

        dest.execute(
            AlembicVersion.__table__.insert(), [{"version_num": TARGET_ALEMBIC_VERSION}]
        )

    dest_engine.dispose()
    logger.info(f"Exported publish catalog to {dest_path}: {summary}")
    return summary


# Tables whose surrogate id nothing references: match rows by natural key and
# never import the artifact's id, so re-imports can't collide with rows the
# server already assigned ids to.
# assetname is here because of unique_group_build_digest: a rebuilt build-node
# catalog renumbers assetname.id, so importing the artifact's id can land a
# build_digest the server already holds under a different id and fail the whole
# transaction on the index. Matching on the natural key never imports the id.
_NATURAL_KEYS = {
    "stagedasset": ("asset_digest", "mode"),
    "assetname": ("asset_group_id", "name"),
}
# Server-owned columns an import must never overwrite on an existing row.
_PRESERVE = {"stagedasset": ("download_count",)}


def import_publish_catalog(engine, source: str) -> dict[str, int]:
    """
    Upsert a publish-catalog artifact (local path or http(s) URL) into
    ``engine``'s catalog.

    Insert-or-update only, in one transaction; nothing is ever deleted here.
    Deleting through the ORM would fire the before_delete listeners that queue
    *filesystem* removals (db/events.py) -- and a metadata mirror has no
    business touching files. Rows retracted upstream therefore linger until
    cleaned up by hand; for a nightly-regenerated public catalog that trade is
    fine. The artifact's schema stamp must match this code's, so a stale
    server refuses new artifacts instead of mangling them.
    """
    if source.startswith(("http://", "https://")):
        import httpx

        with tempfile.NamedTemporaryFile(suffix=".sqlite", delete=False) as tmp:
            with httpx.stream("GET", source, follow_redirects=True, timeout=60) as resp:
                resp.raise_for_status()
                for chunk in resp.iter_bytes():
                    tmp.write(chunk)
            local_path = tmp.name
    else:
        local_path = source

    src_engine = create_engine(f"sqlite:///{local_path}")
    summary: dict[str, int] = {}
    try:
        with src_engine.connect() as src:
            stamp = src.execute(select(AlembicVersion.__table__.c.version_num)).scalar()
            if stamp != TARGET_ALEMBIC_VERSION:
                raise ValueError(
                    f"Publish catalog schema {stamp!r} does not match this "
                    f"server's schema {TARGET_ALEMBIC_VERSION!r}; refusing to import"
                )
            with engine.begin() as dest:
                for model in _ORDERED_MODELS:
                    table = model.__table__
                    name = model.__tablename__
                    key_cols = _NATURAL_KEYS.get(
                        name, tuple(c.name for c in table.primary_key.columns)
                    )
                    preserve = _PRESERVE.get(name, ())
                    drop_id = name in _NATURAL_KEYS

                    rows = _rows(src, table)
                    for row in rows:
                        key = {k: row[k] for k in key_cols}
                        existing = dest.execute(select(table).filter_by(**key)).first()
                        if existing:
                            update = {
                                k: v
                                for k, v in row.items()
                                if k not in key_cols
                                and k not in preserve
                                and not (drop_id and k == "id")
                            }
                            if update:
                                dest.execute(table.update().filter_by(**key).values(**update))
                        else:
                            insert_row = dict(row)
                            if drop_id:
                                insert_row.pop("id", None)
                            dest.execute(table.insert().values(**insert_row))
                    summary[name] = len(rows)
    finally:
        src_engine.dispose()
        if local_path != source:
            Path(local_path).unlink(missing_ok=True)

    logger.info(f"Imported publish catalog from {source}: {summary}")
    return summary
