"""
The SQL half of the publish catalog: reading a build node's publishable rows,
writing them to a fresh SQLite file, and upserting such a file into a server.

These copy raw tables between databases. Normal code reads and writes remotes
and their links through ``RemoteManager`` (``rgc.remote``). What to publish,
and how rows are rewritten on the way, is ``refgenie.publish_catalog``.
"""

from pathlib import Path

from sqlalchemy import create_engine, func, select
from sqlalchemy.engine import Engine
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
    SeekKey,
    Species,
    SpeciesName,
    StagedAsset,
    Store,
)

# FK-safe order for both export inserts and import upserts.
ORDERED_MODELS = [
    Configuration,
    Remote,
    Store,
    Genome,
    Alias,
    Species,
    SpeciesName,
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

# Tables whose surrogate id nothing references: match rows by natural key and
# never import the artifact's id, which a rebuilt build node may have renumbered.
_NATURAL_KEYS = {
    "stagedasset": ("asset_digest", "mode"),
    "assetname": ("asset_group_id", "name"),
}
# Tables whose id other tables DO reference, matched by natural key too. The
# server keeps its own id, and artifact rows pointing at the table are
# translated to it (_REMAP_FKS).
_REMAPPED_KEYS = {
    "assetgroup": ("genome_digest", "name"),
}
# {table: {column: referenced table}} for columns holding a _REMAPPED_KEYS id.
_REMAP_FKS = {
    "asset": {"asset_group_id": "assetgroup"},
    "assetname": {"asset_group_id": "assetgroup"},
}
# Server-owned columns an import must never overwrite on an existing row.
_PRESERVE = {"stagedasset": ("download_count",), "configuration": ("plugin_settings",)}


#: Tables the import rewrites per taxon instead of upserting, so a species
#: name removed upstream leaves the server too.
_PER_TAXON = ("species", "speciesname")


def _rows(conn, table, where=None) -> list[dict]:
    """All rows of a table as plain dicts, optionally filtered."""
    stmt = select(table)
    if where is not None:
        stmt = stmt.where(where)
    return [dict(r._mapping) for r in conn.execute(stmt)]


def _replace_species(src, dest) -> dict[str, int]:
    """Rewrite, in ``dest``, every taxon the artifact ``src`` carries. Other taxa stay."""
    species, names = Species.__table__, SpeciesName.__table__
    species_rows, name_rows = _rows(src, species), _rows(src, names)
    ids = [r["taxon_id"] for r in species_rows]
    if ids:
        dest.execute(names.delete().where(names.c.taxon_id.in_(ids)))
        dest.execute(species.delete().where(species.c.taxon_id.in_(ids)))
        dest.execute(species.insert(), species_rows)
        if name_rows:
            dest.execute(names.insert(), name_rows)
    return {"species": len(species_rows), "speciesname": len(name_rows)}


def _disable_absent_stores(dest, store_rows: list[dict]) -> list[str]:
    """Set ``enabled=False`` on ``dest`` stores the artifact does not carry.

    An artifact with no store rows leaves the server's stores alone rather
    than disabling every one.

    Returns:
        The names of the stores disabled.
    """
    if not store_rows:
        return []
    store = Store.__table__
    kept = [r["name"] for r in store_rows]
    disabled = [
        r.name
        for r in dest.execute(
            select(store.c.name).where(store.c.enabled == True, store.c.name.not_in(kept))  # noqa: E712
        )
    ]
    if disabled:
        dest.execute(store.update().where(store.c.name.in_(disabled)).values(enabled=False))
    return disabled


def read_publishable(engine: Engine) -> tuple[set[tuple[str, str]], dict[str, list[dict]]]:
    """
    The rows a publish catalog is built from, as plain dicts.

    Assets, their names and seek keys are restricted to the published ones
    (at least one pushed link); asset groups to the groups of those assets.
    Every other table is read whole. The ``alias`` rows are the SQL half of
    the alias space only.

    Returns:
        The pushed ``(asset_digest, mode)`` pairs, and rows by table name.
    """
    t = {m.__tablename__: m.__table__ for m in ORDERED_MODELS}
    with engine.connect() as src:
        pushed_pairs = {
            (r.asset_digest, r.mode)
            for r in src.execute(
                select(t["remoteassetlink"].c.asset_digest, t["remoteassetlink"].c.mode).where(
                    t["remoteassetlink"].c.pushed == True  # noqa: E712
                )
            )
        }
        published = {digest for digest, _ in pushed_pairs}
        published_group_ids = {
            r.asset_group_id
            for r in src.execute(
                select(t["asset"].c.asset_group_id).where(t["asset"].c.digest.in_(published))
            )
        }
        rows = {
            name: _rows(src, t[name])
            for name in (
                "configuration",
                "store",
                "genome",
                "alias",
                "species",
                "speciesname",
                "assetclass",
                "assetclassseekkey",
                "assetclasslink",
                "recipe",
                "recipeassetclassesinputs",
                "assetlink",
                "stagedasset",
            )
        }
        rows["assetgroup"] = _rows(
            src, t["assetgroup"], t["assetgroup"].c.id.in_(published_group_ids)
        )
        for name in ("asset", "seekkey"):
            column = t[name].c.digest if name == "asset" else t[name].c.asset_digest
            rows[name] = _rows(src, t[name], column.in_(published))
        rows["assetname"] = _rows(src, t["assetname"], t["assetname"].c.asset_digest.in_(published))
    return pushed_pairs, rows


def write_catalog(dest_path: Path, table_rows: dict[str, list[dict]]) -> dict[str, int]:
    """
    Write ``table_rows`` to a fresh SQLite file at ``dest_path``, stamped at the head.

    Every table in :data:`ORDERED_MODELS` must have an entry, possibly empty.

    Returns:
        A table -> row-count summary.
    """
    dest_engine = create_engine(f"sqlite:///{dest_path}")
    SQLModel.metadata.create_all(dest_engine)
    summary: dict[str, int] = {}
    try:
        with dest_engine.begin() as dest:
            for model in ORDERED_MODELS:
                name = model.__tablename__
                rows = table_rows[name]
                if rows:
                    dest.execute(model.__table__.insert(), rows)
                summary[name] = len(rows)
            dest.execute(
                AlembicVersion.__table__.insert(), [{"version_num": TARGET_ALEMBIC_VERSION}]
            )
    finally:
        dest_engine.dispose()
    return summary


def _free_id(dest, table, wanted: int) -> int:
    """``wanted`` if no row of ``table`` has that id yet, else one past the highest.

    Not left to the database: rows written with explicit ids leave a Postgres id
    sequence behind the rows it numbers, so ``nextval`` can hand out a taken id.
    """
    if dest.execute(select(table.c.id).where(table.c.id == wanted)).first() is None:
        return wanted
    return (dest.execute(select(func.max(table.c.id))).scalar() or 0) + 1


def _check_stamp(src, action: str) -> None:
    """Raise ValueError unless the artifact ``src`` is stamped at this code's schema."""
    stamp = src.execute(select(AlembicVersion.__table__.c.version_num)).scalar()
    if stamp != TARGET_ALEMBIC_VERSION:
        raise ValueError(
            f"Publish catalog schema {stamp!r} does not match this "
            f"server's schema {TARGET_ALEMBIC_VERSION!r}; refusing to {action}"
        )


def import_catalog(engine: Engine, local_path: str) -> dict[str, int]:
    """
    Upsert the publish catalog SQLite file at ``local_path`` into ``engine``.

    Insert-or-update only, in one transaction; nothing is ever deleted, except
    that each species taxon in the artifact has its row and names rewritten.
    Stores are matched by name, and a server store the artifact does not carry
    is disabled (genomes keep a foreign key to it, so it is never deleted).

    Returns:
        A table -> row-count summary.

    Raises:
        ValueError: If the file's schema stamp is not this code's.
    """
    src_engine = create_engine(f"sqlite:///{local_path}")
    summary: dict[str, int] = {}
    try:
        with src_engine.connect() as src:
            _check_stamp(src, "import")
            with engine.begin() as dest:
                # {table: {artifact id: server id}} for _REMAPPED_KEYS tables.
                id_maps: dict[str, dict[int, int]] = {t: {} for t in _REMAPPED_KEYS}
                for model in ORDERED_MODELS:
                    table = model.__table__
                    name = model.__tablename__
                    if name in _PER_TAXON:
                        continue
                    natural = _NATURAL_KEYS.get(name) or _REMAPPED_KEYS.get(name)
                    key_cols = natural or tuple(c.name for c in table.primary_key.columns)
                    preserve = _PRESERVE.get(name, ())
                    drop_id = natural is not None
                    remapped = name in _REMAPPED_KEYS
                    fks = _REMAP_FKS.get(name, {})

                    rows = _rows(src, table)
                    for row in rows:
                        for col, ref in fks.items():
                            if row.get(col) is not None:
                                row[col] = id_maps[ref].get(row[col], row[col])
                        key = {k: row[k] for k in key_cols}
                        existing = dest.execute(select(table).filter_by(**key)).first()
                        if remapped:
                            id_maps[name][row["id"]] = (
                                existing.id if existing else _free_id(dest, table, row["id"])
                            )
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
                            if remapped:
                                insert_row["id"] = id_maps[name][row["id"]]
                            elif drop_id:
                                insert_row["id"] = _free_id(dest, table, row["id"])
                            dest.execute(table.insert().values(**insert_row))
                    summary[name] = len(rows)
                    if name == "store":
                        summary["store_disabled"] = len(_disable_absent_stores(dest, rows))
                summary.update(_replace_species(src, dest))
    finally:
        src_engine.dispose()
    return summary


def prune_catalog(engine: Engine, local_path: str, dry_run: bool = False) -> dict[str, list[str]]:
    """
    Delete, from ``engine``, the assets the publish catalog at ``local_path`` retracted.

    :func:`import_catalog` never deletes, so this removes assets of genomes the
    artifact carries whose digests it lacks, with their dependent rows and any
    asset group left empty. Uses Core deletes, so no files are removed. Run it
    after importing the same artifact: an older one would retract newer assets.

    Raises:
        ValueError: If the file's schema stamp is not this code's.

    Returns:
        ``{"asset": [...], "assetgroup": [...]}``: what was (or, with ``dry_run``,
        would be) deleted, as asset digests and ``genome_digest/group`` names.
    """
    asset, group = Asset.__table__, AssetGroup.__table__
    src_engine = create_engine(f"sqlite:///{local_path}")
    try:
        with src_engine.connect() as src:
            _check_stamp(src, "prune")
            genomes = {r.digest for r in src.execute(select(Genome.__table__.c.digest))}
            kept_assets = {r.digest for r in src.execute(select(asset.c.digest))}
            kept_groups = {
                (r.genome_digest, r.name)
                for r in src.execute(select(group.c.genome_digest, group.c.name))
            }
    finally:
        src_engine.dispose()

    with engine.begin() as dest:
        stale_assets = [
            r.digest
            for r in dest.execute(
                select(asset.c.digest)
                .join(group, asset.c.asset_group_id == group.c.id)
                .where(group.c.genome_digest.in_(genomes))
            )
            if r.digest not in kept_assets
        ]
        live_group_ids = {
            r.asset_group_id
            for r in dest.execute(
                select(asset.c.asset_group_id).where(asset.c.digest.not_in(stale_assets))
            )
        }
        stale_groups = [
            r
            for r in dest.execute(
                select(group.c.id, group.c.genome_digest, group.c.name).where(
                    group.c.genome_digest.in_(genomes)
                )
            )
            if (r.genome_digest, r.name) not in kept_groups and r.id not in live_group_ids
        ]
        if not dry_run and stale_assets:
            links, names = AssetLink.__table__, AssetName.__table__
            seeks, staged = SeekKey.__table__, StagedAsset.__table__
            remote_links = RemoteAssetLink.__table__
            dest.execute(remote_links.delete().where(remote_links.c.asset_digest.in_(stale_assets)))
            dest.execute(staged.delete().where(staged.c.asset_digest.in_(stale_assets)))
            dest.execute(seeks.delete().where(seeks.c.asset_digest.in_(stale_assets)))
            dest.execute(names.delete().where(names.c.asset_digest.in_(stale_assets)))
            dest.execute(
                links.delete().where(
                    links.c.parent_digest.in_(stale_assets) | links.c.child_digest.in_(stale_assets)
                )
            )
            dest.execute(asset.delete().where(asset.c.digest.in_(stale_assets)))
        if not dry_run and stale_groups:
            dest.execute(group.delete().where(group.c.id.in_([g.id for g in stale_groups])))
    return {
        "asset": stale_assets,
        "assetgroup": [f"{g.genome_digest}/{g.name}" for g in stale_groups],
    }
