"""Accessors for the standard built world and catalog, session and disk snapshots."""

from contextlib import contextmanager
from pathlib import Path
from unittest.mock import patch

from tests.support.factories import ASSET, GENOME, GROUP


def fasta_asset(r, name=ASSET, group=GROUP, genome=GENOME):
    """The built fasta asset ``name`` for ``genome``."""
    return r.asset.get(
        genome_digest=r.alias.resolve(genome),
        asset_group_name=group,
        asset_name=name,
    )


def genome_digest(r, genome=GENOME) -> str:
    return r.alias.resolve(genome)


def only_asset(r):
    """The single asset in a one-asset world."""
    return list(r.asset.list_assets())[0]


def content_dir(r) -> Path:
    """The digest-addressed content directory of the one built asset."""
    return r.genome_folder / only_asset(r).path


def build_flag(r, group=GROUP, asset=ASSET, genome=GENOME) -> Path:
    """The completion flag at the path advertised to snakemake."""
    template = r.build.target_template(group, asset)
    return Path(str(template).replace("{genome_name}", genome))


def build_dir(r, asset_name=ASSET, group=GROUP, genome=GENOME) -> Path:
    from refgenie.utils.paths import get_build_dir

    return get_build_dir(
        genome_folder=r.genome_folder,
        genome_name=genome,
        asset_group_name=group,
        asset_name=asset_name,
    )


def boom(*_args, **_kwargs):
    """Monkeypatch target that kills a write path mid-flight."""
    raise RuntimeError("killed before the commit")


def assets_rows(r):
    """Every Asset ORM row."""
    from sqlmodel import Session, select

    from refgenie.db.tables import Asset

    with Session(r.database_engine) as session:
        return session.exec(select(Asset)).unique().all()


def asset_name_rows(r):
    """Every AssetName ORM row."""
    from sqlmodel import Session, select

    from refgenie.db.tables import AssetName

    with Session(r.database_engine) as session:
        return session.exec(select(AssetName)).all()


def db_snapshot(r) -> dict:
    """Every table removal is supposed to touch, as plain comparable values."""
    from sqlmodel import Session, select

    from refgenie.db.tables import AssetGroup, Genome

    with Session(r.database_engine) as session:
        groups = sorted(
            (g.genome_digest, g.name) for g in session.exec(select(AssetGroup)).unique().all()
        )
        genomes = sorted(g.digest for g in session.exec(select(Genome)).unique().all())
    return {
        "genomes": genomes,
        "groups": groups,
        "assets": sorted(a.digest for a in assets_rows(r)),
        "asset_names": sorted(n.name for n in asset_name_rows(r)),
    }


@contextmanager
def session_depth():
    """Record how many manager sessions are open at once; yields a dict whose
    ``"max"`` is the peak. A peak above 1 means ``_database_session`` was
    re-entered from inside an open block (see ``ResourceManager``)."""
    import refgenie.managers.base as base

    real = base.database_session
    depth = {"now": 0, "max": 0}

    @contextmanager
    def counting(engine):
        depth["now"] += 1
        depth["max"] = max(depth["max"], depth["now"])
        try:
            with real(engine) as session:
                yield session
        finally:
            depth["now"] -= 1

    with patch.object(base, "database_session", counting):
        yield depth


def disk_snapshot(r) -> list[str]:
    """Every path under the genome folder, relative and sorted."""
    root = r.genome_folder
    return sorted(p.relative_to(root).as_posix() for p in root.rglob("*"))
