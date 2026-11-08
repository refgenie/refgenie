"""Constants and builders for engines, catalogs, assets and Refgenie worlds."""

import hashlib
import shutil
from pathlib import Path

TESTS_DATA_DIR = Path(__file__).parents[1] / "data"


RCRSD_FASTA = TESTS_DATA_DIR / "rCRSd.fa"


DEMO_FASTA = TESTS_DATA_DIR / "demo.fa"


# The standard built world: one rCRSd genome with one fasta asset named "test".
GENOME = "rCRSd"


GROUP = "fasta"


ASSET = "test"


#: The in-process ASGI server URLs used by ``serve_refgenie``.
MOCK_SERVER_URL = "http://mock-refgenie-server"


MOCK_MIRROR_URL = "http://mock-refgenie-mirror"


#: Sentinel: leave the key out of a payload / the attribute unset on a mock.
#: An *unset* MagicMock attribute returns a truthy MagicMock, which is
#: behaviourally different from ``[]`` -- some tests depend on that difference.
OMIT = object()


def fake_digest(label: str):
    """A well-formed ``GenomeDigest`` made deterministically from ``label``.

    The GA4GH sha512t24u of the label: 32 URL-safe base64 characters. For
    tests that need a genome digest but have no FASTA to take one from.
    """
    import base64

    from refgenie.models import GenomeDigest

    raw = hashlib.sha512(label.encode()).digest()[:24]
    return GenomeDigest(base64.urlsafe_b64encode(raw).decode("ascii"))


def make_engine():
    """A fresh in-memory SQLite engine (StaticPool, thread-check off).

    Plain function so module-level constants and multi-engine tests
    (server + client pairs, dual-world tests) can build engines directly;
    single-engine tests should use the ``engine`` fixture instead.
    """
    from sqlmodel import create_engine
    from sqlmodel.pool import StaticPool

    return create_engine(
        "sqlite://",
        connect_args={"check_same_thread": False},
        poolclass=StaticPool,
        echo=False,
    )


def register_fasta(r, fixtures_path=None):
    """Register fasta asset class and recipe from test data files.

    Idempotent -- safe to call multiple times on the same database.
    """
    from refgenie.exceptions import AssetClassExistsError, RecipeExistsError

    fixtures_path = Path(fixtures_path) if fixtures_path is not None else TESTS_DATA_DIR
    try:
        r.asset_class.add(fixtures_path / "fasta_asset_class.yaml")
    except AssetClassExistsError:
        pass
    try:
        r.recipe.add(fixtures_path / "fasta_asset_recipe.yaml")
    except RecipeExistsError:
        pass


def samtools_available() -> bool:
    """Check if samtools is available on PATH.

    No samtools? Run the tests inside the bulker crate to get it:
        bulker exec bulker_manifest.yaml -- python3 -m pytest -m e2e
    """
    return shutil.which("samtools") is not None


def build_rcrsd(r, asset_name=ASSET):
    """Build the standard fasta asset for the rCRSd genome."""
    return r.build.run(
        recipe_name="fasta",
        genome_alias=GENOME,
        asset_group_name=GROUP,
        asset_name=asset_name,
    )


def stage_rcrsd(r, asset_name=ASSET):
    """Stage the built rCRSd fasta asset (default serving_modes); return it."""
    genome_digest = r.alias.resolve(GENOME)
    built_asset = r.asset.get(
        genome_digest=genome_digest,
        asset_group_name=GROUP,
        asset_name=asset_name,
    )
    r.stage.create(
        asset=built_asset,
        genome_folder=r.genome_folder,
        genome_stage_folder=r.genome_stage_folder,
    )
    return built_asset


def make_built_refgenie(root, fixtures_path=None, *, asset_name=ASSET, build=True, stage=False):
    """A Refgenie on a real filesystem under ``root`` with the rCRSd genome.

    Factory form of the refgenie_fs/refgenie_built/staged_refgenie fixture
    chain, for tests that need two independent worlds or a named root.
    """
    from refgenie import Refgenie

    fixtures_path = Path(fixtures_path) if fixtures_path is not None else TESTS_DATA_DIR
    root = Path(root)
    r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    r.database.init(
        genome_folder=root / "genomes",
        genome_stage_folder=root / "archives",
    )
    register_fasta(r, fixtures_path)
    r.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=[GENOME],
        description="rCRSd genome",
    )
    if build:
        build_rcrsd(r, asset_name=asset_name)
    if stage:
        stage_rcrsd(r, asset_name=asset_name)
    return r


def make_server_rgc(tmp_path, fixtures_path=None, *, genomes=()):
    """A Refgenie over a fresh in-memory catalog with the fasta definitions.

    ``genomes`` is an iterable of ``(digest, alias_names)`` to register.
    """
    from refgenie import Refgenie

    rgc = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rgc.database.init(genome_folder=Path(tmp_path) / "genomes")
    register_fasta(rgc, fixtures_path)
    for digest, alias_names in genomes:
        rgc.genome.add(digest=digest, description="Test genome", alias_names=list(alias_names))
    return rgc


def seed_search_genomes(rgc):
    """Add two genomes (digests "a"*32 and "b"*32) whose searchable fields are all
    distinct, with aliases hg38 and mm39, and return the Refgenie.

    Aliases go through the mode-selected manager, not the SQL Alias table:
    ``refgenie_minimal`` is local (store-backed) mode.
    """
    from sqlmodel import Session

    from refgenie.db.tables import Genome

    with Session(rgc.database_engine) as session:
        session.add(
            Genome(
                digest="a" * 32,
                description="A Peculiar Description",
                species_name="Homo sapiens",
                taxon_id=9606,
                assembly_source="UCSC",
                assembly_accession="GCA_000001405.15",
            )
        )
        session.add(
            Genome(
                digest="b" * 32,
                description="unrelated",
                species_name="Mus musculus",
                taxon_id=10090,
                assembly_source="Ensembl",
                assembly_accession="GCA_000001635.9",
            )
        )
        session.commit()
    rgc.alias.add("hg38", "a" * 32)
    rgc.alias.add("mm39", "b" * 32)
    return rgc


def seed_species(rgc, taxon_id, scientific_name, common_name=None, synonyms=()):
    """Write one species-list row, as ``store sync`` or a pull would."""
    from refgenie.utils.species_file import SpeciesEntry

    rgc.species.replace(
        [
            SpeciesEntry(
                taxon_id=taxon_id,
                scientific_name=scientific_name,
                common_name=common_name,
                synonyms=list(synonyms),
            )
        ]
    )


def create_asset_class(session, name, serving_modes, version="0.1.0"):
    """Create an AssetClass with the given serving_modes."""
    from refgenie.db.tables import AssetClass

    ac = AssetClass(
        name=name,
        version=version,
        description=f"Test {name} asset class",
        serving_modes=list(serving_modes),
    )
    session.add(ac)
    session.commit()
    session.refresh(ac)
    return ac


def create_asset_on_disk(genome_folder, genome_digest, group_name, asset_name, files) -> Path:
    """Create asset files on disk and return the asset directory path.

    ``files`` is a ``{name: str | bytes}`` mapping.
    """
    asset_path = Path(genome_folder) / "data" / genome_digest / group_name / asset_name
    asset_path.mkdir(parents=True, exist_ok=True)
    for fname, content in files.items():
        if isinstance(content, bytes):
            (asset_path / fname).write_bytes(content)
        else:
            (asset_path / fname).write_text(content)
    return asset_path


def seed_asset(
    session,
    asset_group_id,
    asset_path,
    *,
    name="default",
    description="Test asset",
    size=0,
    is_default=True,
    seed_asset_name=True,
) -> str:
    """Insert a content-addressed Asset (plus its AssetName row); return the digest.

    The content digest is not computed by an insert-time handler, so callers
    that build a directory and seed the DB directly must compute it here.
    ``seed_asset_name=False`` creates an Asset with *no* AssetName.
    """
    from refgenie.db.tables import Asset, AssetName
    from refgenie.utils.digest import get_dir_digest

    digest = get_dir_digest(asset_path)
    session.add(
        Asset(
            digest=digest,
            name=name,
            description=description,
            asset_group_id=asset_group_id,
            size=size,
            path=str(asset_path),
        )
    )
    if seed_asset_name:
        session.add(
            AssetName(
                name=name,
                asset_group_id=asset_group_id,
                asset_digest=digest,
                is_default=is_default,
            )
        )
    session.commit()
    return digest


def seed_staged_asset(
    rgc,
    *,
    genome_digest,
    asset_group_name,
    asset_class,
    files,
    asset_name="default",
    stage_modes=(),
    tarball_bytes=None,
    description="Test asset",
    group_description=None,
    size=0,
    is_default=True,
    seed_asset_name=True,
) -> str:
    """Create an asset group + on-disk asset + DB rows; return the asset digest.

    ``asset_class`` is either an existing class name (looked up) or a
    ``(name, serving_modes[, version])`` tuple to synthesize one.
    ``stage_modes`` is a subset of ``{"file", "archive"}``. ``"file"`` creates
    the stage symlink plus a row with ``directory_contents``; ``"archive"``
    writes ``tarball_bytes`` if given (and records their real sha-256), or
    records a placeholder ``tarball_digest``/``tarball_size``.
    """
    from sqlmodel import Session

    from refgenie.db.tables import AssetGroup, StagedAsset

    engine = rgc.database_engine
    genome_folder = Path(str(rgc.genome_folder))
    stage_folder = Path(str(rgc.genome_stage_folder))
    contents = sorted(files)

    asset_path = create_asset_on_disk(
        genome_folder, genome_digest, asset_group_name, asset_name, files
    )

    with Session(engine) as session:
        if isinstance(asset_class, str):
            asset_class_id = rgc.asset_class.get(asset_class).id
        else:
            ac_name, serving_modes, *rest = asset_class
            asset_class_id = create_asset_class(
                session, ac_name, serving_modes, version=rest[0] if rest else "0.1.0"
            ).id

        group = AssetGroup(
            name=asset_group_name,
            description=group_description or f"Test {asset_group_name} asset group",
            genome_digest=genome_digest,
            asset_class_id=asset_class_id,
        )
        session.add(group)
        session.commit()
        session.refresh(group)

        digest = seed_asset(
            session,
            group.id,
            asset_path,
            name=asset_name,
            description=description,
            size=size,
            is_default=is_default,
            seed_asset_name=seed_asset_name,
        )

        if "file" in stage_modes:
            link_dir = stage_folder / genome_digest / asset_group_name / asset_name
            link_dir.parent.mkdir(parents=True, exist_ok=True)
            link_dir.symlink_to(asset_path)
            session.add(
                StagedAsset(
                    asset_digest=digest,
                    mode="file",
                    directory_contents=contents,
                    tarball_digest=None,
                    tarball_size=None,
                )
            )
            session.commit()

        if "archive" in stage_modes:
            if tarball_bytes is not None:
                archive_dir = stage_folder / genome_digest / asset_group_name
                archive_dir.mkdir(parents=True, exist_ok=True)
                (archive_dir / f"{digest}.tgz").write_bytes(tarball_bytes)
                tarball_digest = hashlib.sha256(tarball_bytes).hexdigest()
                tarball_size = len(tarball_bytes)
            else:
                tarball_digest = "ab" * 32
                tarball_size = 2048
            session.add(
                StagedAsset(
                    asset_digest=digest,
                    mode="archive",
                    directory_contents=contents,
                    tarball_digest=tarball_digest,
                    tarball_size=tarball_size,
                )
            )
            session.commit()

    return digest


def add_asset_from_files(
    r,
    *,
    rel_dir,
    files,
    asset_class_name,
    asset_group_name,
    asset_name,
    custom_seek_keys=None,
    genome=GENOME,
    contents="payload\n",
):
    """Materialize ``files`` under ``genome_folder/rel_dir`` and register them.

    ``{genome}`` in a filename is expanded to the genome digest, matching how
    assets are named on disk. Nested paths are created with ``parents=True``.
    """
    digest = r.alias.resolve(genome)
    asset_dir = r.genome_folder / rel_dir
    asset_dir.mkdir(parents=True, exist_ok=True)
    for name in files:
        path = asset_dir / name.replace("{genome}", digest)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents)
    return r.asset.content.add(
        genome_digest=digest,
        asset_group_name=asset_group_name,
        asset_name=asset_name,
        path=Path(rel_dir),
        asset_class_name=asset_class_name,
        custom_seek_keys=custom_seek_keys,
    )


def stage_copy(r, asset, dirname, mutate=b"") -> tuple[Path, Path]:
    """Copy an asset's content dir to a staging dir under the genome folder.

    Appending ``mutate`` bytes to the .fa makes the content digest differ.
    Returns ``(absolute_staging_path, path_relative_to_genome_folder)``; the
    relative form is what ``content.add`` expects.
    """
    staging = r.genome_folder / dirname
    shutil.copytree(r.genome_folder / asset.path, staging, symlinks=True)
    if mutate:
        fa = next(staging.glob("*.fa"))
        fa.write_bytes(fa.read_bytes() + mutate)
    return staging, Path(dirname)


def make_command_values(**overrides):
    """A minimal BuildCommandValues, with every field overridable."""
    from refgenie.models import BuildCommandValues

    defaults = dict(
        genome_digest="a" * 32,
        asset_group_name="test_group",
        genome_folder=Path("/tmp/genomes"),
        params={},
        files={},
        assets=None,
        custom_seek_keys={},
    )
    defaults.update(overrides)
    return BuildCommandValues(**defaults)


def build_on_disk_store(path, fasta, alias, namespace="name"):
    """Build (or add to) an on-disk RefgetStore from a FASTA with one curated alias.

    Returns the collection digest.
    """
    from refget.store import RefgetStore

    store = RefgetStore.on_disk(str(path))
    meta, _ = store.add_sequence_collection_from_fasta(str(fasta))
    store.add_collection_alias(namespace, alias, meta.digest)
    return meta.digest
