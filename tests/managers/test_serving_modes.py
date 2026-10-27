"""Tests for serving modes: validation and resolution on the AssetClass and
Asset models (unit tier), and StageManager.create()/remove() honoring them
against real genome folders and tarballs (component tier, marked per class).

Also here: DRS router serving-mode-aware responses -- the DRS endpoints return
structurally different responses depending on the asset's serving mode (file,
archive, both, none) and whether the request is asset-level or file-level
(component tier, marked per class).
"""

import hashlib
from pathlib import Path
from unittest.mock import MagicMock

import pytest
from pydantic import ValidationError
from sqlmodel import Session, select
from yaml import safe_load as yload

from refgenie import Refgenie
from refgenie.db.tables import Asset, AssetClass, AssetClassPublic, AssetPublic, StagedAsset
from tests.helpers import make_server_client, only_asset, requires_server, server_with_asset

requires_server()

from refgenie.server.routers.ga4gh_drs import (  # noqa: E402  (follows the guard)
    parse_refgenie_drs_object_id,
)
from refgenie.server.schemas import AssetResponse  # noqa: E402  (follows the guard)


# ---- Model validation tests (no database needed) ----


@pytest.mark.parametrize(
    "kwargs, expected",
    [
        ({}, ["file"]),  # default when unspecified
        ({"serving_modes": ["archive"]}, ["archive"]),
        ({"serving_modes": ["file", "archive"]}, ["file", "archive"]),
    ],
)
def test_asset_class_serving_modes_value(kwargs, expected):
    """serving_modes defaults to ['file'] and otherwise retains what was given."""
    ac = AssetClass(name="x", version="0.1", description=None, **kwargs)
    assert ac.serving_modes == expected


def test_asset_class_rejects_empty_serving_modes():
    """AssetClassPublic with empty serving_modes should raise ValidationError.

    Note: SQLModel table models (AssetClass) skip Pydantic validators.
    Validation is enforced via the Public (non-table) model.
    """
    with pytest.raises(ValidationError, match="serving_modes must not be empty"):
        AssetClassPublic(name="x", version="0.1", description=None, serving_modes=[])


def test_asset_class_rejects_invalid_serving_mode():
    """AssetClassPublic with invalid serving_mode value should raise ValidationError."""
    with pytest.raises(ValidationError, match="Invalid serving modes"):
        AssetClassPublic(name="x", version="0.1", description=None, serving_modes=["invalid"])


@pytest.mark.parametrize(
    "serving_modes, serves_files, serves_archive, is_metadata_only",
    [
        (["file"], True, False, False),
        (["archive"], False, True, False),
        (["none"], False, False, True),
        (["file", "archive"], True, True, False),
    ],
)
def test_asset_class_serving_mode_properties(
    serving_modes, serves_files, serves_archive, is_metadata_only
):
    """serves_files / serves_archive / is_metadata_only reflect the serving_modes."""
    ac = AssetClass(name="x", version="0.1", description=None, serving_modes=serving_modes)
    assert ac.serves_files is serves_files
    assert ac.serves_archive is serves_archive
    assert ac.is_metadata_only is is_metadata_only


# ---- Asset override tests (model-level, no database needed) ----


def test_asset_override_rejects_empty():
    """AssetPublic with empty serving_modes_override should raise ValidationError.

    Note: SQLModel table models (Asset) skip Pydantic validators.
    Validation is enforced via the Public (non-table) model.
    """
    with pytest.raises(ValidationError, match="serving_modes_override must not be empty when set"):
        AssetPublic(
            name="test", description=None, digest="a" * 64, size=0, serving_modes_override=[]
        )


def test_asset_override_rejects_invalid():
    """AssetPublic with invalid serving_modes_override should raise ValidationError."""
    with pytest.raises(ValidationError, match="Invalid serving modes"):
        AssetPublic(
            name="test", description=None, digest="a" * 64, size=0, serving_modes_override=["bad"]
        )


# ---- Database round-trip tests ----


def test_serving_modes_persisted_in_db(refgenie_with_fasta):
    """Fasta asset class loaded from YAML should have serving_modes=['file', 'archive']."""
    fasta_ac = refgenie_with_fasta.asset_class.get("fasta")
    assert fasta_ac.serving_modes == ["file", "archive"]


def test_metadata_only_asset_class_from_yaml(refgenie_minimal, fixtures_path):
    """metadata_only_asset_class.yaml should have serving_modes=['none'] and is_metadata_only."""
    refgenie_minimal.asset_class.add(fixtures_path / "metadata_only_asset_class.yaml")
    meta_ac = refgenie_minimal.asset_class.get("metadata_only")
    assert meta_ac.serving_modes == ["none"]
    assert meta_ac.is_metadata_only is True


# ---- to_yaml round-trip test ----


def test_to_yaml_includes_serving_modes(refgenie_with_fasta):
    """to_yaml() output should include serving_modes with correct value."""
    fasta_ac = refgenie_with_fasta.asset_class.get("fasta")
    yaml_str = fasta_ac.to_yaml()
    parsed = yload(yaml_str)
    assert "serving_modes" in parsed
    assert parsed["serving_modes"] == ["file", "archive"]


# ---- AssetResponse serialization tests ----
# AssetResponse is the response model for the assets endpoints; it must surface the
# *resolved* serving_modes (read from the Asset ORM property via from_attributes) so
# the puller learns the real serving mode.


def _mock_asset_group(serving_modes, class_name="fasta"):
    mock_asset_class = MagicMock()
    mock_asset_class.serving_modes = serving_modes
    mock_asset_class.name = class_name
    mock_asset_group = MagicMock()
    mock_asset_group.asset_class = mock_asset_class
    mock_asset_group.name = "default"
    mock_asset_group.genome_digest = "d" * 64
    return mock_asset_group


@pytest.mark.parametrize(
    "override, group_serving_modes, expected_modes, expected_class_name",
    [
        # override wins over the class default
        (["archive"], ["file"], ["archive"], "fasta"),
        # no override -> resolved from the asset class
        (None, ["file", "archive"], ["file", "archive"], "fasta"),
        # no asset_group -> ['archive'] default, no class name
        (None, None, ["archive"], None),
    ],
)
def test_asset_response_serializes_serving_modes(
    override, group_serving_modes, expected_modes, expected_class_name
):
    """AssetResponse surfaces the *resolved* serving_modes and asset_class_name."""
    a = Asset(name="test", description=None, digest="a" * 64, serving_modes_override=override)
    if group_serving_modes is not None:
        a.asset_group = _mock_asset_group(serving_modes=group_serving_modes)

    resp = AssetResponse.model_validate(a, from_attributes=True)
    assert resp.serving_modes == expected_modes
    assert resp.asset_class_name == expected_class_name


# ---- StageManager.create()/remove() respecting serving modes (component) ----


def _set_serving_modes(r: Refgenie, modes: list[str]):
    """Set serving_modes on the fasta AssetClass."""
    with Session(r.database_engine) as session:
        ac = session.exec(select(AssetClass).where(AssetClass.name == "fasta")).first()
        ac.serving_modes = modes
        session.add(ac)
        session.commit()


def _get_staged(r: Refgenie, asset_digest: str) -> list[StagedAsset]:
    """Get all StagedAsset records for an asset digest."""
    with Session(r.database_engine) as session:
        return list(
            session.exec(select(StagedAsset).where(StagedAsset.asset_digest == asset_digest)).all()
        )


def _stage_paths(r: Refgenie, asset: Asset) -> tuple[Path, Path]:
    """(tarball_path, file_mode_link_dir) for an asset in the stage folder.

    The tarball is content-addressed by digest; the file-mode dir of per-file
    symlinks is still name-based.
    """
    base = Path(r.genome_stage_folder) / asset.asset_group.genome_digest / asset.asset_group.name
    return base / f"{asset.digest}.tgz", base / asset.name


class TestStageCreateServingModes:
    """stage.create() produces exactly the StagedAsset records and on-disk
    artifacts that the asset class's serving_modes call for."""

    pytestmark = pytest.mark.component

    @pytest.mark.parametrize(
        "serving_modes, expected_modes, expect_tarball, expect_link_dir",
        [
            (["archive"], {"archive"}, True, False),
            (["file"], {"file"}, False, True),
            (["archive", "file"], {"archive", "file"}, True, True),
            (["none"], set(), False, False),
            ([], {"archive"}, True, False),  # empty -> default is archive
        ],
        ids=["archive", "file", "both", "none", "default"],
    )
    def test_stage_create(
        self, refgenie_built, serving_modes, expected_modes, expect_tarball, expect_link_dir
    ):
        r = refgenie_built
        _set_serving_modes(r, serving_modes)
        asset = only_asset(r)

        result = r.stage.create(asset, r.genome_folder, r.genome_stage_folder)

        staged = _get_staged(r, asset.digest)
        assert {sa.mode for sa in staged} == expected_modes
        assert len(staged) == len(expected_modes)
        if not expected_modes:
            assert result == []

        tarball_path, link_path = _stage_paths(r, asset)

        if "archive" in expected_modes:
            sa = next(sa for sa in staged if sa.mode == "archive")
            assert sa.tarball_digest is not None
            assert sa.tarball_size is not None
            assert sa.tarball_size > 0
        if "file" in expected_modes:
            sa = next(sa for sa in staged if sa.mode == "file")
            assert sa.tarball_digest is None
            assert sa.tarball_size is None
            assert sa.directory_contents is not None
            assert len(sa.directory_contents) > 0
            # Each entry is a symlink resolving to the matching file in the
            # (digest-addressed) data directory.
            expected_data_dir = Path(r.genome_folder) / asset.path
            children = list(link_path.iterdir())
            assert children
            for child in children:
                assert child.is_symlink()
                assert child.resolve() == (expected_data_dir / child.name).resolve()

        assert tarball_path.is_file() is expect_tarball
        assert link_path.is_dir() is expect_link_dir


class TestStageArchiveAtomicity:
    """_stage_archive() writes to a .partial sibling and swaps it onto the
    final path only once it is fully written and checksummed, so the final
    (StagedAsset-gated) path never holds a partial tarball."""

    pytestmark = pytest.mark.component

    def test_checksum_failure_leaves_no_tarball_or_partial(self, refgenie_built, monkeypatch):
        r = refgenie_built
        _set_serving_modes(r, ["archive"])
        asset = only_asset(r)

        def _boom(_path):
            raise RuntimeError("checksum blew up")

        monkeypatch.setattr("refgenie.managers.stage.checksum", _boom)

        with pytest.raises(RuntimeError, match="checksum blew up"):
            r.stage.create(asset, r.genome_folder, r.genome_stage_folder)

        tarball_path, _ = _stage_paths(r, asset)
        partial_path = tarball_path.with_name(tarball_path.name + ".partial")
        assert not tarball_path.exists()
        assert not partial_path.exists()
        assert _get_staged(r, asset.digest) == []

    def test_success_leaves_no_partial_file(self, refgenie_built):
        r = refgenie_built
        _set_serving_modes(r, ["archive"])
        asset = only_asset(r)

        r.stage.create(asset, r.genome_folder, r.genome_stage_folder)

        tarball_path, _ = _stage_paths(r, asset)
        partial_path = tarball_path.with_name(tarball_path.name + ".partial")
        assert tarball_path.is_file()
        assert not partial_path.exists()


class TestUnstageServingModes:
    """stage.remove() cleans up every staged artifact and record while leaving
    the original (digest-addressed) data untouched."""

    pytestmark = pytest.mark.component

    @pytest.mark.parametrize(
        "serving_modes",
        [["archive"], ["file"], ["archive", "file"]],
        ids=["archive", "file", "both"],
    )
    def test_unstage(self, refgenie_built, serving_modes):
        r = refgenie_built
        _set_serving_modes(r, serving_modes)
        asset = only_asset(r)
        r.stage.create(asset, r.genome_folder, r.genome_stage_folder)

        tarball_path, link_path = _stage_paths(r, asset)
        if "archive" in serving_modes:
            assert tarball_path.is_file()
        if "file" in serving_modes:
            assert link_path.is_dir()

        # Remember the original data path (content is digest-addressed).
        original_data_path = Path(r.genome_folder) / asset.path
        assert original_data_path.exists()

        r.stage.remove(asset_digest=asset.digest)

        assert not tarball_path.exists()
        assert not link_path.exists()
        assert _get_staged(r, asset.digest) == []
        # Original files untouched
        assert original_data_path.exists()


# ---------------------------------------------------------------------------
# DRS router serving-mode-aware responses
#
# Component tier (marked per class): these build real genome folders, asset
# files and archives on disk alongside SQLite. Deselected from the bare
# `pytest` inner loop; run with `pytest -m component`.
# ---------------------------------------------------------------------------

# Placeholder tarball bytes written to disk for archive-mode staged assets.
# StagedAsset.tarball_digest is the sha-256 of THESE bytes -- the true
# "downloaded bytes" digest, distinct from the asset identity digest used for
# the object id/alias/S3 key. StageManager computes it when it writes the
# tarball; these fixtures seed it directly, as the real writer would.
TARBALL_BYTES = b"\x1f\x8b" + b"\x00" * 20
TARBALL_DIGEST = hashlib.sha256(TARBALL_BYTES).hexdigest()

GENOME_DIGEST = "d" * 32


# --- Parsing tests ---


@pytest.mark.component
class TestParseDrsObjectId:
    """Tests for parse_refgenie_drs_object_id()."""

    def test_asset_level(self):
        """digest-only ID returns (digest, None)."""
        digest = "a" * 64
        result = parse_refgenie_drs_object_id(digest)
        assert result == (digest, None)

    def test_file_level(self):
        """{digest}:{filename} returns (digest, filename)."""
        digest = "b" * 64
        result = parse_refgenie_drs_object_id(f"{digest}:test.fa")
        assert result == (digest, "test.fa")

    @pytest.mark.parametrize("object_id", ["abc", ""])
    def test_malformed_id(self, object_id):
        """Short or empty ID raises HTTPException (400)."""
        from fastapi import HTTPException

        with pytest.raises(HTTPException) as exc_info:
            parse_refgenie_drs_object_id(object_id)
        assert exc_info.value.status_code == 400

    def test_colon_in_filename(self):
        """{digest}:path/to:file returns (digest, 'path/to:file')."""
        digest = "c" * 64
        result = parse_refgenie_drs_object_id(f"{digest}:path/to:file")
        assert result == (digest, "path/to:file")


# --- The four serving modes, as one seeded server ---

_ASSET_SPECS = [
    {
        "key": "archive",
        "asset_group_name": "archive_asset",
        "asset_class": ("archive_type", ["archive"]),
        "files": {"archive_file.dat": "archive data content"},
        "stage_modes": ("archive",),
        "tarball_bytes": TARBALL_BYTES,
        "description": "Archive-only test asset",
    },
    {
        "key": "file",
        "asset_group_name": "file_asset",
        "asset_class": ("file_type", ["file"]),
        "files": {"test.fa": ">chr1\nACGT\n", "test.fa.fai": "chr1\t4\t6\t4\t5\n"},
        "stage_modes": ("file",),
        "description": "File-only test asset",
    },
    {
        "key": "combo",
        "asset_group_name": "combo_asset",
        "asset_class": ("combo_type", ["file", "archive"]),
        "files": {"data.bin": b"binary data"},
        "stage_modes": ("file", "archive"),
        "tarball_bytes": TARBALL_BYTES,
        "description": "Combo test asset",
    },
    {
        # Metadata-only: an Asset row with no AssetName and nothing staged.
        "key": "none",
        "asset_group_name": "none_asset",
        "asset_class": ("none_type", ["none"]),
        "files": {"metadata.txt": "just metadata"},
        "seed_asset_name": False,
        "description": "Metadata-only test asset",
    },
]


@pytest.fixture
def drs_client_with_assets(tmp_path):
    """TestClient with assets in every serving mode, for DRS testing."""
    rgc, digests = server_with_asset(
        tmp_path,
        genome_digest=GENOME_DIGEST,
        alias="testgenome",
        assets=_ASSET_SPECS,
    )
    with make_server_client(rgc) as client:
        client.digests = digests
        yield client


# --- Metadata tests ---


@pytest.mark.component
class TestDrsMetadataArchiveMode:
    """DRS metadata for archive-mode assets."""

    def test_archive_mode_local_only(self, drs_client_with_assets):
        """Archive asset returns simple object with access_id 'archive', no contents."""
        digest = drs_client_with_assets.digests["archive"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}")
        assert response.status_code == 200
        data = response.json()

        assert data["id"] == digest
        assert data["contents"] is None  # Not a bundle
        assert len(data["access_methods"]) == 1
        assert data["access_methods"][0]["access_id"] == "archive"
        assert "archive" in data["access_methods"][0]["access_url"]["url"]
        # Checksum must be the downloaded-tarball byte digest, NOT the identity
        # digest (which is the object id/alias/S3 key).
        assert data["checksums"][0]["type"] == "sha-256"
        assert data["checksums"][0]["checksum"] == TARBALL_DIGEST
        assert data["checksums"][0]["checksum"] != digest


@pytest.mark.component
class TestDrsMetadataFileMode:
    """DRS metadata for file-mode assets."""

    def test_file_mode_local_only(self, drs_client_with_assets):
        """File asset returns bundle with contents listing files, access_id 'file'."""
        digest = drs_client_with_assets.digests["file"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}")
        assert response.status_code == 200
        data = response.json()

        assert data["id"] == digest
        # Should be a bundle (contents populated)
        assert data["contents"] is not None
        assert len(data["contents"]) == 2
        content_names = [c["name"] for c in data["contents"]]
        assert "test.fa" in content_names
        assert "test.fa.fai" in content_names

        # Each content item should have a synthetic DRS ID
        for item in data["contents"]:
            assert item["id"].startswith(digest + ":")
            assert item["drs_uri"] is not None

        # Access methods should include file
        assert len(data["access_methods"]) == 1
        assert data["access_methods"][0]["access_id"] == "file"


@pytest.mark.component
class TestDrsMetadataComboMode:
    """DRS metadata for combo (file + archive) assets."""

    def test_combo_mode(self, drs_client_with_assets):
        """Combo asset returns bundle with contents AND both access methods."""
        digest = drs_client_with_assets.digests["combo"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}")
        assert response.status_code == 200
        data = response.json()

        # Should be a bundle (has file mode)
        assert data["contents"] is not None
        assert len(data["contents"]) == 1
        assert data["contents"][0]["name"] == "data.bin"

        # Should have both access methods
        access_ids = {m["access_id"] for m in data["access_methods"]}
        assert "archive" in access_ids
        assert "file" in access_ids


@pytest.mark.component
class TestDrsMetadataNoneMode:
    """DRS metadata for none-mode (metadata-only) assets."""

    def test_none_mode(self, drs_client_with_assets):
        """None asset returns simple object, empty access_methods, null contents."""
        digest = drs_client_with_assets.digests["none"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}")
        assert response.status_code == 200
        data = response.json()

        assert data["contents"] is None
        assert data["access_methods"] == []
        # None mode serves no bytes and has no tarball -> no checksum advertised
        # (advertising the identity digest here would be a false byte checksum).
        assert data["checksums"] == []


@pytest.mark.component
class TestDrsMetadataFileLevel:
    """DRS metadata for file-level requests ({digest}:{filename})."""

    def test_file_level_object(self, drs_client_with_assets):
        """File-level request returns simple object with file-mode access methods only."""
        digest = drs_client_with_assets.digests["file"]
        object_id = f"{digest}:test.fa"
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{object_id}")
        assert response.status_code == 200
        data = response.json()

        assert data["id"] == object_id
        assert data["name"] == "test.fa"
        assert data["contents"] is None  # Simple object, not a bundle
        assert data["description"] == "File: test.fa"

        # Only file access methods
        assert len(data["access_methods"]) == 1
        assert data["access_methods"][0]["access_id"] == "file"

    def test_file_level_not_in_contents(self, drs_client_with_assets):
        """File not in StagedAsset.directory_contents returns 404."""
        digest = drs_client_with_assets.digests["file"]
        object_id = f"{digest}:nonexistent.txt"
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{object_id}")
        assert response.status_code == 404

    def test_file_level_wrong_mode(self, drs_client_with_assets):
        """File-level request on archive-only asset returns 404."""
        digest = drs_client_with_assets.digests["archive"]
        object_id = f"{digest}:archive_file.dat"
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{object_id}")
        assert response.status_code == 404

    def test_file_level_combo_asset(self, drs_client_with_assets):
        """File-level request on combo asset works (combo supports file mode)."""
        digest = drs_client_with_assets.digests["combo"]
        object_id = f"{digest}:data.bin"
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{object_id}")
        assert response.status_code == 200
        data = response.json()
        assert data["name"] == "data.bin"
        assert data["contents"] is None


@pytest.mark.component
class TestDrsChecksumSourcing:
    """
    The DRS ``checksums`` field must describe the bytes a client downloads
    through an access method (the tarball's sha-256 = StagedAsset.tarball_digest),
    NOT the content-identity digest (asset.digest = object id / alias / S3 key).
    """

    def test_combo_checksum_is_tarball_digest(self, drs_client_with_assets):
        """Combo bundle advertises the archive tarball digest (its downloadable bytes)."""
        digest = drs_client_with_assets.digests["combo"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}")
        data = response.json()
        assert data["checksums"] == [{"checksum": TARBALL_DIGEST, "type": "sha-256"}]
        assert data["checksums"][0]["checksum"] != digest

    def test_file_only_checksum_omitted(self, drs_client_with_assets):
        """File-only asset has no tarball digest -> no checksum advertised."""
        digest = drs_client_with_assets.digests["file"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}")
        data = response.json()
        assert data["checksums"] == []

    def test_file_level_checksum_omitted(self, drs_client_with_assets):
        """File-level object (single file) has no per-file byte digest -> omitted."""
        digest = drs_client_with_assets.digests["file"]
        object_id = f"{digest}:test.fa"
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{object_id}")
        data = response.json()
        assert data["checksums"] == []


# --- Access URL tests ---


@pytest.mark.component
class TestDrsAccessUrl:
    """Tests for the DRS access URL endpoint."""

    def test_local_archive(self, drs_client_with_assets):
        """/access/archive returns URL ending in /bytes."""
        digest = drs_client_with_assets.digests["archive"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}/access/archive")
        assert response.status_code == 200
        data = response.json()
        assert data["url"].endswith("/bytes")
        assert "archive" in data["url"]

    def test_local_file(self, drs_client_with_assets):
        """/access/file for file-level ID returns URL ending in /bytes."""
        digest = drs_client_with_assets.digests["file"]
        object_id = f"{digest}:test.fa"
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{object_id}/access/file")
        assert response.status_code == 200
        data = response.json()
        assert data["url"].endswith("/bytes")

    def test_remote_no_remote_configured(self, drs_client_with_assets):
        """Remote access_id with no matching remote returns 404."""
        digest = drs_client_with_assets.digests["archive"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}/access/s3:1:archive")
        assert response.status_code == 404

    def test_malformed_access_id(self, drs_client_with_assets):
        """Malformed remote access_id returns 400."""
        digest = drs_client_with_assets.digests["archive"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}/access/s3:bad")
        assert response.status_code == 400

    def test_unknown_access_id(self, drs_client_with_assets):
        """Completely unknown access_id returns 404."""
        digest = drs_client_with_assets.digests["archive"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}/access/bogus")
        assert response.status_code == 404


# --- Remote redirects (pushed links on http/https remotes) ---

CDN = "https://cdn.example.com/assets"


@pytest.fixture
def pushed_combo_client(tmp_path):
    """The combo asset pushed, in both modes, to an https remote and an s3 remote."""
    from refgenie.db.tables import RemoteType

    rgc, digests = server_with_asset(
        tmp_path, genome_digest=GENOME_DIGEST, alias="testgenome", assets=_ASSET_SPECS
    )
    digest = digests["combo"]
    https = rgc.remote.add(name="cdn", type=RemoteType.https, prefix=CDN)
    s3 = rgc.remote.add(name="bucket", type=RemoteType.s3, prefix="s3://bucket/assets")
    for remote in (https, s3):
        for mode in ("archive", "file"):
            rgc.remote.link(remote.id, digest, mode, pushed=True)
    with make_server_client(rgc) as client:
        client.digest, client.https, client.s3 = digest, https, s3
        yield client


@pytest.mark.component
class TestRemoteRedirects:
    """DRS remote access methods and v4 download redirects come from ``rgc.remote``."""

    def test_drs_lists_one_method_per_pushed_https_link(self, pushed_combo_client):
        c = pushed_combo_client
        methods = c.get(f"/ga4gh/drs/objects/{c.digest}").json()["access_methods"]
        remote = {m["access_id"]: m["access_url"]["url"] for m in methods if ":" in m["access_id"]}
        assert set(remote) == {f"https:{c.https.id}:archive", f"https:{c.https.id}:file"}
        assert remote[f"https:{c.https.id}:archive"] == (
            f"{CDN}/{GENOME_DIGEST}/combo_asset/{c.digest}.tgz"
        )
        assert all(url.startswith(CDN + "/") for url in remote.values())

    def test_drs_access_url_for_remote(self, pushed_combo_client):
        c = pushed_combo_client
        base = f"/ga4gh/drs/objects/{c.digest}/access"
        resp = c.get(f"{base}/https:{c.https.id}:archive")
        assert resp.status_code == 200
        assert resp.json()["url"] == f"{CDN}/{GENOME_DIGEST}/combo_asset/{c.digest}.tgz"
        # An s3 remote has no https URL to hand out.
        assert c.get(f"{base}/s3:{c.s3.id}:archive").status_code == 404

    def test_v4_archive_download_redirects(self, pushed_combo_client):
        c = pushed_combo_client
        resp = c.get(f"/v4/archives/{c.digest}/download", follow_redirects=False)
        assert resp.status_code == 307
        assert resp.headers["location"] == f"{CDN}/{GENOME_DIGEST}/combo_asset/{c.digest}.tgz"

    def test_v4_file_download_redirects(self, pushed_combo_client):
        c = pushed_combo_client
        resp = c.get(f"/v4/assets/{c.digest}/files/data.bin", follow_redirects=False)
        assert resp.status_code == 307
        location = resp.headers["location"]
        assert location.startswith(f"{CDN}/{GENOME_DIGEST}/combo_asset/")
        assert location.endswith("/data.bin")


# --- Bytes tests ---


@pytest.mark.component
class TestDrsBytes:
    """Tests for the DRS bytes endpoint (local serving)."""

    def test_archive_bytes(self, drs_client_with_assets):
        """/access/archive/bytes returns tarball."""
        digest = drs_client_with_assets.digests["archive"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}/access/archive/bytes")
        assert response.status_code == 200
        assert response.headers.get("content-type") in (
            "application/gzip",
            "application/x-gzip",
        )

    def test_file_bytes_redirect(self, drs_client_with_assets):
        """/access/file/bytes redirects to v4 file endpoint."""
        digest = drs_client_with_assets.digests["file"]
        object_id = f"{digest}:test.fa"
        response = drs_client_with_assets.get(
            f"/ga4gh/drs/objects/{object_id}/access/file/bytes",
            follow_redirects=False,
        )
        assert response.status_code == 307
        assert f"/v4/assets/{digest}/files/test.fa" in response.headers["location"]

    def test_file_bytes_requires_filename(self, drs_client_with_assets):
        """/access/file/bytes without file_path in object_id returns 400."""
        digest = drs_client_with_assets.digests["file"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}/access/file/bytes")
        assert response.status_code == 400

    def test_rejects_remote_access_id(self, drs_client_with_assets):
        """/access/s3:1:archive/bytes returns 400."""
        digest = drs_client_with_assets.digests["archive"]
        response = drs_client_with_assets.get(
            f"/ga4gh/drs/objects/{digest}/access/s3:1:archive/bytes"
        )
        assert response.status_code == 400


# --- Publication gate tests (security review finding F1) ---
#
# A file on disk at the predictable stage path is never, by itself, evidence
# that an asset is published: the StagedAsset row is. These tests cover the
# gap the pre-fix code had -- a tarball with no row, a row whose mode is no
# longer in the asset's effective serving_modes, and a pushed remote link for
# an unpublished mode -- across every route that serves bytes or hands out an
# access URL (helpers.require_published).

PUBLICATION_CDN = "https://cdn.example.com/publication-gate"


@pytest.fixture
def drs_client_publication(tmp_path):
    """TestClient with the usual serving-mode assets, plus an orphan tarball
    (StagedAsset row deleted after staging+pushing, tarball and pushed link
    left behind) and a stale archive row (mode no longer in serving_modes)."""
    from sqlmodel import Session, select

    from refgenie.db.tables import RemoteType, StagedAsset

    extra_specs = [
        {
            "key": "orphan",
            "asset_group_name": "orphan_asset",
            "asset_class": ("orphan_type", ["archive"]),
            "files": {"orphan_file.dat": "orphan data content"},
            "stage_modes": ("archive",),
            "tarball_bytes": TARBALL_BYTES,
            "description": "Orphan tarball test asset",
        },
        {
            "key": "stale",
            "asset_group_name": "stale_asset",
            "asset_class": ("stale_type", ["file"]),
            "files": {"stale_file.dat": "stale data content"},
            "stage_modes": ("archive",),
            "tarball_bytes": TARBALL_BYTES,
            "description": "Stale-mode test asset (archive row, file-only class)",
        },
    ]
    rgc, digests = server_with_asset(
        tmp_path,
        genome_digest=GENOME_DIGEST,
        alias="testgenome",
        assets=_ASSET_SPECS + extra_specs,
    )

    orphan_digest = digests["orphan"]
    remote = rgc.remote.add(name="cdn", type=RemoteType.https, prefix=PUBLICATION_CDN)
    # Link while the StagedAsset row still exists (link() requires it), then
    # delete the row -- leaving the tarball and the "pushed" link behind. This
    # is the exact F1 scenario: publication evidence gone, artifacts remain.
    rgc.remote.link(remote.id, orphan_digest, "archive", pushed=True)
    # A pushed file-mode link, so file-level remote URLs can be checked.
    rgc.remote.link(remote.id, digests["file"], "file", pushed=True)
    with Session(rgc.database_engine) as session:
        sa = session.exec(
            select(StagedAsset).where(
                StagedAsset.asset_digest == orphan_digest, StagedAsset.mode == "archive"
            )
        ).one()
        session.delete(sa)
        session.commit()

    with make_server_client(rgc) as client:
        client.digests = digests
        client.remote = remote
        yield client


@pytest.mark.component
class TestDrsPublicationGate:
    """Every DRS byte-serving/access-URL route, and the v4 siblings they
    redirect or delegate to, refuse to serve a mode with no StagedAsset row."""

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_orphan_archive_access_url_404(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["orphan"]
        assert c.get(f"{prefix}/objects/{digest}/access/archive").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_orphan_archive_bytes_404(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["orphan"]
        assert c.get(f"{prefix}/objects/{digest}/access/archive/bytes").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_orphan_remote_access_url_404(self, drs_client_publication, prefix):
        """A pushed remote link for an unpublished mode is not honored either."""
        c = drs_client_publication
        digest = c.digests["orphan"]
        url = f"{prefix}/objects/{digest}/access/https:{c.remote.id}:archive"
        assert c.get(url).status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_orphan_metadata_lists_no_access_methods(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["orphan"]
        data = c.get(f"{prefix}/objects/{digest}").json()
        assert data["access_methods"] == []

    def test_orphan_v4_download_404(self, drs_client_publication):
        c = drs_client_publication
        digest = c.digests["orphan"]
        assert c.get(f"/v4/archives/{digest}/download").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_metadata_only_archive_url_404(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["none"]
        assert c.get(f"{prefix}/objects/{digest}/access/archive").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_metadata_only_file_url_404(self, drs_client_publication, prefix):
        """File-level id so the 400 (missing file_path) rule doesn't shadow
        the 404 (mode not served) rule."""
        c = drs_client_publication
        object_id = f"{c.digests['none']}:metadata.txt"
        assert c.get(f"{prefix}/objects/{object_id}/access/file").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_metadata_only_archive_bytes_404(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["none"]
        assert c.get(f"{prefix}/objects/{digest}/access/archive/bytes").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_metadata_only_file_bytes_404(self, drs_client_publication, prefix):
        c = drs_client_publication
        object_id = f"{c.digests['none']}:metadata.txt"
        assert c.get(f"{prefix}/objects/{object_id}/access/file/bytes").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_file_only_archive_url_404(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["file"]
        assert c.get(f"{prefix}/objects/{digest}/access/archive").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_file_only_archive_bytes_404(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["file"]
        assert c.get(f"{prefix}/objects/{digest}/access/archive/bytes").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_file_only_file_url_200(self, drs_client_publication, prefix):
        c = drs_client_publication
        object_id = f"{c.digests['file']}:test.fa"
        response = c.get(f"{prefix}/objects/{object_id}/access/file")
        assert response.status_code == 200
        assert response.json()["url"].endswith("/bytes")

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_file_only_missing_file_404(self, drs_client_publication, prefix):
        c = drs_client_publication
        object_id = f"{c.digests['file']}:missing.txt"
        assert c.get(f"{prefix}/objects/{object_id}/access/file").status_code == 404
        assert c.get(f"{prefix}/objects/{object_id}/access/file/bytes").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_file_only_no_file_path_400(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["file"]
        assert c.get(f"{prefix}/objects/{digest}/access/file").status_code == 400

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_valid_archive_url_200(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["archive"]
        response = c.get(f"{prefix}/objects/{digest}/access/archive")
        assert response.status_code == 200
        assert response.json()["url"].endswith("/bytes")

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_valid_archive_bytes_200(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["archive"]
        response = c.get(f"{prefix}/objects/{digest}/access/archive/bytes")
        assert response.status_code == 200
        assert response.content == TARBALL_BYTES

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_valid_archive_file_url_404(self, drs_client_publication, prefix):
        """File-level id so the 400 (missing file_path) rule doesn't shadow
        the 404 (mode not served) rule -- this asset serves archive only."""
        c = drs_client_publication
        object_id = f"{c.digests['archive']}:whatever"
        assert c.get(f"{prefix}/objects/{object_id}/access/file").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_stale_archive_bytes_404(self, drs_client_publication, prefix):
        """StagedAsset(archive) row exists, but the class serves file only."""
        c = drs_client_publication
        digest = c.digests["stale"]
        assert c.get(f"{prefix}/objects/{digest}/access/archive/bytes").status_code == 404

    def test_v4_archives_lists_only_served_archives(self, drs_client_publication):
        c = drs_client_publication
        listed = {r["asset_digest"] for r in c.get("/v4/archives").json()["items"]}
        assert c.digests["archive"] in listed
        assert c.digests["stale"] not in listed
        assert c.digests["orphan"] not in listed

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_metadata_only_file_object_metadata_404(self, drs_client_publication, prefix):
        c = drs_client_publication
        object_id = f"{c.digests['none']}:metadata.txt"
        assert c.get(f"{prefix}/objects/{object_id}").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_remote_file_access_url_points_at_file(self, drs_client_publication, prefix):
        c = drs_client_publication
        object_id = f"{c.digests['file']}:test.fa"
        resp = c.get(f"{prefix}/objects/{object_id}/access/https:{c.remote.id}:file")
        assert resp.status_code == 200
        assert resp.json()["url"].endswith("/test.fa")

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_file_object_remote_method_points_at_file(self, drs_client_publication, prefix):
        c = drs_client_publication
        object_id = f"{c.digests['file']}:test.fa"
        methods = c.get(f"{prefix}/objects/{object_id}").json()["access_methods"]
        remote = [m for m in methods if m["access_id"] == f"https:{c.remote.id}:file"]
        assert len(remote) == 1
        assert remote[0]["access_url"]["url"].endswith("/test.fa")

    def test_stale_v4_download_404(self, drs_client_publication):
        c = drs_client_publication
        digest = c.digests["stale"]
        assert c.get(f"/v4/archives/{digest}/download").status_code == 404

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_bad_remote_mode_400(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["archive"]
        assert c.get(f"{prefix}/objects/{digest}/access/https:1:bogus").status_code == 400

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_remote_mode_none_400(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["archive"]
        assert c.get(f"{prefix}/objects/{digest}/access/https:1:none").status_code == 400

    @pytest.mark.parametrize("prefix", ["/ga4gh/drs", "/v4/ga4gh/drs"])
    def test_remote_bad_id_400(self, drs_client_publication, prefix):
        c = drs_client_publication
        digest = c.digests["archive"]
        url = f"{prefix}/objects/{digest}/access/https:notanint:archive"
        assert c.get(url).status_code == 400


# --- DRS spec compliance ---


@pytest.mark.component
class TestDrsSpecCompliance:
    """Verify response shapes match GA4GH DRS 1.2.0 basics."""

    def test_self_uri_uses_drs_scheme(self, drs_client_with_assets):
        """self_uri uses drs:// scheme."""
        digest = drs_client_with_assets.digests["archive"]
        response = drs_client_with_assets.get(f"/ga4gh/drs/objects/{digest}")
        data = response.json()
        assert data["self_uri"].startswith("drs://")
