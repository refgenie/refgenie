"""Local file access ``/v1/files/*`` (``refgenie/server/local/files.py``).

Listing, download and preview of an asset's files on disk; the same-origin
policy that keeps the localhost bridge (the public site) away from file bytes;
and the path-safety rules (absolute paths, ``..``, symlink escapes).

Tier: the pure helpers are ``unit``; everything over a real built asset is
``component``.
"""

import gzip
import os

import pytest

from tests.helpers import (
    EVIL_ORIGIN,
    PUBLIC_ORIGIN,
    make_local_client,
    make_server_app,
    requires_dash,
    requires_server,
    route_keys,
    stub_rgc,
)

requires_server()
requires_dash()

DEV_ORIGIN = "http://localhost:5173"


class TestPureHelpers:
    pytestmark = pytest.mark.unit

    @pytest.mark.parametrize(
        "bad",
        [
            "",
            "/etc/passwd",
            "../secret",
            "a/../../b",
            "a/./b",
            "a//b",
            "a/",
            "..",
            "a\\..\\b",
            "C:/Windows",
            "a\x00b",
        ],
    )
    def test_validate_relative_path_rejects(self, bad):
        from fastapi import HTTPException

        from refgenie.server.local.files import validate_relative_path

        with pytest.raises(HTTPException) as info:
            validate_relative_path(bad)
        assert info.value.status_code == 400
        assert info.value.detail["code"] == "invalid_path"

    @pytest.mark.parametrize("good", ["x.fa", "index/x.1.bt2", "a.b/c-d_e.txt"])
    def test_validate_relative_path_accepts(self, good):
        from refgenie.server.local.files import validate_relative_path

        assert str(validate_relative_path(good)) == good

    def test_decode_text(self):
        from refgenie.server.local.files import decode_text

        assert decode_text(b"chr1\t100\n", truncated=False) == "chr1\t100\n"
        assert decode_text(b"a\x00b", truncated=False) is None
        assert decode_text(gzip.compress(b"hello"), truncated=False) is None
        assert decode_text(b"\xff\xfe\xfd", truncated=False) is None
        # A multi-byte character cut at the preview limit is dropped, not fatal.
        assert decode_text("abc\u00e9".encode()[:-1], truncated=True) == "abc"

    def test_capability_is_local_only(self):
        from refgenie.server.const import APP_MODE_LOCAL, APP_MODE_SERVER
        from refgenie.server.main import _capabilities

        assert _capabilities(APP_MODE_LOCAL)["local_files"] is True
        assert _capabilities(APP_MODE_SERVER)["local_files"] is False

    def test_server_mode_has_no_files_routes(self):
        app = make_server_app(stub_rgc())
        assert [p for _, p in route_keys(app) if p.startswith("/v1/files")] == []

    def test_unknown_asset_is_404(self):
        from unittest.mock import patch

        rgc = stub_rgc()
        with (
            patch("refgenie.server.local.files.get_asset_with_group", return_value=None),
            make_local_client(rgc) as client,
        ):
            response = client.get("/v1/files/assets/nope")
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "asset_not_found"


class TestLocalFiles:
    """Real Refgenie, real built asset, real filesystem."""

    pytestmark = pytest.mark.component

    @pytest.fixture()
    def world(self, refgenie_built, tmp_path):
        rg = refgenie_built
        genome_digest = rg.alias.resolve("rCRSd")
        asset = rg.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name="test"
        )
        asset_dir = rg.genome_folder / asset.path
        # Extra content: a nested directory, a binary file, a gzipped file, a
        # big text file, a colocation-style link inside the genome's content
        # store, and two escaping links.
        (asset_dir / "sub").mkdir()
        (asset_dir / "sub" / "notes.txt").write_text("hello\n")
        (asset_dir / "blob.bin").write_bytes(b"\x00\x01\x02binary")
        (asset_dir / "data.txt.gz").write_bytes(gzip.compress(b"chr1\t1\n"))
        (asset_dir / "big.txt").write_text("x" * (200 * 1024))
        sibling = rg.genome_folder / "data" / genome_digest / "sibling.txt"
        sibling.write_text("from a parent asset\n")
        os.symlink(os.path.relpath(sibling, asset_dir), asset_dir / "coloc.txt")
        outside = tmp_path / "outside-secret.txt"
        outside.write_text("do not leak\n")
        os.symlink(outside, asset_dir / "escape.txt")
        os.symlink(tmp_path, asset_dir / "escape_dir")
        with make_local_client(rg) as client:
            yield client, asset.digest, asset_dir

    # --- listing ----------------------------------------------------------

    def test_listing(self, world):
        client, digest, asset_dir = world
        response = client.get(f"/v1/files/assets/{digest}")
        assert response.status_code == 200
        assert response.headers["cross-origin-resource-policy"] == "same-origin"
        payload = response.json()
        assert payload["asset_digest"] == digest
        assert payload["truncated"] is False
        by_path = {e["path"]: e for e in payload["entries"]}
        assert by_path["sub"]["type"] == "dir"
        assert by_path["sub/notes.txt"] == {
            "path": "sub/notes.txt",
            "type": "file",
            "size": 6,
            "previewable": True,
        }
        assert by_path["blob.bin"]["previewable"] is False
        assert by_path["data.txt.gz"]["previewable"] is False
        assert by_path["coloc.txt"]["previewable"] is True
        # Links that leave the asset's content are not listed at all.
        assert "escape.txt" not in by_path
        assert "escape_dir" not in by_path
        fasta = [p for p in by_path if p.endswith(".fa")]
        assert fasta and by_path[fasta[0]]["previewable"] is True

    # --- download -----------------------------------------------------------

    def test_download_is_an_attachment(self, world):
        client, digest, _ = world
        response = client.get(f"/v1/files/assets/{digest}/raw/sub/notes.txt")
        assert response.status_code == 200
        assert response.content == b"hello\n"
        assert response.headers["content-disposition"].startswith("attachment")
        assert 'filename="notes.txt"' in response.headers["content-disposition"]
        assert response.headers["content-type"] == "application/octet-stream"
        assert response.headers["x-content-type-options"] == "nosniff"
        assert response.headers["cross-origin-resource-policy"] == "same-origin"

    def test_binary_and_gzip_download(self, world):
        client, digest, _ = world
        for name in ("blob.bin", "data.txt.gz"):
            response = client.get(f"/v1/files/assets/{digest}/raw/{name}")
            assert response.status_code == 200, name

    def test_colocation_link_inside_genome_store_is_served(self, world):
        client, digest, _ = world
        response = client.get(f"/v1/files/assets/{digest}/raw/coloc.txt")
        assert response.status_code == 200
        assert response.content == b"from a parent asset\n"

    # --- preview ------------------------------------------------------------

    def test_preview_small_text(self, world):
        client, digest, _ = world
        response = client.get(f"/v1/files/assets/{digest}/preview/sub/notes.txt")
        assert response.status_code == 200
        assert response.json() == {
            "path": "sub/notes.txt",
            "size": 6,
            "content": "hello\n",
            "truncated": False,
            "limit": 64 * 1024,
        }

    def test_preview_large_text_is_truncated(self, world):
        client, digest, _ = world
        payload = client.get(f"/v1/files/assets/{digest}/preview/big.txt").json()
        assert payload["truncated"] is True
        assert payload["size"] == 200 * 1024
        assert len(payload["content"]) == 64 * 1024

    @pytest.mark.parametrize("name", ["blob.bin", "data.txt.gz"])
    def test_preview_refuses_binary_and_gzip(self, world, name):
        client, digest, _ = world
        response = client.get(f"/v1/files/assets/{digest}/preview/{name}")
        assert response.status_code == 409
        assert response.json()["error"]["code"] == "not_previewable"

    # --- path safety --------------------------------------------------------

    @pytest.mark.parametrize("kind", ["raw", "preview"])
    @pytest.mark.parametrize(
        "path",
        ["escape.txt", "escape_dir/outside-secret.txt", "sub", "missing.txt"],
    )
    def test_escapes_and_missing_are_404(self, world, kind, path):
        client, digest, _ = world
        response = client.get(f"/v1/files/assets/{digest}/{kind}/{path}")
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "file_not_found"
        assert "do not leak" not in response.text

    @pytest.mark.parametrize(
        "path",
        [
            "..%2F..%2F..%2Fetc%2Fpasswd",
            "sub/..%2F..%2Fsecret",
            "%2Fetc%2Fpasswd",
            "sub%5C..%5Cnotes.txt",
        ],
    )
    def test_traversal_is_400(self, world, path):
        client, digest, _ = world
        response = client.get(f"/v1/files/assets/{digest}/raw/{path}")
        assert response.status_code == 400
        assert response.json()["error"]["code"] == "invalid_path"

    # --- origin policy --------------------------------------------------------

    @pytest.mark.parametrize("bridge_mode", ["read", "full"])
    @pytest.mark.parametrize("suffix", ["", "/raw/sub/notes.txt", "/preview/sub/notes.txt"])
    def test_bridge_origin_is_refused(self, world, bridge_mode, suffix):
        """The public site may never read file bytes, in any bridge mode."""
        client, digest, _ = world
        client.app.state.local_security.bridge_mode = bridge_mode
        response = client.get(
            f"/v1/files/assets/{digest}{suffix}", headers={"Origin": PUBLIC_ORIGIN}
        )
        assert response.status_code == 403
        assert response.json()["error"]["code"] == "forbidden_origin"
        assert "hello" not in response.text

    def test_bridge_regex_origin_is_refused(self, world):
        client, digest, _ = world
        settings = client.app.state.local_security
        settings.bridge_origin_regex = r"https://[a-z0-9-]+\.refgenie-ui\.pages\.dev"
        response = client.get(
            f"/v1/files/assets/{digest}/raw/sub/notes.txt",
            headers={"Origin": "https://preview-1.refgenie-ui.pages.dev"},
        )
        assert response.status_code == 403

    def test_unlisted_origin_is_refused(self, world):
        client, digest, _ = world
        response = client.get(
            f"/v1/files/assets/{digest}/raw/sub/notes.txt", headers={"Origin": EVIL_ORIGIN}
        )
        assert response.status_code == 403
        assert "access-control-allow-origin" not in response.headers

    @pytest.mark.parametrize("fetch_site", ["cross-site", "same-site"])
    def test_cross_site_embed_without_origin_is_refused(self, world, fetch_site):
        """``<img src>`` / ``<script src>`` on another page sends no Origin."""
        client, digest, _ = world
        response = client.get(
            f"/v1/files/assets/{digest}/raw/sub/notes.txt",
            headers={"Sec-Fetch-Site": fetch_site, "Sec-Fetch-Mode": "no-cors"},
        )
        assert response.status_code == 403

    @pytest.mark.parametrize(
        "headers",
        [
            {},
            {"Origin": "http://localhost"},
            {"Origin": DEV_ORIGIN},
            {"Sec-Fetch-Site": "same-origin"},
            {"Sec-Fetch-Site": "none"},
        ],
        ids=["no-origin", "same-origin", "dev-origin", "fetch-same-origin", "typed-url"],
    )
    def test_local_spa_is_allowed(self, world, headers):
        client, digest, _ = world
        response = client.get(f"/v1/files/assets/{digest}/raw/sub/notes.txt", headers=headers)
        assert response.status_code == 200
        assert response.content == b"hello\n"

    def test_foreign_host_is_refused(self, world):
        """DNS rebinding: the host guard still covers these routes."""
        client, digest, _ = world
        response = client.get(
            f"/v1/files/assets/{digest}/raw/sub/notes.txt", headers={"host": "evil.example"}
        )
        assert response.status_code == 421
