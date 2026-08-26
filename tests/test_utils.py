"""
Tests for ``refgenie.utils.*``: I/O and YAML, checksums, encryption, tarballs,
progress columns, and the confirmation prompts.

The directory-digest spec lives in test_dir_digest.py; colocation symlinks and
staging live in test_asset_content.py.
"""

import os
import tarfile
from subprocess import CalledProcessError
from unittest.mock import MagicMock, patch

import httpx
import pytest

from refgenie import Refgenie
from refgenie.const import USER_AGENT
from refgenie.utils.tarball import (
    copy_asset_dir,
    get_external_symlinks,
    tar,
)
from refgenie.utils.build import checksum
from refgenie.utils.encryption import decrypt_credentials, encrypt_dict
from refgenie.utils.io import coerce_cli_kwargs, read_yaml
from refgenie.utils.console import _DownloadColumn, _TransferSpeedColumn
from tests.helpers import make_engine


# SHA-256 of the rCRSd.fa fixture, computed from its raw bytes.
RCRSD_FA_SHA256 = "970b30744da503041786968946129347605be17513e8ac621f8a72aaae06a86d"


class TestIO:
    """refgenie.utils.io: YAML reading and CLI kwarg coercion."""

    def test_read_yaml_local_path(self, fixtures_path):
        """read_yaml loads a local YAML file given either a Path or a str."""
        yaml_file = fixtures_path / "fasta_asset_class.yaml"
        for arg in (yaml_file, str(yaml_file)):
            content = read_yaml(arg)
            assert isinstance(content, dict)
            assert content  # non-empty

    def test_read_yaml_url_follows_redirect_with_user_agent(self, monkeypatch):
        """
        read_yaml over a URL must carry a real User-Agent and follow a 302
        redirect (Cloudflare-fronted hosts reject requests with no User-Agent).
        """
        captured = {}

        def handler(request: httpx.Request) -> httpx.Response:
            captured["user_agent"] = request.headers.get("user-agent")
            if request.url.path == "/redirect":
                return httpx.Response(302, headers={"Location": "https://example.org/final.yaml"})
            return httpx.Response(200, text="name: test\nversion: 1.0.0\n")

        transport = httpx.MockTransport(handler)
        real_client = httpx.Client

        def client_factory(*args, **kwargs):
            kwargs["transport"] = transport
            return real_client(*args, **kwargs)

        monkeypatch.setattr(httpx, "Client", client_factory)

        content = read_yaml("https://example.org/redirect")
        assert content == {"name": "test", "version": "1.0.0"}
        assert captured["user_agent"] == USER_AGENT

    def test_coerce_cli_kwargs(self):
        """
        coerce_cli_kwargs converts bool/int strings to real types and leaves
        opaque values (memory strings) untouched; empty dict is a no-op.
        """
        result = coerce_cli_kwargs({"multi": "True", "cores": "4", "mem": "4000M"})
        assert result == {"multi": True, "cores": 4, "mem": "4000M"}
        assert coerce_cli_kwargs({}) == {}


class TestChecksum:
    """refgenie.utils.build.checksum: SHA-256 of file content."""

    def test_checksum_matches_known_sha256(self, fixtures_path):
        """checksum returns the SHA-256 hexdigest of the file's bytes."""
        assert checksum(fixtures_path / "rCRSd.fa") == RCRSD_FA_SHA256


class TestEncryption:
    """refgenie.utils.encryption: credential encrypt/decrypt round-trip."""

    def test_round_trip(self):
        creds = {"username": "test_user", "password": "test_pass"}
        encrypted = encrypt_dict(creds)
        assert isinstance(encrypted, str)
        assert encrypted != str(creds)  # not stored in plaintext
        assert decrypt_credentials(encrypted) == creds

    def test_none_and_empty(self):
        assert encrypt_dict(None) is None
        assert decrypt_credentials(None) is None
        assert encrypt_dict({}) == ""
        assert decrypt_credentials("") == {}


class TestArchiving:
    """refgenie.utils.tarball: tarballs, copies, and external symlinks."""

    def test_tar_raises_on_nonzero_exit(self, tmp_path):
        """tar() must fail loudly when the shell pipeline exits non-zero.

        Regression: run(cmd, shell=True) ignored the return code, so a failed
        or partial tarball still produced a StagedAsset pointing at it.
        """
        src = tmp_path / "asset"
        src.mkdir()
        (src / "f.txt").write_text("hello")
        # The output directory does not exist, so the shell redirection fails.
        output = tmp_path / "missing_dir" / "out.tgz"
        with pytest.raises(CalledProcessError):
            tar(src, output)

    def test_copy_asset_dir_raises_on_nonzero_exit(self, tmp_path):
        """copy_asset_dir() must fail loudly when rsync exits non-zero."""
        src = tmp_path / "src"
        src.mkdir()
        (src / "f.txt").write_text("hello")
        dest_parent = tmp_path / "readonly"
        dest_parent.mkdir()
        dest_parent.chmod(0o500)  # unwritable -> rsync exits non-zero
        try:
            with pytest.raises(CalledProcessError):
                copy_asset_dir(src, dest_parent / "target")
        finally:
            dest_parent.chmod(0o700)

    def test_get_external_symlinks(self, tmp_path):
        """Symlinks pointing outside the asset dir are detected; files are not."""
        asset_dir = tmp_path / "bwa_index" / "default"
        asset_dir.mkdir(parents=True)
        (asset_dir / "index.bwt").write_text("mock data")

        external_target = tmp_path / "fasta" / "default" / "genome.fa"
        external_target.parent.mkdir(parents=True)
        external_target.write_text(">chr1\nACGT\n")
        (asset_dir / "genome.fa").symlink_to(os.path.relpath(external_target, asset_dir))

        external = get_external_symlinks(asset_dir)
        assert "genome.fa" in external
        assert "index.bwt" not in external

    @pytest.mark.parametrize("exclude", [True, False])
    def test_tar_honors_exclude_files(self, tmp_path, exclude):
        """tar drops named members (here a real external symlink) only when asked."""
        asset_dir = tmp_path / "bwa_index" / "default"
        asset_dir.mkdir(parents=True)
        (asset_dir / "index.bwt").write_text("mock data")

        # A colocation-style external symlink is the thing we exclude.
        external_target = tmp_path / "parent" / "genome.fa"
        external_target.parent.mkdir(parents=True)
        external_target.write_text(">chr1\nACGT\n")
        (asset_dir / "genome.fa").symlink_to(os.path.relpath(external_target, asset_dir))

        tarball = tmp_path / "test.tgz"
        tar(asset_dir, tarball, exclude_files=["genome.fa"] if exclude else None)

        with tarfile.open(tarball, "r:gz") as tf:
            names = [m.name for m in tf.getmembers()]
        assert any("index.bwt" in n for n in names)
        assert (not any("genome.fa" in n for n in names)) is exclude


class TestProgressColumns:
    """refgenie.utils.console: rich progress column rendering edge cases."""

    def test_download_column_with_total(self):
        """Shows completed/total when total is known."""
        task = MagicMock(completed=1024 * 1024, total=10 * 1024 * 1024)
        result = str(_DownloadColumn.render(task))
        assert "1.0" in result and "10.0" in result and "MB" in result

    def test_download_column_with_none_total(self):
        """Handles None total without crashing."""
        task = MagicMock(completed=5 * 1024 * 1024, total=None)
        result = str(_DownloadColumn.render(task))
        assert "5.0" in result and "MB" in result and "/" not in result

    def test_transfer_speed_column_with_none(self):
        """Shows ? when speed is None."""
        task = MagicMock(speed=None)
        assert "?" in str(_TransferSpeedColumn.render(task))


class TestConfirmBulkPull:
    """The bulk-pull confirmation prompt utility."""

    def test_force_skips_prompt(self):
        from refgenie.utils.io import confirm_bulk_pull

        assert (
            confirm_bulk_pull(asset_count=100, genome_count=5, total_bytes=10 * 1024**3, force=True)
            is True
        )

    def test_prompt_shown_without_force(self):
        from refgenie.utils.io import confirm_bulk_pull

        seen = []
        result = confirm_bulk_pull(
            asset_count=10,
            genome_count=2,
            total_bytes=5 * 1024**3,
            force=False,
            confirm=lambda msg: seen.append(msg) or True,
        )
        assert result is True
        assert len(seen) == 1
        assert "10 asset(s)" in seen[0]
        assert "2 genome(s)" in seen[0]
        assert "5.0 GB" in seen[0]

    def test_user_declines_prompt(self):
        from refgenie.utils.io import confirm_bulk_pull

        assert (
            confirm_bulk_pull(
                asset_count=10,
                genome_count=2,
                total_bytes=5 * 1024**3,
                force=False,
                confirm=lambda msg: False,
            )
            is False
        )

    def test_no_confirmer_refuses_rather_than_reading_stdin(self):
        """A library caller with no terminal must not be prompted."""
        from refgenie.utils.io import confirm_bulk_pull

        with patch("rich.prompt.Confirm.ask", side_effect=AssertionError("read stdin")):
            assert (
                confirm_bulk_pull(
                    asset_count=10, genome_count=2, total_bytes=5 * 1024**3, force=False
                )
                is False
            )

    def test_unknown_size_message(self):
        from refgenie.utils.io import confirm_bulk_pull

        seen = []
        confirm_bulk_pull(
            asset_count=10,
            genome_count=2,
            total_bytes=0,
            force=False,
            confirm=lambda msg: seen.append(msg) or True,
        )
        assert "unknown size" in seen[0]


def _refuse_stdin(*args, **kwargs):
    raise AssertionError("library code read stdin")


class TestConfirmationPrompts:
    """Confirmation prompts never block on stdin and honour injected confirmers.

    Mostly unit tier; the library-entry-point test builds a real genome folder
    and carries its own ``component`` mark.
    """

    def test_large_archive_prompt_reports_the_cutoff_in_gb(self):
        """The confirmation message must state the cutoff in GB, not raw bytes
        with a GB suffix (``10.0GB``, never ``10000000000.0GB``)."""
        from refgenie.utils.build import should_pull_large_archive

        seen = []
        should_pull_large_archive(
            archive_size=50 * 1000**3,
            asset_registry_path="g/fasta:default",
            size_cutoff=10,
            confirm=lambda message: seen.append(message) or False,
        )

        assert len(seen) == 1
        assert "10.0GB" in seen[0], seen[0]
        assert "10000000000" not in seen[0]

    def test_large_archive_confirmation_is_injectable(self, monkeypatch):
        """`should_pull_large_archive` refuses an oversized archive by default
        (never prompting) and honours an injected confirmer."""
        import rich.prompt

        from refgenie.utils.build import should_pull_large_archive

        kwargs = dict(archive_size=50 * 1000**3, asset_registry_path="g/fasta:d", size_cutoff=10)

        monkeypatch.setattr(rich.prompt.Confirm, "ask", _refuse_stdin)
        assert should_pull_large_archive(**kwargs) is False
        assert should_pull_large_archive(**kwargs, confirm=lambda msg: True) is True

    def test_ask_declines_on_eof_instead_of_raising(self, monkeypatch):
        """`ask` is the interactive branch used once the CLI enables prompts.
        A closed stdin (piped input, CI, a subprocess) must degrade to a clean
        refusal, not a traceback."""
        import rich.prompt

        from refgenie.utils.prompt import ask

        monkeypatch.setattr(
            rich.prompt.Confirm, "ask", MagicMock(side_effect=EOFError("EOF when reading a line"))
        )

        assert ask("Would you like to subscribe to the default server?") is False

    def test_default_server_is_the_v4_host(self):
        """Guard against regressing to the legacy v3 refgenomes.databio.org,
        which 404s on /service-info for refgenie1."""
        from refgenie.const import DEFAULT_SERVER_URL

        assert DEFAULT_SERVER_URL == "https://api.refgenie.org"

    @pytest.mark.component
    @pytest.mark.parametrize(
        "call",
        [
            pytest.param(lambda r: r.purge(), id="purge"),
            pytest.param(
                lambda r: r.pull(asset_group_name="fasta", alias_name="x"), id="pull"
            ),
        ],
    )
    def test_library_never_blocks_on_stdin(self, tmp_path, monkeypatch, call):
        """A caller with no terminal must get a refusal, not a prompt."""
        import rich.prompt

        r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
        r.init(genome_folder=tmp_path / "genomes")

        monkeypatch.setattr(rich.prompt.Confirm, "ask", _refuse_stdin)
        assert call(r) is None
