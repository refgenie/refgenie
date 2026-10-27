"""Component tests for refgenie's push/publish path.

Covers the `refgenie remote` handlers (status, remove), the push workflow (per-asset handler, the
content-addressed skip-if-present path, folder-sync strategy, S3 helpers), and
the publish-catalog export/import path (publish_catalog.py), which builds a
real asset on tmp_path, exports a publish catalog, and imports it into a second
(server-like) database.

RemoteManager itself (``rgc.remote``: naming, links, ``pushed_urls``, push
intent from staging) is tested in tests/managers/test_remote.py.

The client-side unit half (store-mode selection, RefgenieserverClient
pagination, download progress) lives in tests/managers/test_server_client.py; remote
*genome-source* unit tests live in tests/managers/test_sources.py.

Tiers: the S3-helper tests are pure `unit` (no Refgenie instance, no
filesystem). Everything else builds real genome folders, tarballs and SQLite
together, so it is marked `component` per class -- there is deliberately no
module-level mark.
"""

import pathlib
import subprocess
from unittest.mock import MagicMock, patch

import pytest
from sqlalchemy import create_engine
from sqlmodel import SQLModel, select

from refgenie import Refgenie
from refgenie.publish_catalog import export_publish_catalog, import_publish_catalog
from refgenie.const import TARGET_ALEMBIC_VERSION
from refgenie.db.tables import (
    AlembicVersion,
    Alias,
    Asset,
    Genome,
    Remote,
    RemoteAssetLink,
    RemoteType,
    StagedAsset,
)
from tests.helpers import (
    fasta_asset,
    make_engine,
    make_server_client,
    register_fasta,
    requires_server,
    seed_staged_asset,
)

HTTPS_PREFIX = "https://example.org/assets"


# ===========================================================================
# Helpers and component fixtures (real genome folder + tarball + SQLite)
# ===========================================================================


@pytest.fixture
def s3_remote_env(staged_refgenie):
    """A staged refgenie plus an s3 remote (with a push_command) holding one
    unpushed archive link. Shared by the remote-status and push test groups."""
    r = staged_refgenie
    asset = fasta_asset(r)
    remote = r.remote.add(
        type=RemoteType.s3,
        prefix="my-bucket",
        name="test-s3",
        push_command="aws s3 cp {local_path} s3://{prefix}/{relative_path}",
    )
    r.remote.link(
        ref=remote.id,
        asset_digest=asset.digest,
        mode="archive",
        pushed=False,
    )
    return r, remote, asset


def _unpushed(r, ref=None):
    """The unpushed links, optionally for one remote."""
    return [link for link, _, _ in r.remote.unpushed(ref)]


def _make_status_cmd(remote=None):
    cmd = MagicMock()
    cmd.remote = remote
    return cmd


def _make_push_cmd(
    dry_run=False, remote=None, genome=None, strategy="per_asset", genome_digest=None
):
    """Create a mock PushModel-like command object."""
    cmd = MagicMock()
    cmd.dry_run = dry_run
    cmd.remote = remote
    cmd.genome = genome
    cmd.genome_digest = genome_digest
    cmd.strategy = strategy
    return cmd


def _add_s3_remote(r, asset):
    """Attach a real s3://-prefixed remote + unpushed archive link; return the remote."""
    s3_remote = r.remote.add(
        type=RemoteType.s3,
        prefix="s3://my-bucket/assets",
        name="s3-real",
        push_command="aws s3 cp {local_path} {prefix}/{relative_path}",
    )
    r.remote.link(
        ref=s3_remote.id,
        asset_digest=asset.digest,
        mode="archive",
        pushed=False,
    )
    return s3_remote


@pytest.fixture
def sync_env(staged_refgenie):
    """A staged refgenie with two remotes for folder-sync push."""
    r = staged_refgenie
    asset = fasta_asset(r)
    remotes = []
    for prefix, desc in (("bucket-a", "sync-s3-a"), ("bucket-b", "sync-s3-b")):
        remote = r.remote.add(
            type=RemoteType.s3,
            prefix=prefix,
            name=desc,
            push_command="aws s3 sync {genome_stage_folder} s3://{prefix}/ --follow-symlinks",
        )
        r.remote.link(
            ref=remote.id,
            asset_digest=asset.digest,
            mode="archive",
            pushed=False,
        )
        remotes.append(remote)
    return r, remotes[0], remotes[1], asset


class TestRemoteStatusHandler:
    """The handle_remote_status CLI handler."""

    pytestmark = pytest.mark.component

    def test_shows_counts(self, s3_remote_env, capsys):
        """Shows the correct pushed/unpushed counts."""
        from refgenie.cli.commands.remote import handle_remote_status

        r, remote, asset = s3_remote_env
        handle_remote_status(_make_status_cmd(), r)
        output = capsys.readouterr().out
        assert "test-s3" in output
        assert "Unpushed: 1" in output or "Unpushed:   1" in output

    def test_no_remotes_does_not_raise(self, tmp_path, capsys):
        """With no remotes configured the handler returns cleanly."""
        from refgenie.cli.commands.remote import handle_remote_status

        r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
        r.database.init(genome_folder=tmp_path / "genomes", genome_stage_folder=tmp_path / "stage")
        handle_remote_status(_make_status_cmd(), r)  # must not raise

    def test_filter_by_remote(self, s3_remote_env, capsys):
        """--remote shows only the named remote."""
        from refgenie.cli.commands.remote import handle_remote_status

        r, remote, asset = s3_remote_env
        r.remote.add(type=RemoteType.https, prefix="cdn.example.com", name="cdn-remote")
        capsys.readouterr()  # drop the "Added remote" log line
        handle_remote_status(_make_status_cmd(remote="test-s3"), r)
        output = capsys.readouterr().out
        assert "test-s3" in output
        assert "cdn-remote" not in output

    def test_lists_unpushed_digests(self, s3_remote_env, capsys):
        """Unpushed assets are listed by digest and mode."""
        from refgenie.cli.commands.remote import handle_remote_status

        r, remote, asset = s3_remote_env
        handle_remote_status(_make_status_cmd(), r)
        output = capsys.readouterr().out
        assert asset.digest in output
        assert "mode=archive" in output

    def test_unknown_remote_exits_not_found(self, s3_remote_env):
        from refgenie.cli.commands.remote import handle_remote_status
        from refgenie.cli.errors import EXIT_NOT_FOUND

        r, _, _ = s3_remote_env
        with pytest.raises(SystemExit) as excinfo:
            handle_remote_status(_make_status_cmd(remote="no-such-remote"), r)
        assert excinfo.value.code == EXIT_NOT_FOUND


# ===========================================================================
# Component: the push workflow (per-asset, skip-if-present, folder-sync)
# ===========================================================================


class TestHandlePush:
    """Tests for the per-asset handle_push handler."""

    pytestmark = pytest.mark.component

    def test_dry_run_pushes_nothing(self, s3_remote_env):
        """Dry run lists what would be pushed without running the command."""
        from refgenie.cli.commands.push import handle_push

        r, remote, asset = s3_remote_env
        with patch("subprocess.run") as mock_run:
            handle_push(_make_push_cmd(dry_run=True), r)
            mock_run.assert_not_called()

        links = _unpushed(r)
        assert len(links) == 1
        assert links[0].pushed is False

    def test_success_marks_pushed(self, s3_remote_env, fake_subprocess):
        """A successful push marks the link pushed=True."""
        from refgenie.cli.commands.push import handle_push

        r, remote, asset = s3_remote_env
        handle_push(_make_push_cmd(), r)

        assert len(_unpushed(r)) == 0

    def test_failed_upload_stays_unpushed_and_exits_nonzero(
        self, s3_remote_env, failing_subprocess
    ):
        """A failed upload leaves the link pushed=False for retry and, since
        every upload failed, the command exits non-zero."""
        from refgenie.cli.commands.push import handle_push
        from refgenie.cli.errors import EXIT_GENERAL_ERROR

        r, remote, asset = s3_remote_env
        failing_subprocess(subprocess.CalledProcessError(1, "aws s3 cp", stderr="err"))
        with pytest.raises(SystemExit) as excinfo:
            handle_push(_make_push_cmd(), r)

        assert excinfo.value.code == EXIT_GENERAL_ERROR
        links = _unpushed(r)
        assert len(links) == 1
        assert links[0].pushed is False

    def test_partial_failure_exits_nonzero_but_marks_successes(self, s3_remote_env):
        """One success + one failure: the success is marked pushed, the
        failure stays unpushed, and the command still exits non-zero."""
        from refgenie.cli.commands.push import handle_push

        r, remote, asset = s3_remote_env
        remote2 = r.remote.add(
            type=RemoteType.s3,
            prefix="other-bucket",
            name="other-s3",
            push_command="aws s3 cp {local_path} s3://{prefix}/{relative_path}",
        )
        r.remote.link(
            ref=remote2.id,
            asset_digest=asset.digest,
            mode="archive",
            pushed=False,
        )

        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = [
                MagicMock(stdout="", stderr="", returncode=0),
                subprocess.CalledProcessError(1, "aws s3 cp", stderr="err"),
            ]
            with pytest.raises(SystemExit) as excinfo:
                handle_push(_make_push_cmd(), r)

        assert excinfo.value.code != 0
        # Exactly one link (the failed one) remains unpushed.
        assert len(_unpushed(r)) == 1

    def test_no_push_command_skips(self, staged_refgenie):
        """A remote with no push_command is a skip, not a failure -- nothing
        is uploaded. But since nothing pushed and a skip is all that
        happened, the command still exits non-zero: "0 pushed, 1 skipped,
        exit 0" would be the worst outcome, silently telling the caller the
        push succeeded."""
        from refgenie.cli.commands.push import handle_push
        from refgenie.cli.errors import EXIT_GENERAL_ERROR

        r = staged_refgenie
        remote = r.remote.add(
            type=RemoteType.s3,
            prefix="my-bucket",
            name="no-push-cmd",
            push_command=None,
        )
        r.remote.link(
            ref=remote.id,
            asset_digest=fasta_asset(r).digest,
            mode="archive",
            pushed=False,
        )

        with patch("subprocess.run") as mock_run:
            with pytest.raises(SystemExit) as excinfo:
                handle_push(_make_push_cmd(), r)
            mock_run.assert_not_called()

        assert excinfo.value.code == EXIT_GENERAL_ERROR
        assert len(_unpushed(r)) == 1

    def test_nothing_to_push_is_noop(self, s3_remote_env):
        """When all links are already pushed, handle_push does nothing and does not raise."""
        from refgenie.cli.commands.push import handle_push

        r, remote, asset = s3_remote_env
        r.remote.mark_pushed(ref=remote.id, asset_digest=asset.digest, mode="archive")

        with patch("subprocess.run") as mock_run:
            handle_push(_make_push_cmd(), r)
            mock_run.assert_not_called()

    def test_filter_by_remote(self, s3_remote_env, fake_subprocess):
        """--remote pushes only to the named remote."""
        from refgenie.cli.commands.push import handle_push

        r, remote, asset = s3_remote_env
        remote2 = r.remote.add(
            type=RemoteType.s3,
            prefix="other-bucket",
            name="other-s3",
            push_command="aws s3 cp {local_path} s3://{prefix}/{relative_path}",
        )
        r.remote.link(
            ref=remote2.id,
            asset_digest=asset.digest,
            mode="archive",
            pushed=False,
        )

        handle_push(_make_push_cmd(remote="other-s3"), r)

        # First remote untouched; second remote pushed.
        first = _unpushed(r, remote.id)
        assert len(first) == 1
        assert first[0].pushed is False
        assert len(_unpushed(r, remote2.id)) == 0

    def test_unknown_remote_exits_not_found(self, s3_remote_env):
        """--remote naming no remote is an error, not "nothing to push"."""
        from refgenie.cli.commands.push import handle_push
        from refgenie.cli.errors import EXIT_NOT_FOUND

        r, _, _ = s3_remote_env
        with patch("subprocess.run") as mock_run, pytest.raises(SystemExit) as excinfo:
            handle_push(_make_push_cmd(remote="no-such-remote"), r)
        mock_run.assert_not_called()
        assert excinfo.value.code == EXIT_NOT_FOUND
        assert len(_unpushed(r)) == 1

    def test_filter_by_genome(self, s3_remote_env, fake_subprocess):
        """--genome pushes only that genome's assets; other genomes stay unpushed."""
        from refgenie.cli.commands.push import handle_push

        r, remote, asset = s3_remote_env
        # A second genome with its own staged asset linked to the same remote.
        # If --genome were ignored, this link would be pushed too.
        other_genome = "other_genome_digest_42"
        r.genome.add(digest=other_genome, description="other", alias_names=["other"])
        other_digest = seed_staged_asset(
            r,
            genome_digest=other_genome,
            asset_group_name="fasta",
            asset_class="fasta",
            files={"other.fa": ">chr1\nACGT\n"},
            stage_modes=("archive",),
        )
        r.remote.link(
            ref=remote.id,
            asset_digest=other_digest,
            mode="archive",
            pushed=False,
        )

        handle_push(_make_push_cmd(genome="rCRSd"), r)

        # Only the filtered genome's link was pushed; the other genome's remains.
        remaining = _unpushed(r)
        assert [link.asset_digest for link in remaining] == [other_digest]

    def test_command_template_substitution(self, s3_remote_env, fake_subprocess):
        """{local_path}, {relative_path}, and {prefix} are all substituted."""
        from refgenie.cli.commands.push import handle_push

        r, remote, asset = s3_remote_env
        handle_push(_make_push_cmd(), r)
        executed_cmd = fake_subprocess.call_args[0][0]

        assert "{local_path}" not in executed_cmd
        assert "{relative_path}" not in executed_cmd
        assert "{prefix}" not in executed_cmd
        assert "my-bucket" in executed_cmd
        assert executed_cmd.startswith("aws s3 cp")

    def test_relative_path_is_content_addressed(self, s3_remote_env, fake_subprocess):
        """The pushed key is the content digest ({digest}.tgz), not the asset name."""
        from refgenie.cli.commands.push import handle_push

        r, remote, asset = s3_remote_env
        handle_push(_make_push_cmd(), r)
        executed_cmd = fake_subprocess.call_args[0][0]

        assert f"{asset.digest}.tgz" in executed_cmd
        assert "test.tgz" not in executed_cmd


class TestSkipIfPresent:
    """Content-addressed skip-if-present + dry-run behavior for S3 remotes."""

    pytestmark = pytest.mark.component

    def test_present_object_marked_pushed_without_upload(self, s3_remote_env):
        """A present content-addressed object is marked pushed without re-uploading."""
        from refgenie.cli.commands.push import handle_push

        r, _, asset = s3_remote_env
        s3_remote = _add_s3_remote(r, asset)

        with (
            patch(
                "refgenie.cli.commands.push._s3_head_object",
                return_value={"ContentLength": 1},
            ) as mock_head,
            patch("subprocess.run") as mock_run,
        ):
            handle_push(_make_push_cmd(remote="s3-real"), r)
            mock_run.assert_not_called()
            mock_head.assert_called()

        assert len(_unpushed(r, s3_remote.id)) == 0

    def test_absent_object_uploads(self, s3_remote_env, fake_subprocess):
        """An absent object falls through to a real upload, which is then
        verified before the link is marked pushed."""
        from refgenie.cli.commands.push import handle_push

        r, _, asset = s3_remote_env
        s3_remote = _add_s3_remote(r, asset)

        with (
            patch("refgenie.cli.commands.push._s3_head_object", return_value=None),
            patch("refgenie.cli.commands.push._verify_s3_object", return_value=True),
        ):
            handle_push(_make_push_cmd(remote="s3-real"), r)
            fake_subprocess.assert_called_once()

        assert len(_unpushed(r, s3_remote.id)) == 0

    def test_dry_run_present_makes_no_db_write(self, s3_remote_env):
        """A dry run over an already-present object writes nothing to the DB."""
        from refgenie.cli.commands.push import handle_push

        r, _, asset = s3_remote_env
        s3_remote = _add_s3_remote(r, asset)

        with (
            patch(
                "refgenie.cli.commands.push._s3_head_object",
                return_value={"ContentLength": 1},
            ),
            patch("subprocess.run") as mock_run,
        ):
            handle_push(_make_push_cmd(dry_run=True, remote="s3-real"), r)
            mock_run.assert_not_called()

        links = _unpushed(r, s3_remote.id)
        assert len(links) == 1
        assert links[0].pushed is False


class TestPushVerification:
    """Post-upload verification: an object landing at the right key with the
    wrong size must not be treated as a successful push."""

    pytestmark = pytest.mark.component

    def test_wrong_size_after_upload_not_marked_pushed(self, s3_remote_env, fake_subprocess):
        """The object is absent before the push (so the real upload runs),
        then present at the key afterward but with the wrong size: the link
        must stay unpushed and the command must exit non-zero."""
        from sqlmodel import select

        from refgenie.cli.commands.push import handle_push
        from refgenie.cli.errors import EXIT_GENERAL_ERROR

        r, _, asset = s3_remote_env
        s3_remote = _add_s3_remote(r, asset)

        with r.database._database_session as session:
            staged = session.exec(
                select(StagedAsset).where(
                    StagedAsset.asset_digest == asset.digest, StagedAsset.mode == "archive"
                )
            ).first()
            expected_size = staged.tarball_size
        assert expected_size is not None

        with (
            patch(
                "refgenie.cli.commands.push._s3_head_object",
                side_effect=[None, {"ContentLength": expected_size + 1}],
            ),
            pytest.raises(SystemExit) as excinfo,
        ):
            handle_push(_make_push_cmd(remote="s3-real"), r)

        assert excinfo.value.code == EXIT_GENERAL_ERROR
        links = _unpushed(r, s3_remote.id)
        assert len(links) == 1
        assert links[0].pushed is False


class TestFolderSync:
    """Tests for the folder-sync push strategy."""

    pytestmark = pytest.mark.component

    def test_one_sync_per_remote(self, sync_env, fake_subprocess):
        """Folder sync executes exactly one sync command per remote."""
        from refgenie.cli.commands.push import handle_push

        r, _, _, _ = sync_env
        handle_push(_make_push_cmd(strategy="folder_sync"), r)
        assert fake_subprocess.call_count == 2

    def test_success_marks_all_pushed(self, sync_env, fake_subprocess):
        """A successful sync marks every link for that remote pushed."""
        from refgenie.cli.commands.push import handle_push

        r, _, _, _ = sync_env
        handle_push(_make_push_cmd(strategy="folder_sync"), r)

        assert len(_unpushed(r)) == 0

    def test_failure_leaves_all_unpushed(self, sync_env, failing_subprocess):
        """A failed sync leaves every link for that remote unpushed and the
        command exits non-zero."""
        from refgenie.cli.commands.push import handle_push

        r, _, _, _ = sync_env
        failing_subprocess(subprocess.CalledProcessError(1, "aws s3 sync", stderr="fail"))
        with pytest.raises(SystemExit) as excinfo:
            handle_push(_make_push_cmd(strategy="folder_sync"), r)

        assert excinfo.value.code != 0
        assert len(_unpushed(r)) == 2

    def test_missing_staged_file_not_marked_pushed(self, staged_refgenie, fake_subprocess):
        """A sync that exits 0 must only mark links whose staged file actually
        exists locally: a silently-skipped (missing) file stays unpushed and
        the command exits non-zero."""
        from refgenie.cli.commands.push import handle_push
        from refgenie.utils.staging import staged_archive_path

        r = staged_refgenie
        asset = fasta_asset(r)
        remote = r.remote.add(
            type=RemoteType.s3,
            prefix="bucket-x",
            name="sync-s3-x",
            push_command="aws s3 sync {genome_stage_folder} s3://{prefix}/ --follow-symlinks",
        )
        for mode in ("archive", "file"):
            r.remote.link(ref=remote.id, asset_digest=asset.digest, mode=mode, pushed=False)

        # Delete the staged archive from disk; the file-mode symlink remains.
        genome_digest = r.alias.resolve("rCRSd")
        archive = staged_archive_path(
            pathlib.Path(str(r.genome_stage_folder)), genome_digest, "fasta", asset.digest
        )
        archive.unlink()

        with pytest.raises(SystemExit) as excinfo:
            handle_push(_make_push_cmd(strategy="folder_sync"), r)

        assert excinfo.value.code != 0
        # Only the archive link (missing file) remains unpushed.
        links = _unpushed(r, remote.id)
        assert [link.mode for link in links] == ["archive"]

    def test_command_template_substitution(self, sync_env, fake_subprocess):
        """{genome_stage_folder} and {prefix} are substituted in each sync command."""
        from refgenie.cli.commands.push import handle_push

        r, _, _, _ = sync_env
        handle_push(_make_push_cmd(strategy="folder_sync"), r)

        assert fake_subprocess.call_count == 2
        for call in fake_subprocess.call_args_list:
            executed_cmd = call[0][0]
            assert "{genome_stage_folder}" not in executed_cmd
            assert "{prefix}" not in executed_cmd
            assert "bucket-a" in executed_cmd or "bucket-b" in executed_cmd


# ===========================================================================
# Unit: S3 helpers (no Refgenie instance, no filesystem)
# ===========================================================================


class TestS3Helpers:
    """Unit coverage for the S3 prefix/head-object/verification helpers."""

    def test_parse_s3_prefix(self):
        from refgenie.cli.commands.push import _parse_s3_prefix

        assert _parse_s3_prefix("s3://bucket/assets") == ("bucket", "assets")
        assert _parse_s3_prefix("s3://bucket/a/b/") == ("bucket", "a/b")
        assert _parse_s3_prefix("s3://bucket") == ("bucket", "")
        assert _parse_s3_prefix("my-bucket") is None
        assert _parse_s3_prefix("https://cdn.example.com/x") is None

    @pytest.mark.parametrize(
        "run_kwargs, expected",
        [
            ({"return_value": MagicMock(stdout='{"ContentLength": 42}')}, {"ContentLength": 42}),
            ({"side_effect": subprocess.CalledProcessError(255, "aws")}, None),
            ({"side_effect": FileNotFoundError()}, None),
        ],
        ids=["object_present", "nonzero_exit", "no_aws_cli"],
    )
    def test_s3_head_object(self, run_kwargs, expected):
        """A present object returns its parsed JSON response (which carries
        ContentLength); an absent object, a nonzero exit, or a missing aws
        CLI all return None."""
        from refgenie.cli.commands.push import _s3_head_object

        with patch("subprocess.run", **run_kwargs):
            assert _s3_head_object("bucket", "key") == expected

    def test_verify_s3_object_matching_size(self):
        """A present object whose size matches verifies successfully."""
        from refgenie.cli.commands.push import _verify_s3_object

        with patch(
            "refgenie.cli.commands.push._s3_head_object", return_value={"ContentLength": 100}
        ):
            assert _verify_s3_object("bucket", "key", 100) is True

    def test_verify_s3_object_wrong_size(self):
        """A present object with the wrong size does NOT verify -- a partial
        or stale upload must not be treated as a successful push."""
        from refgenie.cli.commands.push import _verify_s3_object

        with patch(
            "refgenie.cli.commands.push._s3_head_object", return_value={"ContentLength": 99}
        ):
            assert _verify_s3_object("bucket", "key", 100) is False

    def test_verify_s3_object_absent(self):
        """No object at the key at all does not verify."""
        from refgenie.cli.commands.push import _verify_s3_object

        with patch("refgenie.cli.commands.push._s3_head_object", return_value=None):
            assert _verify_s3_object("bucket", "key", 100) is False

    def test_verify_s3_object_no_expected_size_falls_back_to_presence(self):
        """mode="file" staged assets carry no tarball_size; verification
        then falls back to a presence-only check rather than silently
        treating "no expected size" as verified."""
        from refgenie.cli.commands.push import _verify_s3_object

        with patch(
            "refgenie.cli.commands.push._s3_head_object", return_value={"ContentLength": 12345}
        ):
            assert _verify_s3_object("bucket", "key", None) is True


# ===========================================================================
# Component: publish-catalog export / import (publish_catalog.py)
# ===========================================================================


@pytest.fixture
def build_env(tmp_path, fixtures_path):
    """A build-node-like catalog: built + staged asset, one pushed link."""
    r = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    r.database.init(
        genome_folder=tmp_path / "genomes",
        genome_stage_folder=tmp_path / "stage",
    )
    register_fasta(r, fixtures_path)
    r.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=["rCRSd"],
        description="rCRSd genome",
    )
    r.build.run(
        recipe_name="fasta",
        genome_alias="rCRSd",
        asset_group_name="fasta",
        asset_name="test",
    )
    genome_digest = r.alias.resolve("rCRSd")
    asset = r.asset.get(
        genome_digest=genome_digest,
        asset_group_name="fasta",
        asset_name="test",
    )
    r.stage.create(
        asset=asset,
        genome_folder=r.genome_folder,
        genome_stage_folder=r.genome_stage_folder,
    )
    remote = r.remote.add(
        type=RemoteType.s3,
        prefix="my-bucket/assets",
        name="build-side s3",
        push_command="aws s3 cp {local_path} s3://{prefix}/{relative_path}",
    )
    # archive was pushed; file mode never was
    r.remote.link(ref=remote.id, asset_digest=asset.digest, mode="archive", pushed=True)
    r.remote.link(ref=remote.id, asset_digest=asset.digest, mode="file", pushed=False)
    return r, genome_digest, asset


def _artifact_engine(path):
    return create_engine(f"sqlite:///{path}")


class TestExport:
    pytestmark = pytest.mark.component

    def test_export_contains_published_only(self, build_env, tmp_path):
        r, genome_digest, asset = build_env
        dest = tmp_path / "publish.sqlite"
        summary = export_publish_catalog(r, dest, HTTPS_PREFIX)

        assert summary["genome"] == 1
        assert summary["asset"] == 1
        eng = _artifact_engine(dest)
        with eng.connect() as c:
            genomes = c.execute(select(Genome.__table__)).all()
            assert [g.digest for g in genomes] == [genome_digest]

            # aliases materialized from the (store-backed) manager into SQL rows
            aliases = c.execute(select(Alias.__table__)).all()
            assert [(a.name, a.genome_digest) for a in aliases] == [("rCRSd", genome_digest)]

            # only the pushed (archive) pair survives; local path stripped
            assets = c.execute(select(Asset.__table__)).all()
            assert [a.digest for a in assets] == [asset.digest]
            assert assets[0].path is None
            staged = c.execute(select(StagedAsset.__table__)).all()
            assert [(s.asset_digest, s.mode) for s in staged] == [(asset.digest, "archive")]
            assert staged[0].download_count == 0

            # one synthesized https remote; pushed link re-targeted at it
            remotes = c.execute(select(Remote.__table__)).all()
            assert len(remotes) == 1
            assert remotes[0].type == RemoteType.https
            assert remotes[0].prefix == HTTPS_PREFIX
            links = c.execute(select(RemoteAssetLink.__table__)).all()
            assert [(x.remote_id, x.asset_digest, x.mode, x.pushed) for x in links] == [
                (1, asset.digest, "archive", True)
            ]

            stamp = c.execute(select(AlembicVersion.__table__.c.version_num)).scalar()
            assert stamp == TARGET_ALEMBIC_VERSION


class TestImport:
    pytestmark = pytest.mark.component

    @pytest.fixture
    def artifact(self, build_env, tmp_path):
        r, genome_digest, asset = build_env
        dest = tmp_path / "publish.sqlite"
        export_publish_catalog(r, dest, HTTPS_PREFIX)
        return dest, genome_digest, asset

    @pytest.fixture
    def server_engine(self):
        eng = make_engine()
        SQLModel.metadata.create_all(eng)
        return eng

    def test_import_populates_server_catalog(self, artifact, server_engine):
        dest, genome_digest, asset = artifact
        summary = import_publish_catalog(server_engine, str(dest))
        assert summary["genome"] == 1

        with server_engine.connect() as c:
            assert c.execute(select(Genome.__table__.c.digest)).scalars().all() == [genome_digest]
            assert c.execute(select(Alias.__table__.c.name)).scalars().all() == ["rCRSd"]

    def test_import_yields_https_download_url(self, artifact, server_engine):
        """The redirect query the server runs must resolve to the https mirror."""
        dest, genome_digest, asset = artifact
        import_publish_catalog(server_engine, str(dest))

        from sqlmodel import Session

        from refgenie.db.tables import Configuration
        from refgenie.managers import RemoteManager

        with Session(server_engine) as session:
            stage_folder = session.scalar(select(Configuration.genome_stage_folder))
        tarball = pathlib.Path(stage_folder) / genome_digest / "fasta" / f"{asset.digest}.tgz"
        [(_, _, url)] = RemoteManager(server_engine).pushed_urls(asset.digest, {"archive": tarball})
        assert url == f"{HTTPS_PREFIX}/{genome_digest}/fasta/{asset.digest}.tgz"

    def test_reimport_is_idempotent_and_preserves_download_count(self, artifact, server_engine):
        dest, genome_digest, asset = artifact
        import_publish_catalog(server_engine, str(dest))

        # server-side activity between imports
        with server_engine.begin() as c:
            c.execute(
                StagedAsset.__table__.update()
                .where(StagedAsset.__table__.c.asset_digest == asset.digest)
                .values(download_count=7)
            )

        import_publish_catalog(server_engine, str(dest))
        with server_engine.connect() as c:
            staged = c.execute(select(StagedAsset.__table__)).all()
            assert len(staged) == 1  # merged by natural key, not duplicated
            assert staged[0].download_count == 7  # server-owned counter kept
            assert c.execute(select(Genome.__table__)).all().__len__() == 1

    def test_import_refuses_schema_mismatch(self, artifact, server_engine):
        dest, _, _ = artifact
        eng = _artifact_engine(dest)
        with eng.begin() as c:
            c.execute(AlembicVersion.__table__.update().values(version_num="bogus"))
        eng.dispose()
        with pytest.raises(ValueError, match="refusing to import"):
            import_publish_catalog(server_engine, str(dest))

    def test_imported_catalog_serves_genomes_endpoint(self, artifact, tmp_path):
        """End to end: a server app bound to the imported catalog lists the genome."""
        requires_server()

        eng = make_engine()
        rgc = Refgenie(database_engine=eng, suppress_migrations=True)
        rgc.database.init(genome_folder=tmp_path / "server_genomes")

        dest, genome_digest, asset = artifact
        import_publish_catalog(eng, str(dest))

        with make_server_client(rgc) as client:
            resp = client.get("/v4/genomes")
            assert resp.status_code == 200
            digests = [g["digest"] for g in resp.json()["items"]]
            assert digests == [genome_digest]
