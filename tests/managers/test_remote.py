"""RemoteManager (``rgc.remote``): push targets and their asset links.

Covers the one naming rule (id or name), add/upsert/remove, link/unlink,
``unpushed``/``mark_pushed``, status, ``pushed_urls``, and the push intent
``StageManager.create(push_to=...)`` records through it. The push workflow
and the ``refgenie remote`` handlers are in ``tests/cli/test_remote_push.py``.

Every test here stages a real asset, so the module is ``component``.
"""

import pathlib

import pytest

from refgenie.db.tables import RemoteType
from refgenie.exceptions import MissingRemoteError
from refgenie.models import GenomeAlias, GenomeDigest
from refgenie.utils.staging import staged_archive_path
from tests.helpers import fasta_asset

pytestmark = pytest.mark.component

CDN = "https://cdn.example.com/data"


def _unpushed(r, ref=None):
    return [link for link, _, _ in r.remote.unpushed(ref)]


def _stage_paths(r, asset):
    """mode -> the local staged path of the fixture's rCRSd fasta asset."""
    genome_digest = r.alias.resolve(GenomeAlias("rCRSd"))
    stage = pathlib.Path(str(r.genome_stage_folder))
    return genome_digest, {
        "archive": staged_archive_path(stage, genome_digest, "fasta", asset.digest),
        "file": stage / genome_digest / "fasta" / asset.name,
    }


class TestGet:
    """The one naming rule: an int or a digit string is an id; anything else is a name."""

    def test_by_id_digit_string_and_name(self, staged_refgenie):
        r = staged_refgenie
        remote = r.remote.add(name="cdn", type=RemoteType.https, prefix=CDN)
        assert r.remote.get(remote.id).id == remote.id
        assert r.remote.get(str(remote.id)).id == remote.id
        assert r.remote.get("cdn").id == remote.id

    def test_unknown_raises(self, staged_refgenie):
        r = staged_refgenie
        r.remote.add(name="cdn", type=RemoteType.https, prefix=CDN)
        for ref in ("nope", 999, "999"):
            with pytest.raises(MissingRemoteError):
                r.remote.get(ref)
        assert r.remote.exists("cdn")
        assert not r.remote.exists("nope")


class TestAddUpsertRemove:
    def test_multiple_remotes_of_same_type_allowed(self, staged_refgenie):
        r = staged_refgenie
        r1 = r.remote.add(name="CDN 1", type=RemoteType.https, prefix="https://cdn1.example.com")
        r2 = r.remote.add(name="CDN 2", type=RemoteType.https, prefix="https://cdn2.example.com")
        assert r1.id != r2.id
        assert [x.id for x in r.remote.list_all()] == [r1.id, r2.id]

    @pytest.mark.parametrize("name", ["", "42"])
    def test_add_rejects_names_that_are_not_names(self, staged_refgenie, name):
        with pytest.raises(ValueError):
            staged_refgenie.remote.add(name=name, type=RemoteType.s3, prefix="b")

    def test_add_rejects_duplicate_name(self, staged_refgenie):
        r = staged_refgenie
        r.remote.add(name="cdn", type=RemoteType.https, prefix=CDN)
        with pytest.raises(ValueError, match="already exists"):
            r.remote.add(name="cdn", type=RemoteType.s3, prefix="b")

    def test_upsert_creates_then_updates(self, staged_refgenie):
        r = staged_refgenie
        r1 = r.remote.upsert(name="test_cdn", type=RemoteType.https, prefix="https://cdn1.x")
        assert r1.id is not None
        assert r1.name == "test_cdn"
        r2 = r.remote.upsert(
            name="test_cdn",
            type=RemoteType.https,
            prefix="https://cdn2.x",
            push_command="aws s3 cp {local_path} {prefix}/{relative_path}",
        )
        assert r1.id == r2.id
        assert r2.prefix == "https://cdn2.x"
        assert r2.push_command == "aws s3 cp {local_path} {prefix}/{relative_path}"

    def test_remove_by_name_and_id(self, staged_refgenie):
        r = staged_refgenie
        a = r.remote.add(name="a", type=RemoteType.s3, prefix="a")
        r.remote.add(name="b", type=RemoteType.s3, prefix="b")
        r.remote.remove("b")
        r.remote.remove(a.id)
        assert r.remote.list_all() == []
        with pytest.raises(MissingRemoteError):
            r.remote.remove("b")


class TestLinks:
    def test_link_creates_intent(self, staged_refgenie):
        r = staged_refgenie
        remote = r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN)
        asset = fasta_asset(r)
        link, created = r.remote.link("CDN", asset.digest, "file")
        assert created
        assert link.pushed is False
        assert (link.remote_id, link.asset_digest, link.mode) == (remote.id, asset.digest, "file")

    def test_link_requires_staged_asset(self, staged_refgenie):
        r = staged_refgenie
        r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN)
        with pytest.raises(ValueError, match="StagedAsset not found"):
            r.remote.link("CDN", "nonexistent_digest", "archive")

    def test_link_exist_ok_is_idempotent_and_keeps_pushed(self, staged_refgenie):
        r = staged_refgenie
        r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN)
        asset = fasta_asset(r)
        r.remote.link("CDN", asset.digest, "file", pushed=True)
        with pytest.raises(ValueError, match="already linked"):
            r.remote.link("CDN", asset.digest, "file")
        link, created = r.remote.link("CDN", asset.digest, "file", exist_ok=True)
        assert not created
        assert link.pushed is True
        assert r.remote.unpushed() == []

    def test_unpushed_and_mark_pushed_round_trip(self, staged_refgenie):
        r = staged_refgenie
        remote = r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN)
        asset = fasta_asset(r)
        r.remote.link(remote.id, asset.digest, "file")

        [(link, got_remote, staged)] = r.remote.unpushed()
        assert link.asset_digest == asset.digest
        assert got_remote.id == remote.id
        assert (staged.asset_digest, staged.mode) == (asset.digest, "file")
        # Loaded for callers that derive staged paths after the session closes.
        assert staged.asset.asset_group.name == "fasta"

        r.remote.mark_pushed("CDN", asset.digest, "file")
        assert r.remote.unpushed() == []

    def test_unpushed_filtered_by_remote(self, staged_refgenie):
        r = staged_refgenie
        r1 = r.remote.add(name="CDN 1", type=RemoteType.https, prefix="https://cdn1.x")
        r2 = r.remote.add(name="S3", type=RemoteType.s3, prefix="s3://bucket")
        asset = fasta_asset(r)
        r.remote.link(r1.id, asset.digest, "file")
        r.remote.link(r2.id, asset.digest, "file")

        assert len(_unpushed(r, r1.id)) == 1
        assert len(_unpushed(r, "S3")) == 1
        assert len(_unpushed(r)) == 2
        with pytest.raises(MissingRemoteError):
            r.remote.unpushed("nope")

    def test_unpushed_filtered_by_genome(self, staged_refgenie):
        r = staged_refgenie
        r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN)
        asset = fasta_asset(r)
        r.remote.link("CDN", asset.digest, "file")
        genome_digest = r.alias.resolve(GenomeAlias("rCRSd"))
        assert len(r.remote.unpushed(genome_digest=genome_digest)) == 1
        assert r.remote.unpushed(genome_digest=GenomeDigest("A" * 32)) == []

    def test_unlink(self, staged_refgenie):
        r = staged_refgenie
        r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN)
        asset = fasta_asset(r)
        r.remote.link("CDN", asset.digest, "file")
        r.remote.unlink("CDN", asset.digest, "file")
        assert _unpushed(r) == []


class TestStatus:
    def test_no_remotes(self, staged_refgenie):
        assert staged_refgenie.remote.status() == {}

    def test_remote_without_links(self, staged_refgenie):
        r = staged_refgenie
        remote = r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN)
        status = r.remote.status()
        assert status[remote.id]["remote"].name == "CDN"
        assert status[remote.id]["pushed"] == []
        assert status[remote.id]["unpushed"] == []

    def test_groups_pushed_and_unpushed(self, staged_refgenie):
        r = staged_refgenie
        remote = r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN)
        asset = fasta_asset(r)
        r.remote.link(remote.id, asset.digest, "file")
        status = r.remote.status()
        assert (len(status[remote.id]["unpushed"]), len(status[remote.id]["pushed"])) == (1, 0)

        r.remote.mark_pushed(remote.id, asset.digest, "file")
        status = r.remote.status()
        assert (len(status[remote.id]["unpushed"]), len(status[remote.id]["pushed"])) == (0, 1)

    def test_filtered_by_remote(self, staged_refgenie):
        r = staged_refgenie
        r1 = r.remote.add(name="CDN 1", type=RemoteType.https, prefix="https://cdn1.x")
        r2 = r.remote.add(name="S3", type=RemoteType.s3, prefix="s3://bucket")
        assert list(r.remote.status(r1.id)) == [r1.id]
        assert list(r.remote.status("S3")) == [r2.id]
        with pytest.raises(MissingRemoteError):
            r.remote.status("no-such-remote")


class TestPushedUrls:
    """``pushed_urls``: the download URL on each remote an asset was pushed to."""

    @pytest.mark.parametrize("mode", ["archive", "file"])
    def test_pushed_link_gives_prefix_plus_stage_relative_path(self, staged_refgenie, mode):
        r = staged_refgenie
        remote = r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN + "/")
        asset = fasta_asset(r)
        r.remote.link("CDN", asset.digest, mode, pushed=True)
        genome_digest, paths = _stage_paths(r, asset)

        [(got_remote, got_mode, url)] = r.remote.pushed_urls(asset.digest, {mode: paths[mode]})
        assert (got_remote.id, got_mode) == (remote.id, mode)
        relative = paths[mode].relative_to(pathlib.Path(str(r.genome_stage_folder)))
        assert url == f"{CDN}/{relative}"
        assert genome_digest in url
        if mode == "archive":
            assert url.endswith(f"{asset.digest}.tgz")

    def test_unlinked_or_unpushed_gives_nothing(self, staged_refgenie):
        r = staged_refgenie
        r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN)
        asset = fasta_asset(r)
        _, paths = _stage_paths(r, asset)
        assert r.remote.pushed_urls(asset.digest, paths) == []
        r.remote.link("CDN", asset.digest, "file")
        assert r.remote.pushed_urls(asset.digest, paths) == []

    def test_s3_remote_gives_no_url(self, staged_refgenie):
        r = staged_refgenie
        r.remote.add(name="S3", type=RemoteType.s3, prefix="s3://bucket")
        asset = fasta_asset(r)
        r.remote.link("S3", asset.digest, "archive", pushed=True)
        _, paths = _stage_paths(r, asset)
        assert r.remote.pushed_urls(asset.digest, paths) == []

    def test_url_comes_from_the_linked_remote(self, staged_refgenie):
        r = staged_refgenie
        r.remote.add(name="CDN 1", type=RemoteType.https, prefix="https://cdn1.example.com")
        r.remote.add(name="CDN 2", type=RemoteType.https, prefix="https://cdn2.example.com")
        asset = fasta_asset(r)
        r.remote.link("CDN 1", asset.digest, "file", pushed=True)
        _, paths = _stage_paths(r, asset)
        [(_, _, url)] = r.remote.pushed_urls(asset.digest, paths)
        assert "cdn1.example.com" in url

    def test_https_before_http(self, staged_refgenie):
        r = staged_refgenie
        r.remote.add(name="plain", type=RemoteType.http, prefix="http://plain.example.com")
        r.remote.add(name="secure", type=RemoteType.https, prefix="https://secure.example.com")
        asset = fasta_asset(r)
        r.remote.link("plain", asset.digest, "archive", pushed=True)
        r.remote.link("secure", asset.digest, "archive", pushed=True)
        _, paths = _stage_paths(r, asset)
        urls = [url for _, _, url in r.remote.pushed_urls(asset.digest, paths)]
        assert urls[0].startswith("https://secure.example.com/")
        assert len(urls) == 2

    def test_path_outside_stage_folder_is_skipped(self, staged_refgenie, tmp_path):
        r = staged_refgenie
        r.remote.add(name="CDN", type=RemoteType.https, prefix=CDN)
        asset = fasta_asset(r)
        r.remote.link("CDN", asset.digest, "file", pushed=True)
        assert r.remote.pushed_urls(asset.digest, {"file": tmp_path / "elsewhere"}) == []


class TestStagePushTo:
    """``StageManager.create(push_to=...)`` records push intent through ``rgc.remote``."""

    def _restage(self, r, asset, push_to):
        r.stage.create(
            asset=asset,
            genome_folder=r.genome_folder,
            genome_stage_folder=r.genome_stage_folder,
            push_to=push_to,
        )

    @pytest.mark.parametrize("by", ["name", "id"])
    def test_push_to_finds_remote_by_name_or_id(self, staged_refgenie, by):
        """``--push-to`` matches a remote's name or id, not its prefix, and never
        silently skips."""
        r = staged_refgenie
        remote = r.remote.add(name="my-s3", type=RemoteType.s3, prefix="s3://bucket/assets")
        asset = fasta_asset(r)
        self._restage(r, asset, ["my-s3" if by == "name" else str(remote.id)])
        assert {link.asset_digest for link in _unpushed(r, remote.id)} == {asset.digest}

    def test_unknown_remote_is_skipped(self, staged_refgenie):
        r = staged_refgenie
        asset = fasta_asset(r)
        self._restage(r, asset, ["no-such-remote"])  # must not raise
        assert r.remote.unpushed() == []

    def test_restage_with_push_to_is_idempotent(self, staged_refgenie):
        """A persistent build catalog re-stages already-staged assets on every
        nightly run, so ``stage.create(..., push_to=[...])`` is called
        repeatedly for the same (remote, asset, mode). It must not raise a
        RemoteAssetLink primary-key collision, must not duplicate links, and
        must preserve the pushed state of links already marked pushed."""
        r = staged_refgenie
        remote = r.remote.add(
            name="asset-s3",
            type=RemoteType.s3,
            prefix="s3://bucket/assets",
            push_command="aws s3 sync {genome_stage_folder} {prefix}/",
        )
        asset = fasta_asset(r)

        self._restage(r, asset, ["asset-s3"])
        links = r.remote.status(remote.id)[remote.id]["unpushed"]
        assert len(links) >= 1

        for link in links:
            r.remote.mark_pushed(remote.id, link.asset_digest, link.mode)

        self._restage(r, asset, ["asset-s3"])
        status = r.remote.status(remote.id)[remote.id]
        assert len(status["pushed"]) == len(links)
        assert status["unpushed"] == []
