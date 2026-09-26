"""Hooks fire from real refgenie operations: builds, pulls, adds, removes.

Component tier: every test builds or pulls real assets on disk.
"""

import logging
import threading
from pathlib import Path
from unittest.mock import patch

import pytest

from refgenie import Refgenie
from refgenie.config import config
from refgenie.exceptions import AssetExistsError
from refgenie.models import GenomeAlias
from refgenie.plugins import registry
from refgenie.plugins.hooks import POST_BUILD, POST_PULL, POST_UPDATE, PRE_BUILD, PRE_PULL
from tests.helpers import (
    ASSET,
    GENOME,
    GROUP,
    build_rcrsd,
    fasta_asset,
    requires_dash,
    requires_server,
    stage_copy,
)
from tests.plugins import fake_plugin
from tests.plugins.fake_plugin import BOOM, RECORD, on_every_hook

pytestmark = pytest.mark.component


def hooks() -> list[str]:
    return [hook for hook, _ in fake_plugin.CALLS]


def events(hook: str) -> list:
    return [event for h, event in fake_plugin.CALLS if h == hook]


def actions(event) -> list[str]:
    return [change.action for change in event.changes]


def only_update():
    """The single post_update fired, asserting there was exactly one."""
    updates = events(POST_UPDATE)
    assert len(updates) == 1, f"expected one post_update, got {hooks()}"
    return updates[0]


# ---------------------------------------------------------------------------
# Build, add and the local-state mutators
# ---------------------------------------------------------------------------


def test_build_fires_pre_post_build_then_one_update(refgenie_fs, install_plugins):
    install_plugins(on_every_hook())
    build_rcrsd(refgenie_fs)
    assert hooks() == [PRE_BUILD, POST_BUILD, POST_UPDATE]
    (post_build,) = events(POST_BUILD)
    assert post_build.succeeded is True
    assert (post_build.genome, post_build.asset_group, post_build.asset) == (GENOME, GROUP, ASSET)
    added = [c for c in only_update().changes if c.action == "asset_added"]
    assert [(c.asset_group, c.asset) for c in added] == [(GROUP, ASSET)]


def test_content_add_fires_one_update_after_the_alias_tree_exists(refgenie_built, install_plugins):
    rg = refgenie_built
    _, rel = stage_copy(rg, fasta_asset(rg), "staging", mutate=b">x\nACGT\n")
    install_plugins(
        {POST_UPDATE: [("a_record", RECORD), ("b_seek", "tests.plugins.fake_plugin:seek_added")]}
    )
    rg.asset.content.add(
        "fasta", rel, GROUP, genome_digest=rg.alias.resolve(GENOME), asset_name="second"
    )
    assert actions(only_update()) == ["asset_added"]
    ((_, seeked),) = [call for call in fake_plugin.CALLS if call[0] == "seek"]
    assert Path(seeked).exists()


@pytest.mark.parametrize(
    "operation, expected",
    [
        (
            lambda rg, d: rg.asset.remove(GROUP, ASSET, genome_digest=d, keep_asset_group=True),
            "asset_removed",
        ),
        (lambda rg, d: rg.asset.rename(GROUP, ASSET, "renamed", genome_digest=d), "asset_renamed"),
        (lambda rg, d: rg.set_genome_alias(GenomeAlias("other"), genome_digest=d), "alias_added"),
    ],
    ids=["remove", "rename", "set_genome_alias"],
)
def test_each_mutator_fires_one_update(refgenie_built, install_plugins, operation, expected):
    rg = refgenie_built
    digest = rg.alias.resolve(GENOME)
    install_plugins(on_every_hook())
    operation(rg, digest)
    assert actions(only_update()) == [expected]


def test_batch_updates_groups_operations_into_one_update(refgenie_built, install_plugins):
    """What the CLI does around each command, e.g. `genome init --fasta`."""
    rg = refgenie_built
    digest = rg.alias.resolve(GENOME)
    install_plugins(on_every_hook())
    with rg.batch_updates():
        rg.set_genome_alias(GenomeAlias("other"), genome_digest=digest)
        rg.asset.rename(GROUP, ASSET, "renamed", genome_digest=digest)
    assert actions(only_update()) == ["alias_added", "asset_renamed"]


def test_alias_remove_fires_one_update(refgenie_built, install_plugins):
    rg = refgenie_built
    rg.set_genome_alias(GenomeAlias("other"), genome_digest=rg.alias.resolve(GENOME))
    install_plugins(on_every_hook())
    rg.alias.remove(GenomeAlias("other"))
    update = only_update()
    assert [(c.action, c.alias) for c in update.changes] == [("alias_removed", "other")]


def test_genome_remove_fires_one_update(refgenie_built, install_plugins):
    """The genome's aliases go with it, so they are in the same update."""
    rg = refgenie_built
    digest = rg.alias.resolve(GENOME)
    install_plugins(on_every_hook())
    rg.genome.remove(digest)
    update = only_update()
    assert [c.action for c in update.changes] == ["genome_removed", "alias_removed"]
    assert update.changes[0].genome == digest


def test_set_default_fires_only_when_the_default_changes(refgenie_built, install_plugins):
    rg = refgenie_built
    digest = rg.alias.resolve(GENOME)
    build_rcrsd(rg, asset_name="other")  # same content, a second name; now the default
    assert rg.asset.group.get_default(GROUP, genome_digest=digest) == "other"
    install_plugins(on_every_hook())

    rg.asset.group.set_default(GROUP, "other", genome_digest=digest)
    assert fake_plugin.CALLS == []

    rg.asset.group.set_default(GROUP, ASSET, genome_digest=digest)
    (change,) = only_update().changes
    assert (change.action, change.asset, change.previous) == ("default_changed", ASSET, "other")


# ---------------------------------------------------------------------------
# Failures and the off switches
# ---------------------------------------------------------------------------


def test_a_failing_plugin_never_fails_the_build(refgenie_fs, install_plugins, caplog):
    install_plugins({POST_UPDATE: [("a_boom", BOOM), ("b_record", RECORD)]})
    with caplog.at_level(logging.WARNING, logger="refgenie"):
        asset = build_rcrsd(refgenie_fs)
    assert asset is not None
    assert fasta_asset(refgenie_fs).digest == asset.digest
    assert "refgenie plugin 'a_boom' failed in post_update" in caplog.text
    assert hooks() == [POST_UPDATE]


def test_disable_all_fires_nothing(refgenie_fs, install_plugins, disable_plugins):
    install_plugins(on_every_hook())
    disable_plugins("1")
    build_rcrsd(refgenie_fs)
    assert fake_plugin.CALLS == []


def test_plugins_false_fires_nothing(refgenie_built, install_plugins):
    rg = Refgenie(
        database_engine=refgenie_built.database_engine, suppress_migrations=True, plugins=False
    )
    install_plugins(on_every_hook())
    rg.asset.rename(GROUP, ASSET, "renamed", genome_digest=rg.alias.resolve(GENOME))
    assert fake_plugin.CALLS == []
    assert rg.plugins.enabled is False


def test_disable_by_name_skips_only_that_plugin(
    refgenie_fs, install_plugins, disable_plugins, caplog
):
    install_plugins({POST_UPDATE: [("boom", BOOM), ("record", RECORD)]})
    disable_plugins("boom")
    with caplog.at_level(logging.WARNING, logger="refgenie"):
        build_rcrsd(refgenie_fs)
    assert hooks() == [POST_UPDATE]
    assert "failed in post_update" not in caplog.text


class TestServerMode:
    """Plugins are off in server mode unless the server opts in."""

    @staticmethod
    def _server(rg, **kwargs):
        return Refgenie(
            database_engine=rg.database_engine,
            suppress_migrations=True,
            server_mode=True,
            **kwargs,
        )

    @staticmethod
    def _rename(rg):
        rg.asset.rename(GROUP, ASSET, "renamed", genome_digest=rg.alias.resolve(GENOME))

    def test_off_by_default(self, refgenie_built, install_plugins):
        server = self._server(refgenie_built)
        install_plugins(on_every_hook())
        self._rename(server)
        assert fake_plugin.CALLS == []

    def test_env_opt_in(self, refgenie_built, install_plugins, monkeypatch):
        monkeypatch.setattr(config, "server_plugins", True)
        server = self._server(refgenie_built)
        install_plugins(on_every_hook())
        self._rename(server)
        assert actions(only_update()) == ["asset_renamed"]

    def test_library_opt_in(self, refgenie_built, install_plugins):
        server = self._server(refgenie_built, plugins=True)
        install_plugins(on_every_hook())
        self._rename(server)
        assert actions(only_update()) == ["asset_renamed"]

    def test_disable_beats_opt_in(self, refgenie_built, install_plugins, disable_plugins):
        server = self._server(refgenie_built, plugins=True)
        install_plugins(on_every_hook())
        disable_plugins("1")
        self._rename(server)
        assert fake_plugin.CALLS == []


def test_reads_never_scan_for_plugins(refgenie_built, monkeypatch):
    def no_scan():
        raise AssertionError("entry points were scanned on a read")

    monkeypatch.setattr(registry, "entry_points", no_scan)
    registry.clear_cache()
    try:
        rg = Refgenie(database_engine=refgenie_built.database_engine, suppress_migrations=True)
        digest = rg.alias.resolve(GENOME)
        rg.asset.seek(digest, GROUP, ASSET)
        rg.asset.list_all()
        rg.paths().to_dict()
    finally:
        registry.clear_cache()


# ---------------------------------------------------------------------------
# Pulls: single, repeated, bulk, and from a dash job thread
# ---------------------------------------------------------------------------


class TestPull:
    requires_server()

    @staticmethod
    def _pull(client, url):
        return client.transfer.pull(
            asset_group_name=GROUP,
            genome=GenomeAlias(GENOME),
            force_large=True,
            force_server_urls=[url],
        )

    def test_pull_then_pull_again(self, server_client_world, install_plugins):
        client, _, url = server_client_world
        install_plugins(on_every_hook())

        assert self._pull(client, url) is not None
        assert hooks() == [PRE_PULL, POST_PULL, POST_UPDATE]
        assert events(POST_PULL)[0].succeeded is True
        assert "asset_added" in actions(only_update())

        fake_plugin.CALLS.clear()
        with pytest.raises(AssetExistsError):
            self._pull(client, url)
        assert hooks() == [PRE_PULL, POST_PULL]
        assert events(POST_PULL)[0].succeeded is True

    def test_pull_genomes_fires_per_asset_and_one_update(
        self, server_client_world, install_plugins, fixtures_path
    ):
        client, server, url = server_client_world
        # A second genome on the server, so the bulk pull has two assets.
        server.genome.initialize_genome(
            fasta_file_path=fixtures_path / "t7.fa", alias_names=["t7"], description="t7"
        )
        t7 = server.build.run(
            recipe_name="fasta", genome_alias="t7", asset_group_name=GROUP, asset_name="default"
        )
        server.stage.create(t7, server.genome_folder, server.genome_stage_folder)
        client.servers.subscribe(url)
        install_plugins(on_every_hook())

        # As in TestInitFromRemote: genomes resolve through the server client.
        with patch(
            "refgenie.managers.sources.servers.make_source", side_effect=ConnectionError("no")
        ):
            pulled = client.transfer.pull_genomes(
                [GenomeAlias(GENOME), GenomeAlias("t7")], force=True, force_large=True
            )

        assert len(pulled) == 2
        assert hooks() == [PRE_PULL, POST_PULL, PRE_PULL, POST_PULL, POST_UPDATE]
        changes = only_update().changes
        added = {c.genome for c in changes if c.action == "asset_added"}
        assert added == {client.alias.resolve(GENOME), client.alias.resolve("t7")}
        genomes_added = {c.genome for c in changes if c.action == "genome_added"}
        assert genomes_added == added

    def test_a_dash_pull_job_fires_from_its_worker_thread(
        self, server_client_world, install_plugins
    ):
        requires_dash()
        from refgenie.server.jobs.manager import JobManager
        from refgenie.server.jobs.schemas import JobStatus, PullJobParams

        client, _, url = server_client_world
        install_plugins(on_every_hook())
        manager = JobManager(client)
        try:
            ref = manager.submit_pull(
                PullJobParams(server_url=url, genome_name=GENOME, asset_group_name=GROUP)
            )
            record = manager.wait(ref.job_id, timeout=120)
        finally:
            manager.shutdown(wait=False)

        assert record.status == JobStatus.SUCCEEDED
        assert hooks() == [PRE_PULL, POST_PULL, POST_UPDATE]
        assert threading.current_thread().name not in fake_plugin.THREADS
