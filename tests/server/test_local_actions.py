"""The local actions API ``/v1/actions/*`` (``refgenie/server/local/``).

The HTTP envelope, 202 job submission, the action-header / CORS / host-guard
security stack, server-mode absence, build preflight, the synchronous curation
endpoints over a real Refgenie, and the route-inventory contract.

Tier: ``unit`` by default. The classes that build real genome folders on disk
carry a class-level ``pytest.mark.component``.
"""

import threading
from unittest.mock import patch

import pytest

from tests.helpers import (
    EVIL_ORIGIN,
    PUBLIC_ORIGIN,
    TESTS_DATA_DIR,
    act,
    make_local_app,
    make_local_client,
    make_server_app,
    preflight,
    requires_dash,
    requires_server,
    route_keys,
    stub_rgc,
    web_stub_rgc,
)

requires_server()
requires_dash()


#: Every state-changing route, with a minimal valid body. The security tests
#: parametrize over this list so a newly added route is covered by adding a row.
ROUTES = [
    ("POST", "/v1/actions/pull", {"asset_group": "fasta", "genome": "rCRSd"}),
    ("POST", "/v1/actions/build", {"recipe": "fasta", "genome": "rCRSd", "asset_group": "fasta"}),
    (
        "POST",
        "/v1/actions/build/preflight",
        {"recipe": "fasta", "genome": "rCRSd", "asset_group": "fasta"},
    ),
    ("POST", "/v1/actions/genomes", {"fasta": "/tmp/genome.fa", "aliases": ["g1"]}),
    ("DELETE", "/v1/actions/assets/deadbeef", None),
    ("DELETE", "/v1/actions/genomes/rCRSd", None),
    ("POST", "/v1/actions/aliases", {"alias": "a2", "genome_digest": "d1"}),
    ("DELETE", "/v1/actions/aliases/a2", None),
    ("POST", "/v1/actions/subscriptions", {"server_urls": ["http://s.example"]}),
    ("DELETE", "/v1/actions/subscriptions", {"server_urls": ["http://s.example"]}),
    (
        "POST",
        "/v1/actions/assets/default",
        {"genome_digest": "d1", "asset_group": "fasta", "asset": "test"},
    ),
]


@pytest.mark.parametrize(
    "env_var,value",
    [
        ("REFGENIE_LOCAL_ALLOWED_ORIGINS", '["*"]'),
        ("REFGENIE_BRIDGE_ORIGINS", "*"),
        ("REFGENIE_BRIDGE_ORIGIN_REGEX", ".*"),
    ],
)
def test_wildcard_origin_config_refuses_to_construct(monkeypatch, env_var, value):
    """``allow_origins=['*']`` (or a ``.*`` regex) would defeat the whole design
    in local mode: fail at app-construction time, whichever knob sets it."""
    monkeypatch.setenv(env_var, value)
    with pytest.raises(ValueError, match="allowlist"):
        make_local_app(web_stub_rgc())


class TestActionSecurity:
    """Header, CORS allowlist, host guard, and the server-mode absence."""

    @pytest.fixture()
    def client(self):
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as c:
            yield c

    @pytest.mark.parametrize("method,path,body", ROUTES)
    def test_missing_action_header_is_403_with_envelope(self, client, method, path, body):
        response = client.request(method, path, json=body)
        assert response.status_code == 403, f"{method} {path}"
        payload = response.json()
        assert payload["ok"] is False
        assert payload["error"]["code"] == "missing_action_header"

    @pytest.mark.parametrize("method,path,body", ROUTES)
    def test_with_header_the_request_reaches_the_handler(self, client, method, path, body):
        response = act(client, method, path, json=body)
        assert response.status_code != 403, f"{method} {path}: {response.text}"

    def test_preflight_for_allowed_origin_allows_the_action_header(self, client):
        response = preflight(
            client,
            "/v1/actions/pull",
            PUBLIC_ORIGIN,
            "POST",
            **{"Access-Control-Request-Headers": "x-refgenie-action"},
        )
        assert response.status_code == 200
        assert response.headers["access-control-allow-origin"] == PUBLIC_ORIGIN
        assert "x-refgenie-action" in response.headers["access-control-allow-headers"].lower()

    def test_preflight_for_delete_is_refused(self, client):
        """Destructive verbs are not on the cross-origin surface at all."""
        response = preflight(client, "/v1/actions/assets/deadbeef", PUBLIC_ORIGIN, "DELETE")
        assert response.status_code == 400

    def test_preflight_for_unlisted_origin_gets_no_cors_grant(self, client):
        response = preflight(client, "/v1/actions/pull", EVIL_ORIGIN, "POST")
        assert "access-control-allow-origin" not in response.headers

    def test_origin_allowlist_is_configurable(self, monkeypatch):
        # Bridge off so this isolates REFGENIE_LOCAL_ALLOWED_ORIGINS: with the
        # bridge on (default "read"), refgenie.org would be granted anyway.
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "off")
        monkeypatch.setenv("REFGENIE_LOCAL_ALLOWED_ORIGINS", '["https://example.test"]')
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            granted = preflight(client, "/v1/actions/pull", "https://example.test", "POST")
            refused = preflight(client, "/v1/actions/pull", PUBLIC_ORIGIN, "POST")
        assert granted.headers["access-control-allow-origin"] == "https://example.test"
        assert "access-control-allow-origin" not in refused.headers

    def test_non_loopback_host_is_421_forbidden_host(self, client):
        """The DNS-rebinding guard: a rebound hostname arrives as a foreign Host
        header and must be refused before any handler runs."""
        response = client.get("/service-info", headers={"host": "evil.example"})
        assert response.status_code == 421
        payload = response.json()
        assert payload["ok"] is False
        assert payload["error"]["code"] == "forbidden_host"

    @pytest.mark.parametrize("host", ["127.0.0.1", "localhost:8080", "[::1]:9000"])
    def test_loopback_hosts_pass_the_guard(self, client, host):
        assert client.get("/service-info", headers={"host": host}).status_code == 200

    def test_server_mode_has_no_actions_routes(self):
        """Security-critical: several request models accept server-local
        filesystem paths, so the router must never reach the public server."""
        app = make_server_app(stub_rgc())
        actions_paths = [p for _, p in route_keys(app) if p.startswith("/v1/actions")]
        assert actions_paths == []


class TestPullEndpoint:
    """The submission contract only -- error mapping lives with the runners."""

    def test_pull_returns_202_jobref_and_typed_params(self):
        rgc = web_stub_rgc()
        with make_local_client(rgc) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/pull",
                json={"asset_group": "fasta", "genome": "rCRSd"},
            )
            assert response.status_code == 202
            ref = response.json()
            assert ref["kind"] == "pull"
            assert ref["status"] == "queued"
            assert ref["duplicate"] is False
            assert ref["links"]["events"] == "/v1/jobs/events"
            record = client.app.state.job_manager.wait(ref["job_id"])
        assert record.params == {
            "server_url": None,
            "genome_name": "rCRSd",
            "genome_digest": None,
            "asset_group_name": "fasta",
            "asset_name": None,
            "force": False,
        }

    def test_duplicate_submission_coalesces_with_202(self):
        rgc = web_stub_rgc()
        release = threading.Event()
        asset = rgc.transfer.pull.return_value
        rgc.transfer.pull.side_effect = lambda **kwargs: (release.wait(10), asset)[1]
        body = {"asset_group": "fasta", "genome": "rCRSd"}
        with make_local_client(rgc) as client:
            first = act(client, "POST", "/v1/actions/pull", json=body).json()
            second_response = act(client, "POST", "/v1/actions/pull", json=body)
            release.set()
            client.app.state.job_manager.wait(first["job_id"])
        assert second_response.status_code == 202  # there is no 409 job_in_progress
        second = second_response.json()
        assert second["job_id"] == first["job_id"]
        assert second["duplicate"] is True

    @pytest.mark.parametrize("server_url", [None, "http://mock-server"])
    def test_root_call_contract(self, server_url):
        """server_url -> force_server_urls, and the invisible half: the runner
        must pass force_large=True and an explicit confirmer, or a pull can
        block a worker thread on a prompt nobody sees."""
        rgc = web_stub_rgc()
        body = {"asset_group": "fasta", "genome": "rCRSd"}
        if server_url:
            body["server_url"] = server_url
        with make_local_client(rgc) as client:
            ref = act(client, "POST", "/v1/actions/pull", json=body).json()
            client.app.state.job_manager.wait(ref["job_id"])
        kwargs = rgc.transfer.pull.call_args.kwargs
        assert kwargs["force_server_urls"] == ([server_url] if server_url else None)
        assert kwargs["force_large"] is True
        assert kwargs["confirm"] is not None

    @pytest.mark.parametrize(
        "body",
        [
            {"asset_group": "fasta"},  # neither genome nor genome_digest
            {"asset_group": "fasta", "genome": "g", "genome_digest": "d"},  # both
            {"asset_group": "fasta", "genome": "g", "bogus_field": 1},  # extra=forbid
        ],
    )
    def test_invalid_bodies_are_422_with_envelope(self, body):
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            response = act(client, "POST", "/v1/actions/pull", json=body)
        assert response.status_code == 422
        payload = response.json()
        assert payload["ok"] is False
        assert payload["error"]["code"] == "validation_error"


class TestBuildEndpoint:
    def test_build_returns_202_and_converts_params(self):
        from refgenie.models import BuildParams

        rgc = web_stub_rgc()
        body = {
            "recipe": "bwa_index",
            "genome": "rCRSd",
            "asset_group": "bwa_index",
            "pull_parents": True,
            "params": {"files": {"data": "/tmp/data.txt"}, "params": {"cores": 4}},
        }
        with make_local_client(rgc) as client:
            response = act(client, "POST", "/v1/actions/build", json=body)
            assert response.status_code == 202
            ref = response.json()
            assert ref["kind"] == "build"
            record = client.app.state.job_manager.wait(ref["job_id"])
        assert record.status == "succeeded"
        kwargs = rgc.build.run.call_args.kwargs
        assert kwargs["pull_parents"] is True
        assert isinstance(kwargs["params"], BuildParams)
        assert kwargs["params"].params == {"cores": 4}
        assert str(kwargs["params"].files["data"]) == "/tmp/data.txt"

    def test_preflight_reaches_the_root_and_submits_no_job(self):
        """The router must use ``rgc.build.preflight``, never a private
        ``BuildManager`` step -- patching the stub proves the path."""
        rgc = web_stub_rgc()
        rgc.build.preflight.return_value = {
            "ok": True,
            "errors": [],
            "resolved": {"genome_digest": "genomedigest123", "asset_name": "default"},
        }
        with make_local_client(rgc) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/build/preflight",
                json={"recipe": "fasta", "genome": "rCRSd", "asset_group": "fasta"},
            )
            jobs = client.get("/v1/jobs").json()["items"]
        assert response.status_code == 200
        payload = response.json()
        assert payload["ok"] is True
        assert payload["resolved"]["genome_digest"] == "genomedigest123"
        rgc.build.preflight.assert_called_once()
        assert rgc.build.preflight.call_args.kwargs["recipe_name"] == "fasta"
        assert jobs == []

    def test_preflight_reports_field_scoped_errors_with_200(self):
        rgc = web_stub_rgc()
        rgc.build.preflight.return_value = {
            "ok": False,
            "errors": [
                {"field": "params", "code": "missing_build_input", "message": "Missing 'data'."}
            ],
            "resolved": {},
        }
        with make_local_client(rgc) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/build/preflight",
                json={"recipe": "needsfile", "genome": "rCRSd", "asset_group": "grp"},
            )
        assert response.status_code == 200  # a preflight that found problems succeeded
        payload = response.json()
        assert payload["ok"] is False
        assert payload["errors"][0]["field"] == "params"
        assert payload["errors"][0]["code"] == "missing_build_input"


class TestGenomeInitEndpoint:
    def test_genome_init_returns_202_and_runs_initialize_and_build(self):
        rgc = web_stub_rgc()
        body = {"fasta": "/tmp/genome.fa", "aliases": ["mygenome"], "species": "H. testens"}
        with make_local_client(rgc) as client:
            response = act(client, "POST", "/v1/actions/genomes", json=body)
            assert response.status_code == 202
            ref = response.json()
            assert ref["kind"] == "genome_init"
            record = client.app.state.job_manager.wait(ref["job_id"])
        assert record.status == "succeeded"
        assert record.result.genome_digest == "genomedigest123"
        kwargs = rgc.build.initialize_and_build.call_args.kwargs
        assert kwargs["genome_names"] == ["mygenome"]
        assert kwargs["species_name"] == "H. testens"


class TestSyncActions:
    """Real Refgenie, real filesystem, no mocks."""

    pytestmark = pytest.mark.component

    @pytest.fixture()
    def world(self, refgenie_built):
        with make_local_client(refgenie_built) as client:
            yield client, refgenie_built

    def _built_digest(self, rg):
        genome_digest = rg.alias.resolve("rCRSd")
        return genome_digest, rg.asset.get(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name="test"
        ).digest

    def test_delete_asset_removes_it(self, world):
        """DELETE on an asset removes it from the catalog."""
        client, rg = world
        genome_digest, asset_digest = self._built_digest(rg)
        response = act(client, "DELETE", f"/v1/actions/assets/{asset_digest}")
        assert response.status_code == 200
        payload = response.json()
        assert payload["ok"] is True
        assert payload["data"]["registry_path"]
        assert not rg.asset.exists(
            genome_digest=genome_digest, asset_group_name="fasta", asset_name="test"
        )

    def test_delete_asset_unknown_digest_is_404(self, world):
        client, _ = world
        response = act(client, "DELETE", "/v1/actions/assets/no_such_digest")
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "asset_not_found"

    def test_delete_asset_with_children_is_409_naming_them(self, world, fixtures_path):
        client, rg = world
        genome_digest, parent_digest = self._built_digest(rg)
        rg.genome.initialize_genome(
            fasta_file_path=fixtures_path / "demo.fa",
            alias_names=["demo"],
            description="demo genome",
        )
        rg.build.run(
            recipe_name="fasta", genome_alias="demo", asset_group_name="fasta", asset_name="test"
        )
        child = rg.asset.get(
            genome_digest=rg.alias.resolve("demo"), asset_group_name="fasta", asset_name="test"
        )
        rg.asset.links.set_children(genome_digest, "fasta", "test", [child.digest])

        response = act(client, "DELETE", f"/v1/actions/assets/{parent_digest}")
        assert response.status_code == 409
        error = response.json()["error"]
        assert error["code"] == "conflict"
        assert f"{rg.alias.resolve('demo')}/fasta" in error["message"]

    def test_delete_genome_by_digest_only(self, world):
        """The route takes a digest. An alias there is not a digest: 404."""
        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        response = act(client, "DELETE", "/v1/actions/genomes/rCRSd")
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "genome_not_found"
        assert rg.genome.exists(digest)

        response = act(client, "DELETE", f"/v1/actions/genomes/{digest}")
        assert response.status_code == 200
        assert response.json()["data"]["genome_digest"] == digest
        assert not rg.genome.exists(digest)

    @pytest.mark.parametrize("ref", ["never_heard_of_it", "c" * 32], ids=["malformed", "unknown"])
    def test_delete_genome_unknown_is_404(self, world, ref):
        client, _ = world
        response = act(client, "DELETE", f"/v1/actions/genomes/{ref}")
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "genome_not_found"

    def test_alias_set_resolves_afterwards(self, world):
        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        response = act(
            client,
            "POST",
            "/v1/actions/aliases",
            json={"alias": "rCRSd2", "genome_digest": digest},
        )
        assert response.status_code == 200
        assert rg.alias.resolve("rCRSd2") == digest

    def test_alias_set_for_unknown_digest_is_404_not_a_phantom_genome(self, world):
        """Regression guard for the phantom-genome bug: without the existence
        guard, set_genome_alias silently creates a genome row for any string."""
        client, rg = world
        ghost = "0000000000_no_such_genome"
        response = act(
            client,
            "POST",
            "/v1/actions/aliases",
            json={"alias": "ghost", "genome_digest": ghost},
        )
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "genome_not_found"
        assert not rg.genome.exists(ghost)

    def test_alias_remove(self, world):
        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        rg.set_genome_alias(alias_name="doomed", genome_digest=digest)
        response = act(client, "DELETE", "/v1/actions/aliases/doomed")
        assert response.status_code == 200
        response = act(client, "DELETE", "/v1/actions/aliases/doomed")
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "alias_not_found"

    def test_subscribe_reset_and_unsubscribe(self, world):
        client, rg = world
        response = act(
            client,
            "POST",
            "/v1/actions/subscriptions",
            json={"server_urls": ["http://a.example"]},
        )
        assert response.status_code == 200
        assert "http://a.example" in response.json()["data"]["subscriptions"]

        response = act(
            client,
            "POST",
            "/v1/actions/subscriptions",
            json={"server_urls": ["http://b.example"], "reset": True},
        )
        assert response.json()["data"]["subscriptions"] == ["http://b.example"]

        response = act(
            client,
            "DELETE",
            "/v1/actions/subscriptions",
            json={"server_urls": ["http://b.example"]},
        )
        assert response.status_code == 200
        assert "http://b.example" not in response.json()["data"]["subscriptions"]
        assert rg.servers.subscriptions() == []

    def test_data_channel_add_sync_and_remove(self, world, monkeypatch):
        """The channel row is real; the network is not. The sync loop is fed by
        the sources manager, so stubbing its three fetchers stands in for a
        reachable index that publishes one asset class and one recipe."""
        client, rg = world
        fixtures = TESTS_DATA_DIR
        asset_class = fixtures / "fasta_asset_class.yaml"
        recipe = fixtures / "fasta_asset_recipe.yaml"
        monkeypatch.setattr(rg.sources, "test_channel", lambda name: True)
        monkeypatch.setattr(rg.sources, "iter_asset_classes", lambda name: iter([str(asset_class)]))
        monkeypatch.setattr(rg.sources, "iter_recipes", lambda name: iter([str(recipe)]))

        response = act(
            client,
            "POST",
            "/v1/actions/data_channels",
            json={"name": "chan", "index_address": "https://example.org/index.yaml"},
        )
        assert response.status_code == 200, response.text
        payload = response.json()
        assert "not verified" in payload["message"]
        assert payload["data"]["trusted"] is False
        # `fasta` is already registered by the fixture, so exists_ok skips it.
        assert payload["data"]["sync"]["asset_classes_skipped"] == 1
        assert payload["data"]["sync"]["recipes_skipped"] == 1
        assert payload["data"]["sync"]["errors"] == []
        assert rg.sources.get_channel("chan").type.value == "https"

        listed = client.get("/v1/remote/data_channels").json()["channels"]
        assert [c["name"] for c in listed] == ["chan"]
        assert listed[0]["trusted"] is False
        assert listed[0]["credentials_set"] is False

        response = act(client, "POST", "/v1/actions/data_channels/chan/sync")
        assert response.status_code == 200
        assert response.json()["data"]["sync"]["channel"] == "chan"

        response = act(
            client,
            "POST",
            "/v1/actions/data_channels",
            json={"name": "chan", "index_address": "https://example.org/index.yaml"},
        )
        assert response.status_code == 409

        response = act(client, "DELETE", "/v1/actions/data_channels/chan")
        assert response.status_code == 200
        assert rg.sources.get_channel("chan") is None
        response = act(client, "DELETE", "/v1/actions/data_channels/chan")
        assert response.status_code == 404

    @pytest.mark.parametrize(
        "index_address", ["/etc/passwd", "file:///etc/passwd", "ftp://host/index.yaml"]
    )
    def test_data_channel_add_takes_only_http_urls(self, world, index_address):
        """Over HTTP a local path is an arbitrary local-file read."""
        client, rg = world
        response = act(
            client,
            "POST",
            "/v1/actions/data_channels",
            json={"name": "chan", "index_address": index_address},
        )
        assert response.status_code == 422
        assert "http" in response.json()["error"]["message"]
        assert rg.sources.get_channel("chan") is None

    def test_data_channel_add_without_sync_only_records_it(self, world, monkeypatch):
        client, rg = world
        monkeypatch.setattr(
            rg.sources, "test_channel", lambda name: pytest.fail("sync must not run")
        )
        response = act(
            client,
            "POST",
            "/v1/actions/data_channels",
            json={"name": "chan", "index_address": "http://example.org/index.yaml", "sync": False},
        )
        assert response.status_code == 200
        assert "sync" not in response.json()["data"]
        assert rg.sources.get_channel("chan").type.value == "http"

    def test_data_channel_sync_unreachable_is_502(self, world, monkeypatch):
        client, rg = world
        rg.sources.add_channel(name="chan", type="https", index_address="https://x/index.yaml")
        monkeypatch.setattr(rg.sources, "test_channel", lambda name: False)
        response = act(client, "POST", "/v1/actions/data_channels/chan/sync")
        assert response.status_code == 502
        assert "not accessible" in response.json()["error"]["message"]

    def test_set_default_asset(self, world):
        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        response = act(
            client,
            "POST",
            "/v1/actions/assets/default",
            json={"genome_digest": digest, "asset_group": "fasta", "asset": "test"},
        )
        assert response.status_code == 200
        assert rg.asset.group.get_default("fasta", genome_digest=digest) == "test"

    def test_set_default_asset_unknown_name_is_404(self, world):
        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        response = act(
            client,
            "POST",
            "/v1/actions/assets/default",
            json={"genome_digest": digest, "asset_group": "fasta", "asset": "nope"},
        )
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "asset_not_found"

    def test_genome_gone_but_alias_left_says_so(self, world):
        """The drift ``GenomeManager.remove`` leaves when it is interrupted
        between the catalog commit and the store cleanup: the aliases page still
        lists the name and links it here. A bare ``genome_not_found`` would tell
        the user nothing about which half of their instance is wrong."""
        from sqlmodel import Session, select

        from refgenie.db.tables import Genome

        client, rg = world
        digest = rg.alias.resolve("rCRSd")
        # The catalog half only -- the store keeps naming the digest.
        with Session(rg.database_engine) as session:
            session.delete(session.exec(select(Genome).where(Genome.digest == digest)).one())
            session.commit()
        assert "rCRSd" in rg.alias.get_for_genome(digest)

        response = client.get(f"/v4/genomes/{digest}")
        assert response.status_code == 404
        error = response.json()["error"]
        assert error["code"] == "stale_alias"
        assert "rCRSd" in error["message"]
        assert "refgenie alias remove rCRSd" in error["message"]

    def test_genome_never_known_stays_a_plain_404(self, world):
        """No alias points at it, so there is no drift to explain: the generic
        404 this route has always raised, which the UI renders as the bare
        identifier."""
        client, _ = world
        response = client.get("/v4/genomes/0000000000_no_such_genome")
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "not_found"


class TestPreflightReal:
    """``rgc.build.preflight`` over a real world (no HTTP mocking)."""

    pytestmark = pytest.mark.component

    def test_valid_build_preflights_ok(self, refgenie_fs):
        with make_local_client(refgenie_fs) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/build/preflight",
                json={"recipe": "fasta", "genome": "rCRSd", "asset_group": "fasta"},
            )
        assert response.status_code == 200
        payload = response.json()
        assert payload["ok"] is True, payload
        assert payload["errors"] == []
        assert payload["resolved"]["genome_digest"] == refgenie_fs.alias.resolve("rCRSd")
        assert payload["resolved"]["asset_name"] == "default"

    def test_unknown_genome_and_recipe_are_field_scoped(self, refgenie_fs):
        with make_local_client(refgenie_fs) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/build/preflight",
                json={"recipe": "no_such_recipe", "genome": "no_such_genome", "asset_group": "x"},
            )
        payload = response.json()
        assert response.status_code == 200
        assert payload["ok"] is False
        by_field = {error["field"]: error["code"] for error in payload["errors"]}
        assert by_field["genome"] == "genome_not_found"
        assert by_field["recipe"] == "recipe_not_found"

    def test_missing_required_file_is_field_scoped(self, refgenie_fs, tmp_path):
        recipe_yaml = tmp_path / "needsfile_recipe.yaml"
        recipe_yaml.write_text(
            "\n".join(
                [
                    "name: needsfile",
                    "version: 0.1.0",
                    "output_asset_class: fasta",
                    "description: test recipe requiring an input file",
                    "input_files:",
                    "  data:",
                    "    description: required data file",
                    "input_params: null",
                    "input_assets: null",
                    "docker_image: null",
                    "command_templates:",
                    "  - cp {{values.files.data}} {{values.output_folder}}/",
                    'default_asset: "default"',
                ]
            )
        )
        refgenie_fs.recipe.add(recipe_yaml)
        with make_local_client(refgenie_fs) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/build/preflight",
                json={"recipe": "needsfile", "genome": "rCRSd", "asset_group": "grp"},
            )
        payload = response.json()
        assert response.status_code == 200
        assert payload["ok"] is False
        assert any(
            e["field"] == "params" and e["code"] == "missing_build_input" and "data" in e["message"]
            for e in payload["errors"]
        ), payload["errors"]


PROBE = "refgenie.managers.build.check_output"
BOWTIE2_IMAGE = "docker.io/databio/refgenie"


def tool_only_inside_the_image(image, version):
    """A host like a real server: the recipe's tool exists only in its image.

    Anything not a `docker run` of `image` gets the ENOENT a missing binary
    would give, so a probe that reaches for the host is indistinguishable from
    one on a machine where the tool was never installed.
    """

    def probe(command, *args, **kwargs):
        if isinstance(command, list) and command[:2] == ["docker", "run"] and image in command:
            return version.encode()
        raise FileNotFoundError(2, "No such file or directory: 'bowtie2-build'")

    return probe


class TestPreflightNamesContainerRecipes:
    """Preflight must resolve an asset name for a recipe whose tools live in a
    container -- 15 of the 29 stock recipes (component tier).

    Those recipes name their asset after a tool version read by a shell
    one-liner, and declare the `docker_image` that one-liner needs. Preflight
    must run it in that image: run on the host, any machine without bowtie2
    installed -- which is every server -- gets an empty name back, and the
    build form refuses the build before the user can start it.
    """

    pytestmark = pytest.mark.component

    @staticmethod
    def _with_bowtie2_recipe(rgc, fixtures_path):
        rgc.asset_class.add(fixtures_path / "bowtie2_index_asset_class.yaml")
        rgc.recipe.add(fixtures_path / "bowtie2_index_asset_recipe.yaml")
        return rgc

    def _preflight(self, rgc, recipe="bowtie2_index", group="bowtie2_index"):
        with make_local_client(rgc) as client:
            return act(
                client,
                "POST",
                "/v1/actions/build/preflight",
                json={"recipe": recipe, "genome": "rCRSd", "asset_group": group},
            )

    def test_asset_name_is_read_from_the_recipes_image(self, refgenie_fs, fixtures_path):
        rgc = self._with_bowtie2_recipe(refgenie_fs, fixtures_path)
        with patch(PROBE, side_effect=tool_only_inside_the_image(BOWTIE2_IMAGE, "2.3.0\n")):
            payload = self._preflight(rgc).json()
        assert payload["resolved"]["asset_name"] == "2.3.0", payload
        assert not any(e["field"] == "asset" for e in payload["errors"]), payload["errors"]

    def test_unusable_docker_is_a_field_problem_not_a_crash(self, refgenie_fs, fixtures_path):
        """Preflight is advisory: the form still lets the submit through and
        lets the job report the failure. So say plainly what could not be run."""
        rgc = self._with_bowtie2_recipe(refgenie_fs, fixtures_path)
        with patch(PROBE, side_effect=FileNotFoundError(2, "No such file: 'docker'")):
            response = self._preflight(rgc)
        assert response.status_code == 200
        payload = response.json()
        assert payload["ok"] is False
        problem = next(e for e in payload["errors"] if e["field"] == "asset")
        assert BOWTIE2_IMAGE in problem["message"], problem
        assert "asset_name" not in payload["resolved"]

    def test_an_unnameable_asset_asks_for_a_name_under_its_own_code(
        self, refgenie_fs, fixtures_path
    ):
        """A probe that prints nothing is not a broken build: it is a field
        refgenie could not fill in. The form reads ``asset_name_required`` as
        "make Asset name required", and the same request WITH a name is ok."""
        rgc = self._with_bowtie2_recipe(refgenie_fs, fixtures_path)
        with patch(PROBE, side_effect=tool_only_inside_the_image(BOWTIE2_IMAGE, "")):
            payload = self._preflight(rgc).json()
            problem = next(e for e in payload["errors"] if e["field"] == "asset")
            assert problem["code"] == "asset_name_required", problem
            assert "asset_name" not in payload["resolved"]

            # The fixture genome has no fasta asset, so `ok` stays False for
            # that reason; what matters here is that the NAME problem is gone.
            with make_local_client(rgc) as client:
                named = act(
                    client,
                    "POST",
                    "/v1/actions/build/preflight",
                    json={
                        "recipe": "bowtie2_index",
                        "genome": "rCRSd",
                        "asset_group": "bowtie2_index",
                        "asset": "my_name",
                    },
                ).json()
        assert not any(e["field"] == "asset" for e in named["errors"]), named["errors"]
        assert named["resolved"]["asset_name"] == "my_name"

    def test_a_recipe_named_default_runs_no_command_at_all(self, refgenie_fs):
        """`fasta` declares no custom seek keys, so nothing is probed and no
        container starts -- on a host with no docker at all, it still says
        `default`."""
        with patch(PROBE, side_effect=AssertionError("preflight ran a command")) as probe:
            payload = self._preflight(refgenie_fs, recipe="fasta", group="fasta").json()
        probe.assert_not_called()
        assert payload["ok"] is True, payload
        assert payload["resolved"]["asset_name"] == "default"


class TestActionContract:
    """The surface is exactly ``EXPECTED`` -- a silently added or dropped
    endpoint is a build failure, the lesson of the rotted v1 tree."""

    EXPECTED = {
        ("POST", "/v1/actions/pull"),
        ("POST", "/v1/actions/build"),
        ("POST", "/v1/actions/build/preflight"),
        ("POST", "/v1/actions/genomes"),
        ("DELETE", "/v1/actions/assets/{asset_digest}"),
        ("DELETE", "/v1/actions/genomes/{genome_digest}"),
        ("POST", "/v1/actions/aliases"),
        ("DELETE", "/v1/actions/aliases/{alias_name}"),
        ("POST", "/v1/actions/subscriptions"),
        ("DELETE", "/v1/actions/subscriptions"),
        ("POST", "/v1/actions/data_channels"),
        ("DELETE", "/v1/actions/data_channels/{name}"),
        ("POST", "/v1/actions/data_channels/{name}/sync"),
        ("POST", "/v1/actions/assets/default"),
    }

    def _actions_routes(self):
        app = make_local_app(stub_rgc())
        return {
            (method, path)
            for method, path in route_keys(app)
            if path.startswith("/v1/actions") and method != "HEAD"
        }

    def test_route_inventory_matches_the_contract(self):
        assert self._actions_routes() == self.EXPECTED
