"""The localhost bridge: ``refgenie/server/routers/ping.py`` and the bridge
half of ``refgenie/server/local/security.py``.

``/ping`` and its contract, bridge modes, the bridge CORS allowlist, the
local-network-access preflight, and the cross-origin action policy.
"""

import pytest
from starlette.requests import Request

from refgenie.server.local.security import LocalSecuritySettings, is_bridge_request
from tests.helpers import (
    CAPABILITY_KEY_SET,
    EVIL_ORIGIN,
    PUBLIC_ORIGIN,
    act,
    make_local_client,
    make_server_client,
    preflight,
    requires_dash,
    requires_server,
    stub_rgc,
    web_stub_rgc,
)

requires_server()
requires_dash()


REQUIRED_PING_FIELDS = {
    "service",
    "bridge_version",
    "mode",
    "refgenie_version",
    "api_version",
    "instance_id",
    "instance_label",
    "bridge_mode",
    "action_header",
    "capabilities",
    "bridge",
}


class TestPingContract:
    @pytest.mark.parametrize(
        # bridge_mode: "read" is the local default; the public API has no bridge.
        "mode,bridge_mode,pull,archives",
        [("local", "read", True, False), ("server", "off", False, True)],
    )
    def test_ping_shape(self, mode, bridge_mode, pull, archives):
        make_client = make_local_client if mode == "local" else make_server_client
        rgc = web_stub_rgc() if mode == "local" else stub_rgc()
        with make_client(rgc) as client:
            response = client.get("/ping")
        assert response.status_code == 200
        assert response.headers["cache-control"] == "no-store"
        payload = response.json()
        assert REQUIRED_PING_FIELDS <= set(payload)
        assert payload["service"] == "refgenie"
        assert isinstance(payload["bridge_version"], int)
        assert payload["action_header"] == "X-Refgenie-Action"
        assert payload["mode"] == mode
        assert payload["bridge_mode"] == bridge_mode
        assert payload["bridge"] == {"actions_cross_origin": False}
        assert set(payload["capabilities"]) == CAPABILITY_KEY_SET
        assert payload["capabilities"]["pull"] is pull
        assert payload["capabilities"]["archives"] is archives

    def test_ping_capabilities_match_service_info(self):
        """One shared vocabulary: /ping must emit exactly what /service-info
        emits, with no bridge-specific renames."""
        with make_local_client(web_stub_rgc()) as client:
            ping = client.get("/ping").json()
            service_info = client.get("/service-info").json()
        assert ping["capabilities"] == service_info["refgenie"]["capabilities"]

    def test_ping_reports_full_mode(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "full")
        with make_local_client(web_stub_rgc()) as client:
            payload = client.get("/ping").json()
        assert payload["bridge_mode"] == "full"
        assert payload["bridge"] == {"actions_cross_origin": True}

    def test_ping_still_answers_under_bridge_off(self, monkeypatch):
        """Same-origin callers (the local SPA) use /ping too; off only means no
        cross-origin caller can *read* it."""
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "off")
        with make_local_client(web_stub_rgc()) as client:
            payload = client.get("/ping").json()
        assert payload["bridge_mode"] == "off"

    def test_ping_omits_paths_by_default(self, tmp_path):
        rgc = web_stub_rgc()
        rgc.genome_folder = tmp_path / "genomes"
        with make_local_client(rgc) as client:
            payload = client.get("/ping").json()
        assert payload["instance_label"] == "local refgenie"
        assert str(tmp_path) not in str(payload)

    def test_ping_includes_paths_when_opted_in(self, monkeypatch, tmp_path):
        monkeypatch.setenv("REFGENIE_BRIDGE_EXPOSE_PATHS", "true")
        rgc = web_stub_rgc()
        rgc.genome_folder = tmp_path / "genomes"
        with make_local_client(rgc) as client:
            payload = client.get("/ping").json()
        assert payload["instance_label"] == str(tmp_path / "genomes")

    def test_instance_id_is_stable_across_app_constructions(self, monkeypatch, tmp_path):
        monkeypatch.setenv("REFGENIE_HOME_PATH", str(tmp_path))
        ids = []
        for _ in range(2):
            with make_local_client(web_stub_rgc()) as client:
                ids.append(client.get("/ping").json()["instance_id"])
        assert ids[0] == ids[1]
        assert (tmp_path / "instance_id").read_text().strip() == ids[0]


class TestBridgeCors:
    def test_preflight_from_bridge_origin_is_granted_by_default(self):
        with make_local_client(web_stub_rgc()) as client:
            response = preflight(client, "/ping", PUBLIC_ORIGIN)
        assert response.status_code == 200
        assert response.headers["access-control-allow-origin"] == PUBLIC_ORIGIN

    def test_preflight_from_unlisted_origin_gets_no_grant(self):
        with make_local_client(web_stub_rgc()) as client:
            response = preflight(client, "/ping", EVIL_ORIGIN)
        assert "access-control-allow-origin" not in response.headers

    def test_actual_read_carries_cors_headers_for_bridge_origin(self):
        """Job-status polling is a read; under ``read`` mode the bridge origin
        must be able to read /ping, /v4 and /v1/jobs responses."""
        with make_local_client(web_stub_rgc()) as client:
            for path in ("/ping", "/v1/jobs"):
                response = client.get(path, headers={"Origin": PUBLIC_ORIGIN})
                assert response.status_code == 200, path
                assert response.headers.get("access-control-allow-origin") == PUBLIC_ORIGIN, path

    def test_bridge_off_produces_no_cors_for_the_public_origin(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "off")
        with make_local_client(web_stub_rgc()) as client:
            options = preflight(client, "/ping", PUBLIC_ORIGIN)
            actual = client.get("/ping", headers={"Origin": PUBLIC_ORIGIN})
            same_origin = client.get("/ping")
        assert "access-control-allow-origin" not in options.headers
        assert "access-control-allow-origin" not in actual.headers
        assert same_origin.status_code == 200

    def test_origin_regex_escape_hatch(self, monkeypatch):
        monkeypatch.setenv(
            "REFGENIE_BRIDGE_ORIGIN_REGEX", r"https://[a-z0-9-]+\.refgenie-ui\.pages\.dev"
        )
        with make_local_client(web_stub_rgc()) as client:
            granted = client.get(
                "/ping", headers={"Origin": "https://deadbeef.refgenie-ui.pages.dev"}
            )
            refused = client.get("/ping", headers={"Origin": EVIL_ORIGIN})
        assert (
            granted.headers.get("access-control-allow-origin")
            == "https://deadbeef.refgenie-ui.pages.dev"
        )
        assert "access-control-allow-origin" not in refused.headers


class TestLocalNetworkAccessPreflight:
    @pytest.mark.parametrize(
        "request_header,response_header",
        [
            ("Access-Control-Request-Private-Network", "access-control-allow-private-network"),
            ("Access-Control-Request-Local-Network", "access-control-allow-local-network"),
        ],
    )
    def test_lna_grant_for_allowlisted_origin(self, request_header, response_header):
        """Both header spellings (PNA-era and LNA-era) are feature-detected. The
        status assertion is load-bearing: Starlette 400s any preflight carrying
        Access-Control-Request-Private-Network unless CORSMiddleware is built
        with ``allow_private_network=True``."""
        with make_local_client(web_stub_rgc()) as client:
            response = preflight(client, "/ping", PUBLIC_ORIGIN, **{request_header: "true"})
        assert response.status_code == 200
        assert response.headers.get("access-control-allow-origin") == PUBLIC_ORIGIN
        assert response.headers.get(response_header) == "true"

    def test_lna_grant_withheld_for_unlisted_origin(self):
        with make_local_client(web_stub_rgc()) as client:
            response = preflight(
                client,
                "/ping",
                EVIL_ORIGIN,
                **{"Access-Control-Request-Private-Network": "true"},
            )
        assert response.status_code == 400
        assert "access-control-allow-origin" not in response.headers

    def test_lna_middleware_absent_under_bridge_off(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "off")
        with make_local_client(web_stub_rgc()) as client:
            response = preflight(
                client,
                "/ping",
                PUBLIC_ORIGIN,
                **{"Access-Control-Request-Private-Network": "true"},
            )
        assert "access-control-allow-private-network" not in response.headers


PULL_BODY = {"asset_group": "fasta", "genome": "rCRSd"}


class TestCrossOriginActionPolicy:
    @pytest.mark.parametrize("bridge_mode", ["off", "read", "full"])
    def test_pull_without_action_header_is_403_in_every_mode(self, monkeypatch, bridge_mode):
        """The bridge must never weaken the anti-CSRF header requirement."""
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", bridge_mode)
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            response = client.post("/v1/actions/pull", json=PULL_BODY)
        assert response.status_code == 403
        assert response.json()["error"]["code"] == "missing_action_header"

    def test_cross_origin_pull_is_403_under_read_with_remedy(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "read")
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/pull",
                json=PULL_BODY,
                headers={"Origin": PUBLIC_ORIGIN},
            )
        assert response.status_code == 403
        payload = response.json()
        assert payload["error"]["code"] == "forbidden_origin"
        assert "refgenie dash --bridge full" in payload["error"]["message"]

    def test_cross_origin_pull_is_accepted_under_full(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "full")
        with make_local_client(web_stub_rgc()) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/pull",
                json=PULL_BODY,
                headers={"Origin": PUBLIC_ORIGIN},
            )
            assert response.status_code == 202
            client.app.state.job_manager.wait(response.json()["job_id"])

    def test_same_origin_pull_is_unaffected_by_bridge_mode(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "off")
        with make_local_client(web_stub_rgc()) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/pull",
                json=PULL_BODY,
                headers={"Origin": "http://localhost"},  # == the request's own origin
            )
            assert response.status_code == 202
            client.app.state.job_manager.wait(response.json()["job_id"])

    @pytest.mark.parametrize("bridge_mode", ["off", "read", "full"])
    def test_cross_origin_delete_is_403_in_every_mode(self, monkeypatch, bridge_mode):
        """Destructive verbs are never on the cross-origin surface."""
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", bridge_mode)
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            response = act(
                client, "DELETE", "/v1/actions/assets/deadbeef", headers={"Origin": PUBLIC_ORIGIN}
            )
        assert response.status_code == 403
        assert response.json()["error"]["code"] == "forbidden_origin"

    @pytest.mark.parametrize("bridge_mode", ["read", "full"])
    def test_unlisted_origin_actions_are_403(self, monkeypatch, bridge_mode):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", bridge_mode)
        with make_local_client(web_stub_rgc(), raise_server_exceptions=False) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/pull",
                json=PULL_BODY,
                headers={"Origin": EVIL_ORIGIN},
            )
        assert response.status_code == 403
        assert response.json()["error"]["code"] == "forbidden_origin"

    def test_dev_origin_keeps_full_action_access(self):
        """The Vite dev server serves the *local* SPA; it is trusted like
        same-origin and is not subject to the bridge's pull-only policy."""
        with make_local_client(web_stub_rgc()) as client:
            response = act(
                client,
                "POST",
                "/v1/actions/pull",
                json=PULL_BODY,
                headers={"Origin": "http://localhost:5173"},
            )
            assert response.status_code == 202
            client.app.state.job_manager.wait(response.json()["job_id"])


class TestBridgePullSubscriptionCheck:
    """A bridge-origin pull may only name a subscribed server; same-origin may name any.

    A bridge page cannot subscribe, so the subscription list is the only place
    the user said which servers they trust. A same-origin caller could
    subscribe first anyway, and the local /pull page is how a user confirms a
    pull from an unsubscribed server.
    """

    SUBSCRIBED = "http://s.example"

    @pytest.fixture
    def rgc(self, monkeypatch):
        monkeypatch.setenv("REFGENIE_BRIDGE_MODE", "full")
        rgc = web_stub_rgc()
        rgc.servers.subscriptions.return_value = [self.SUBSCRIBED]
        return rgc

    def _pull(self, client, server_url, origin):
        body = {**PULL_BODY}
        if server_url is not None:
            body["server_url"] = server_url
        return act(client, "POST", "/v1/actions/pull", json=body, headers={"Origin": origin})

    def test_bridge_pull_from_unsubscribed_server_is_403(self, rgc):
        with make_local_client(rgc, raise_server_exceptions=False) as client:
            response = self._pull(client, "http://169.254.169.254", PUBLIC_ORIGIN)
            assert response.status_code == 403
            assert response.json()["error"]["code"] == "server_not_subscribed"
            assert client.app.state.job_manager.list() == []
        rgc.transfer.pull.assert_not_called()

    def test_bridge_pull_from_subscribed_server_uses_stored_form(self, rgc):
        with make_local_client(rgc) as client:
            response = self._pull(client, self.SUBSCRIBED + "/", PUBLIC_ORIGIN)
            assert response.status_code == 202
            client.app.state.job_manager.wait(response.json()["job_id"])
        assert rgc.transfer.pull.call_args.kwargs["force_server_urls"] == [self.SUBSCRIBED]

    def test_bridge_pull_without_server_url_is_unchanged(self, rgc):
        with make_local_client(rgc) as client:
            response = self._pull(client, None, PUBLIC_ORIGIN)
            assert response.status_code == 202
            client.app.state.job_manager.wait(response.json()["job_id"])
        assert rgc.transfer.pull.call_args.kwargs["force_server_urls"] is None

    def test_same_origin_pull_accepts_any_server_url(self, rgc):
        """Deliberate split: the local SPA can subscribe anyway."""
        with make_local_client(rgc) as client:
            response = self._pull(client, "http://other.example/", "http://localhost")
            assert response.status_code == 202
            client.app.state.job_manager.wait(response.json()["job_id"])
        assert rgc.transfer.pull.call_args.kwargs["force_server_urls"] == ["http://other.example/"]


class _App:
    def __init__(self, settings):
        self.state = type("State", (), {"local_security": settings})()


def _request(origin):
    headers = [(b"host", b"localhost")]
    if origin is not None:
        headers.append((b"origin", origin.encode()))
    return Request(
        {
            "type": "http",
            "method": "POST",
            "scheme": "http",
            "path": "/v1/actions/pull",
            "headers": headers,
            "query_string": b"",
            "server": ("localhost", 80),
            "app": _App(LocalSecuritySettings()),
        }
    )


@pytest.mark.parametrize(
    "origin,expected",
    [
        (None, False),
        ("http://localhost", False),
        (PUBLIC_ORIGIN, True),
        ("http://localhost:5173", False),
    ],
)
def test_is_bridge_request(origin, expected):
    assert is_bridge_request(_request(origin)) is expected
