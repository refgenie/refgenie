"""Local-mode remote browsing (``refgenie/server/routers/remote.py``).

``server_url`` must name a subscription: these are GET routes with no action
header guard that an allowlisted bridge origin can read, so an arbitrary URL
would let a public page aim this machine at any host. "No outbound request" is
checked as "no client constructed", since building a client is what fetches
the remote ``/openapi.json``.
"""

import logging
from unittest.mock import MagicMock, patch

import pytest

from tests.helpers import (
    PUBLIC_ORIGIN,
    fake_digest,
    make_local_app,
    make_server_rgc,
    requires_dash,
    requires_server,
    serve_refgenie,
)

requires_server()
requires_dash()

from fastapi.testclient import TestClient  # noqa: E402  (must follow the extras guard)

SUBSCRIBED = "http://s.example"
OTHER = "http://t.example"
DIGEST_S = fake_digest("remote_s")
DIGEST_T = fake_digest("remote_t")
CLIENT_PATH = "refgenie.managers.sources.servers.RefgenieserverClient"


@pytest.fixture
def rgc(tmp_path, fixtures_path):
    """The local refgenie whose /v1/remote routes are under test."""
    return make_server_rgc(tmp_path / "local", fixtures_path)


@pytest.fixture
def client(rgc):
    # base_url: the local app's Host-header guard admits loopback names only.
    with TestClient(
        make_local_app(rgc), base_url="http://localhost", raise_server_exceptions=False
    ) as tc:
        yield tc


@pytest.fixture
def remote_s(tmp_path, fixtures_path):
    return make_server_rgc(tmp_path / "s", fixtures_path, genomes=[(DIGEST_S, ["s_genome"])])


@pytest.fixture
def remote_t(tmp_path, fixtures_path):
    return make_server_rgc(tmp_path / "t", fixtures_path, genomes=[(DIGEST_T, ["t_genome"])])


UNSUBSCRIBED_REQUESTS = [
    ("/v1/remote/genomes", {"server_url": "http://169.254.169.254"}),
    (
        "/v1/remote/assets",
        {"genome_digest": str(DIGEST_S), "server_url": "http://169.254.169.254"},
    ),
]


class TestUnsubscribedServerRefused:
    @pytest.mark.parametrize("origin", [None, PUBLIC_ORIGIN])
    @pytest.mark.parametrize("path,params", UNSUBSCRIBED_REQUESTS)
    def test_refused_before_any_client_is_built(self, rgc, client, path, params, origin):
        rgc.servers.subscribe(SUBSCRIBED)
        headers = {"Origin": origin} if origin else {}
        with patch(CLIENT_PATH, MagicMock()) as make_client:
            response = client.get(path, params=params, headers=headers)
        assert response.status_code == 404
        assert response.json()["error"]["code"] == "server_not_subscribed"
        assert "169.254" not in response.text
        make_client.assert_not_called()
        assert rgc.servers.clients == {}


class TestSubscribedServers:
    def test_subscribed_url_with_trailing_slash(self, rgc, client, remote_s):
        rgc.servers.subscribe(SUBSCRIBED)
        with serve_refgenie(rgc, remote_s, SUBSCRIBED):
            response = client.get("/v1/remote/genomes", params={"server_url": SUBSCRIBED + "/"})
        assert response.status_code == 200
        rows = response.json()
        assert [row["genome_digest"] for row in rows] == [DIGEST_S]
        assert all(row["server_url"] == SUBSCRIBED for row in rows)

    def test_no_server_url_lists_every_subscription(self, rgc, client, remote_s, remote_t):
        rgc.servers.subscribe([SUBSCRIBED, OTHER])
        with serve_refgenie(rgc, remote_s, SUBSCRIBED), serve_refgenie(rgc, remote_t, OTHER):
            response = client.get("/v1/remote/genomes")
        assert response.status_code == 200
        assert {row["genome_digest"] for row in response.json()} == {DIGEST_S, DIGEST_T}


class TestErrorsStayInTheLog:
    SECRET = "secret-internal-host:5432 refused"

    @pytest.mark.parametrize(
        "path,params,method",
        [
            ("/v1/remote/genomes", {}, "list_genomes"),
            ("/v1/remote/assets", {"genome_digest": str(DIGEST_S)}, "list_assets_for_genome"),
        ],
    )
    def test_generic_502(self, rgc, client, monkeypatch, caplog, path, params, method):
        def boom(*args, **kwargs):
            raise RuntimeError(self.SECRET)

        monkeypatch.setattr(rgc.servers, method, boom)
        with caplog.at_level(logging.ERROR, logger="refgenie.server.routers.remote"):
            response = client.get(path, params=params)
        assert response.status_code == 502
        assert response.json()["error"]["code"] == "remote_unavailable"
        assert "secret-internal-host" not in response.text
        assert "secret-internal-host" in caplog.text

    def test_unreachable_server_reports_fixed_error(self, rgc, client):
        rgc.servers.subscribe(SUBSCRIBED)
        with patch(CLIENT_PATH, side_effect=ConnectionError(self.SECRET)):
            response = client.get("/v1/remote/servers")
        assert response.status_code == 200
        [server] = response.json()["servers"]
        assert server["reachable"] is False
        assert server["error"] == "unreachable"
        assert "secret-internal-host" not in response.text
