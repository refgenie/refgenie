"""The app factory and the SPA: ``refgenie/server/main.py`` and ``spa.py``.

``create_app(mode=...)`` construction, mode isolation, the shared JSON
surface, the ``/service-info`` bootstrap, the root namespace reserved for the
SPA, route/operationId uniqueness and the SPA-route/API non-collision guard,
SPA serving (index, hashed assets, cache headers, deep-link fallback,
missing-bundle 503), and the OpenAPI <-> TypeScript / SKILL.md drift guards.

Apps are always built through the real ``create_app`` factory (via the
``make_*_app`` helpers), never a hand-assembled ``FastAPI()``; the one
exception is the ``openapi_schema`` fixture, which includes the routers on a
bare app to read their schema.
"""

import re
from pathlib import Path

import pytest

from tests.helpers import (
    CAPABILITY_KEY_SET,
    assert_no_duplicate_routes,
    assert_unique_operation_ids,
    fake_digest,
    make_local_app,
    make_server_app,
    make_server_client,
    make_server_rgc,
    requires_dash,
    requires_server,
    requires_web_assets,
    route_keys,
    stub_rgc,
)

WEB_DIGEST = fake_digest("web_digest")

requires_server()
requires_dash()

from fastapi import FastAPI  # noqa: E402  (must follow the extras guard)
from fastapi.testclient import TestClient  # noqa: E402

from refgenie.server.const import (  # noqa: E402
    APP_MODE_LOCAL,
    APP_MODE_SERVER,
    SPA_CLIENT_ROUTES,
)
from refgenie.server.main import CAPABILITY_KEYS, _capabilities, create_app  # noqa: E402
from refgenie.server.routers import catalog as catalog_router  # noqa: E402
from refgenie.server.routers import version4 as version4_router  # noqa: E402

#: The repository root, for the frontend files the drift guards read.
_REPO_ROOT = Path(__file__).resolve().parents[2]


# ===========================================================================
# One factory, two modes: create_app(mode="server"|"local")
# ===========================================================================


@pytest.fixture
def rgc(tmp_path, fixtures_path):
    """A catalog holding one genome, shared by the app-modes and SPA tests."""
    return make_server_rgc(tmp_path, fixtures_path, genomes=[(WEB_DIGEST, ["rCRSd"])])


@pytest.fixture
def server_client(rgc):
    with TestClient(make_server_app(rgc), raise_server_exceptions=True) as client:
        yield client


@pytest.fixture
def local_client(rgc):
    # base_url: the local app's Host-header guard admits loopback names only,
    # so the default "testserver" host would 421 everything.
    with TestClient(
        make_local_app(rgc), base_url="http://localhost", raise_server_exceptions=True
    ) as client:
        yield client


class TestAppConstruction:
    """Both modes build, and neither reaches for the caller's real config."""

    @pytest.mark.parametrize("mode", [APP_MODE_SERVER, APP_MODE_LOCAL])
    def test_construction_does_not_touch_the_real_home(self, mode, rgc, tmp_path, monkeypatch):
        monkeypatch.setenv("HOME", str(tmp_path / "fake_home"))
        app = create_app(mode=mode, refgenie_instance=rgc)
        assert app.state.mode == mode
        assert not (tmp_path / "fake_home").exists()

    def test_unknown_mode_is_rejected(self, rgc):
        with pytest.raises(ValueError):
            create_app(mode="dashboard", refgenie_instance=rgc)


class TestSharedSurface:
    """The JSON API is the same in both modes, at the same prefix."""

    @pytest.mark.parametrize("client_name", ["server_client", "local_client"])
    def test_v4_genomes_lists_the_seeded_genome(self, client_name, request):
        client = request.getfixturevalue(client_name)
        response = client.get("/v4/genomes")
        assert response.status_code == 200
        assert [g["digest"] for g in response.json()["items"]] == [WEB_DIGEST]


class TestServiceInfoBootstrap:
    """/service-info is how one SPA bundle serves both modes."""

    def test_reports_its_mode(self, server_client, local_client):
        assert server_client.get("/service-info").json()["refgenie"]["mode"] == "server"
        assert local_client.get("/service-info").json()["refgenie"]["mode"] == "local"

    def test_ga4gh_fields_stay_at_the_top_level(self, server_client):
        data = server_client.get("/service-info").json()
        assert data["id"] == "org.refgenie.api"
        assert data["type"]["group"] == "org.refgenie"
        assert "mode" not in data

    def test_local_mode_omits_the_seqcol_block(self, local_client, server_client):
        assert "seqcol" not in local_client.get("/service-info").json()
        assert "seqcol" in server_client.get("/service-info").json()

    def test_capability_key_set_is_identical_across_modes(self, server_client, local_client):
        """One vocabulary, one shape: the SPA reads the same keys in both modes."""
        server_caps = server_client.get("/service-info").json()["refgenie"]["capabilities"]
        local_caps = local_client.get("/service-info").json()["refgenie"]["capabilities"]
        assert set(server_caps) == set(local_caps) == CAPABILITY_KEY_SET
        assert all(isinstance(v, bool) for v in server_caps.values())

    #: Commands are local-only, bulk data is server-only, and the two
    #: definition writes are off everywhere -- the modes are mirror images.
    COMMAND_KEYS = ("pull", "build", "delete", "jobs", "remote_browse", "data_channels")
    DATA_KEYS = ("downloads", "archives", "seqcol", "drs")

    @pytest.mark.parametrize("mode,commands_on", [(APP_MODE_SERVER, False), (APP_MODE_LOCAL, True)])
    def test_capability_matrix(self, mode, commands_on):
        caps = _capabilities(mode)
        assert all(caps[key] is commands_on for key in self.COMMAND_KEYS)
        assert all(caps[key] is not commands_on for key in self.DATA_KEYS)
        assert not caps["recipes_write"] and not caps["asset_classes_write"]


#: (path, status in server mode, status in local mode).
_MODE_ISOLATION = [
    ("/v1/remote/genomes", 404, 200),
    ("/v1/remote/servers", 404, 200),
    ("/v1/remote/data_channels", 404, 200),
    ("/v1/jobs", 404, 200),
    ("/ga4gh/drs/service-info", 200, 404),
    ("/v4/ga4gh/drs/service-info", 200, 404),
    ("/seqcol/service-info", 200, 404),
    ("/data_channel/", 200, 404),
    ("/v4/archives", 200, 404),
]


@pytest.mark.parametrize("path,server_status,local_status", _MODE_ISOLATION)
def test_mode_isolation(server_client, local_client, path, server_status, local_status):
    # /v1/actions/* is absent from this table on purpose: the actions router has
    # no GET routes and the SPA catch-all answers GET for unmatched paths in
    # both modes, so a GET probe cannot tell the modes apart. Its mode isolation
    # is asserted by route inventory (TestActionContract,
    # TestActionSecurity.test_server_mode_has_no_actions_routes) and the
    # combined server-mode surface check below.
    assert server_client.get(path).status_code == server_status, f"server mode: {path}"
    assert local_client.get(path).status_code == local_status, f"local mode: {path}"


class TestModeSurfaces:
    def test_foreign_host_is_accepted_in_server_mode(self):
        """The public API is proxied under real hostnames; the DNS-rebinding
        guard must not exist there."""
        with make_server_client(stub_rgc()) as client:
            response = client.get("/ping", headers={"host": "evil.example.com"})
        assert response.status_code == 200


class TestRootNamespaceIsReservedForTheSpa:
    """The JSON API lives at /v4 and /v1 -- never at the root, where a second
    copy of the catalog router would shadow the SPA's pages in both apps."""

    @pytest.mark.parametrize("client_name", ["server_client", "local_client"])
    def test_root_genomes_is_not_the_api(self, client_name, request):
        """GET /genomes is an SPA page, not a JSON listing, in both modes."""
        client = request.getfixturevalue(client_name)
        response = client.get("/genomes")
        assert response.headers["content-type"].startswith("text/html")

    @pytest.mark.parametrize("client_name", ["server_client", "local_client"])
    def test_unmatched_api_path_is_a_json_404(self, client_name, request):
        """An API typo must not come back as a 200 HTML document -- and the body
        is the standard envelope, because the SPA catch-all builds this response
        itself, bypassing the handlers."""
        client = request.getfixturevalue(client_name)
        response = client.get("/v4/does_not_exist")
        assert response.status_code == 404
        assert response.headers["content-type"].startswith("application/json")
        payload = response.json()
        assert payload["ok"] is False
        assert payload["error"]["code"] == "not_found"


class TestRouteHygiene:
    """Guards over the local app, matching the ones test_server.py runs."""

    def test_local_app_has_no_duplicate_routes(self, rgc):
        assert_no_duplicate_routes(make_local_app(rgc), "local app")

    def test_local_app_has_unique_operation_ids(self, rgc):
        assert_unique_operation_ids(make_local_app(rgc), "local app")

    @pytest.mark.parametrize("segment", SPA_CLIENT_ROUTES)
    def test_spa_client_routes_do_not_collide_with_the_api(self, rgc, segment):
        """No API route may shadow an SPA page, in either mode."""
        for app, label in ((make_server_app(rgc), "server"), (make_local_app(rgc), "local")):
            assert ("GET", f"/{segment}") not in set(route_keys(app)), f"{label}: /{segment}"


#: The frontend route table. Present in a source checkout; absent from a wheel.
_ROUTES_TSX = _REPO_ROOT / "frontend" / "src" / "app" / "routes.tsx"


@pytest.mark.skipif(not _ROUTES_TSX.is_file(), reason="frontend source not present")
def test_spa_route_mirror():
    """``SPA_CLIENT_ROUTES`` (Python) and the frontend route table must agree.
    Two dumb regex passes over ``routes.tsx``, on purpose -- no TS parser."""
    source = _ROUTES_TSX.read_text()

    const_block = re.search(r"export const SPA_CLIENT_ROUTES = \[(.*?)\]", source, re.DOTALL)
    assert const_block, "routes.tsx no longer exports SPA_CLIENT_ROUTES"
    ts_const = set(re.findall(r"'/([a-z0-9-]+)'", const_block.group(1)))
    assert ts_const == set(SPA_CLIENT_ROUTES), (
        "frontend/src/app/routes.tsx SPA_CLIENT_ROUTES drifted from "
        "refgenie/server/const.py:\n"
        f"  frontend only: {sorted(ts_const - set(SPA_CLIENT_ROUTES))}\n"
        f"  backend only:  {sorted(set(SPA_CLIENT_ROUTES) - ts_const)}"
    )

    registered = re.findall(r"path:\s*'([^']+)'", source)
    segments = {path.split("/")[0] for path in registered} - {"", "*"}
    assert segments == set(SPA_CLIENT_ROUTES), (
        "registered frontend routes drifted from SPA_CLIENT_ROUTES:\n"
        f"  registered only: {sorted(segments - set(SPA_CLIENT_ROUTES))}\n"
        f"  tuple only:      {sorted(set(SPA_CLIENT_ROUTES) - segments)}"
    )


# ===========================================================================
# SPA serving (refgenie.server.spa)
# ===========================================================================

#: Matches a hashed JS asset reference the way vite emits it.
_HASHED_JS_RE = re.compile(r"/_app/[A-Za-z0-9._-]+\.js")


@requires_web_assets
class TestSpaServing:
    """Everything here needs the frontend built into ``refgenie/server/webui/``."""

    @pytest.mark.parametrize("client_name", ["server_client", "local_client"])
    def test_spa_index_served(self, client_name, request):
        """Both modes serve the SPA shell at ``/`` -- there is no mode-specific root."""
        client = request.getfixturevalue(client_name)
        response = client.get("/")
        assert response.status_code == 200
        assert response.headers["content-type"].startswith("text/html")
        assert 'id="root"' in response.text

    def test_hashed_asset_served(self, server_client):
        """The hash is derived from the index, never hard-coded -- it rots on rebuild."""
        index = server_client.get("/").text
        match = _HASHED_JS_RE.search(index)
        assert match, f"no hashed JS asset referenced in index.html: {index!r}"
        response = server_client.get(match.group(0))
        assert response.status_code == 200
        assert response.content
        assert "javascript" in response.headers["content-type"]

    def test_asset_cache_headers(self, server_client):
        index_response = server_client.get("/")
        assert "no-cache" in index_response.headers["cache-control"]

        match = _HASHED_JS_RE.search(index_response.text)
        assert match
        asset_response = server_client.get(match.group(0))
        assert "immutable" in asset_response.headers["cache-control"]

    def test_client_route_falls_back_to_index(self, local_client):
        """Deep links are the first thing users hit and the first thing to break."""
        index = local_client.get("/").text
        deep_link = local_client.get("/genomes/abc123")
        assert deep_link.status_code == 200
        assert deep_link.text == index

    def test_unknown_asset_404s(self, server_client):
        """The registration-order regression test: a missing hashed asset must
        never fall back to index.html with a 200."""
        response = server_client.get("/_app/nope-deadbeef.js")
        assert response.status_code == 404
        assert "text/html" not in response.headers["content-type"]

    def test_service_info_bootstrap(self, server_client):
        """The SPA reads its bootstrap from /service-info (no injected global)."""
        data = server_client.get("/service-info").json()
        refgenie = data["refgenie"]
        assert refgenie["mode"] == "server"
        assert refgenie["api_base"] == "/v4"
        assert refgenie["root_path"] == ""
        assert "web_ui" in refgenie
        assert set(refgenie["capabilities"]) == set(CAPABILITY_KEYS)
        assert all(isinstance(v, bool) for v in refgenie["capabilities"].values())

    def test_root_path_rewrites_base_href(self, rgc):
        """One build, every deployment: root_path is a construction-time rewrite,
        not a rebuild. TestClient talks to the app post-reverse-proxy."""
        app = create_app(mode=APP_MODE_SERVER, root_path="/refgenie", refgenie_instance=rgc)
        with TestClient(app, raise_server_exceptions=True) as client:
            response = client.get("/")
            info = client.get("/service-info").json()
        assert '<base href="/refgenie/" />' in response.text
        assert info["refgenie"]["root_path"] == "/refgenie"

    def test_route_guard_single_spa_mount(self, rgc):
        """Exactly one catch-all, registered last -- Starlette matches routes in
        registration order, so anything after it would be dead."""
        app = create_app(mode=APP_MODE_SERVER, refgenie_instance=rgc)
        catch_alls = [
            r for r in app.router.routes if getattr(r, "path", None) == "/{full_path:path}"
        ]
        assert len(catch_alls) == 1
        assert app.router.routes[-1] is catch_alls[0]


def test_missing_assets_returns_503(rgc, monkeypatch):
    """Absence must be tolerated, never fatal. Not marked ``requires_web_assets``:
    this is the Node-free contributor path, exercised regardless of whether this
    checkout has a bundle. Patch the resolver itself so the "no bundle anywhere"
    case is exercised even on a machine that has built the frontend."""
    import refgenie.server.main as server_main

    monkeypatch.setattr(server_main, "resolve_web_dist", lambda explicit=None: None)

    app = create_app(mode=APP_MODE_LOCAL, refgenie_instance=rgc)
    with TestClient(app, base_url="http://localhost", raise_server_exceptions=True) as client:
        index_response = client.get("/")
        api_response = client.get("/v4/genomes")

    assert index_response.status_code == 503
    assert "npm --prefix frontend run build" in index_response.text
    assert api_response.status_code == 200


def test_skill_md_is_served_as_markdown(rgc, tmp_path):
    """A .md in the bundle root must come back as markdown, not as the SPA shell.

    ``spa_catch_all`` prefers a real file and lets ``mimetypes`` pick the type.
    The status code proves nothing here: a missing file falls through to
    index.html with a 200, so the content type is the assertion that matters.
    Builds its own tiny bundle, so this runs without ``npm run build``.
    """
    (tmp_path / "index.html").write_text("<!doctype html><base href='/' /><div id='root'></div>")
    (tmp_path / "SKILL.md").write_text("---\nname: refgenie-api\n---\n# x\n")

    app = create_app(mode=APP_MODE_SERVER, refgenie_instance=rgc, web_dist=tmp_path)
    with TestClient(app, raise_server_exceptions=True) as client:
        response = client.get("/SKILL.md")

    assert response.status_code == 200
    assert response.headers["content-type"].startswith("text/markdown")
    assert response.text.startswith("---")


@requires_web_assets
def test_built_bundle_carries_skill_md(server_client):
    """``frontend/public/*`` ships at the bundle root, so a real build serves the
    committed SKILL.md -- not the SPA shell with a 200.

    ``requires_web_assets`` proves a bundle exists, not that it is current, so a
    failure here is usually a stale bundle rather than a missing file.
    """
    response = server_client.get("/SKILL.md")
    hint = "Rebuild it: npm --prefix frontend run build"
    assert response.headers["content-type"].startswith("text/markdown"), (
        f"/SKILL.md came back as {response.headers['content-type']}, which means the "
        f"bundle in refgenie/server/webui/ has no SKILL.md and the SPA shell answered. {hint}"
    )
    assert "name: refgenie-api" in response.text, (
        f"The bundle's SKILL.md is not the one in frontend/public/. {hint}"
    )


# ===========================================================================
# OpenAPI <-> TypeScript wire-type drift guard (frontend/src/types/api.ts),
# and the same guard for the hand-written frontend/public/SKILL.md
# ===========================================================================

#: Property key sets, mirroring frontend/src/types/api.ts one interface at a time.
EXPECTED_PROPERTIES = {
    "GenomeResponse": {
        "digest",
        "aliases",
        "description",
        "asset_count",
        "species_name",
        "common_name",
        "taxon_id",
        "assembly_source",
        "assembly_accession",
    },
    "GenomeDetailResponse": {
        "digest",
        "description",
        "species_name",
        "common_name",
        "taxon_id",
        "assembly_source",
        "assembly_accession",
        "assembly_level",
        "remote_url",
        "taxon_uri",
        "fhr",
    },
    "AssetGroupPublic": {"id", "name", "description", "genome_digest", "asset_class_id"},
    "AssetResponse": {
        "digest",
        "name",
        "description",
        "recipe_id",
        "asset_group_id",
        "size",
        "serving_modes_override",
        "colocate",
        "serving_modes",
        "asset_class_name",
        "asset_group_name",
        "genome_digest",
        "names",
        "seek_keys",
        "is_default",
    },
    "SeekKeyResponse": {"name", "value", "description", "type", "size"},
    "AssetClassPublic": {"id", "name", "version", "description", "serving_modes"},
    "RecipePublic": {
        "id",
        "name",
        "version",
        "description",
        "output_asset_class_id",
        "command_templates",
        "input_params",
        "input_files",
        "input_assets",
        "docker_image",
        "custom_seek_keys",
        "default_asset",
        "inherent",
    },
    "StagedAssetPublic": {
        "asset_digest",
        "mode",
        "directory_contents",
        "build_commands",
        "download_count",
        "tarball_digest",
        "tarball_size",
    },
    "AliasPublic": {"name", "genome_digest"},
    "AliasResponse": {"alias", "digest", "source", "collection", "fhr"},
    "PaginationMeta": {"offset", "limit", "total"},
    "ArchiveRecord": {
        "digest",
        "asset_digest",
        "tarball_digest",
        "size",
        "directory_contents",
        "build_commands",
        "download_count",
    },
    "DatabaseSummaryResponse": {"genomes", "asset_groups", "assets"},
}

#: Endpoints the SPA's service layer calls on the shared read API.
EXPECTED_SHARED_PATHS = {
    "/genomes",
    "/genomes/{digest}",
    "/asset_groups",
    "/asset_groups/{id}",
    "/assets",
    "/assets/{digest}",
    "/assets/{asset_digest}/files",
    "/asset_classes",
    "/asset_classes/{id}",
    "/recipes",
    "/recipes/{id}",
    "/configurations",
    "/configurations/{id}",
    "/staged_assets",
    "/staged_assets/{id}",
    "/relationships/{asset_digest}",
    "/aliases",
    "/aliases/{name}",
}

#: Server-only reads the SPA gates on capabilities.archives / .downloads.
EXPECTED_VERSION4_PATHS = {
    "/archives",
    "/archives/{asset_digest}/download",
    "/assets/{asset_digest}/files/{file_path}",
    "/summary",
}


@pytest.fixture(scope="module")
def openapi_schema():
    # GA4GH DRS is included here (not mounted) in the real app, so it belongs in
    # the schema the SKILL.md guard below checks against. /seqcol and /mcp are
    # mounted sub-applications and never appear in any OpenAPI document.
    from refgenie.server.routers.ga4gh_drs import router as ga4gh_drs_router

    app = FastAPI()
    app.include_router(catalog_router.router, prefix="/v4")
    app.include_router(version4_router.router, prefix="/v4")
    app.include_router(ga4gh_drs_router, prefix="/v4/ga4gh/drs")
    return app.openapi()


@pytest.mark.parametrize("model_name", sorted(EXPECTED_PROPERTIES))
def test_wire_model_properties_match_frontend_types(openapi_schema, model_name):
    schemas = openapi_schema["components"]["schemas"]
    assert model_name in schemas, (
        f"{model_name} is no longer in the OpenAPI schema. "
        "frontend/src/types/api.ts declares it; update both."
    )
    actual = set(schemas[model_name].get("properties", {}))
    expected = EXPECTED_PROPERTIES[model_name]
    assert actual == expected, (
        f"{model_name} properties drifted from frontend/src/types/api.ts.\n"
        f"  added on the server:   {sorted(actual - expected)}\n"
        f"  missing on the server: {sorted(expected - actual)}\n"
        "Update frontend/src/types/api.ts (and whatever renders the field), "
        "then update EXPECTED_PROPERTIES here."
    )


@pytest.mark.parametrize(
    "expected,service_module",
    [
        (EXPECTED_SHARED_PATHS, "resources/*.ts"),
        (EXPECTED_VERSION4_PATHS, "resources/serverInfo.ts"),
    ],
    ids=["shared", "server_only"],
)
def test_paths_the_spa_calls_are_stable(openapi_schema, expected, service_module):
    paths = {p.removeprefix("/v4") for p in openapi_schema["paths"]}
    missing = expected - paths
    assert not missing, (
        f"Endpoints the SPA calls disappeared: {sorted(missing)}. "
        f"Update frontend/src/services/{service_module}."
    )


#: The agent-facing capability doc served at /SKILL.md on every refgenie origin.
#: It is hand-written, so it goes stale silently -- every /v4 path it names must
#: still exist in the schema. The doc's own convention: {braces} are path
#: parameters; <ANGLE_BRACKETS> are values the reader substitutes, and the regex
#: stops at the '<' so they are never mistaken for routes.
SKILL_MD = _REPO_ROOT / "frontend" / "public" / "SKILL.md"

_V4_PATH_RE = re.compile(r"/v4/[A-Za-z0-9_{}/.-]*")


def _normalize_api_path(path: str) -> str:
    return re.sub(r"\{[^}]*\}", "{}", path.rstrip("/.-"))


def _skill_md_v4_paths() -> set[str]:
    return {
        p
        for p in (_normalize_api_path(m) for m in _V4_PATH_RE.findall(SKILL_MD.read_text()))
        if p != "/v4"
    }


def test_skill_md_exists_and_is_a_skill_file():
    """A missing SKILL.md is invisible: /SKILL.md is not an API prefix, so the
    SPA catch-all answers 200 with index.html rather than 404. Assert presence
    positively -- no negative test will ever fire."""
    assert SKILL_MD.is_file(), f"{SKILL_MD} is missing; /SKILL.md would serve the SPA shell"
    text = SKILL_MD.read_text()
    assert text.startswith("---\n"), "SKILL.md must open with YAML frontmatter"
    head = text.split("---", 2)[1]
    assert "name:" in head and "description:" in head


def test_skill_md_paths_still_exist(openapi_schema):
    """Every /v4 path SKILL.md names still resolves to a real route.

    A worked example spells a parameter out (``/v4/aliases/hg38``) where a
    reference entry writes the template (``/v4/aliases/{name}``). Both are
    claims that the same route exists, so a concrete final segment is accepted
    when the one-parameter template above it does.
    """
    actual = {_normalize_api_path(p) for p in openapi_schema["paths"]}
    missing = sorted(
        p
        for p in _skill_md_v4_paths()
        if p not in actual and f"{p.rpartition('/')[0]}/{{}}" not in actual
    )
    assert not missing, (
        f"frontend/public/SKILL.md documents endpoints that no longer exist: "
        f"{missing}. Update SKILL.md."
    )


@pytest.mark.parametrize(
    "enum_name,members",
    [
        ("SeekKeyType", {"file", "directory", "prefix", "string", "json"}),
        ("SearchOperator", {"eq", "contains", "starts_with", "ends_with"}),
    ],
)
def test_enums_match_the_frontend_unions(openapi_schema, enum_name, members):
    assert set(openapi_schema["components"]["schemas"][enum_name]["enum"]) == members
