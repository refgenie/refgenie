"""Local-action stubs, preflight requests, job results and the phase vocabulary."""

from unittest.mock import MagicMock

from tests.support.server import stub_rgc

#: Every state-changing web route carries this anti-CSRF header.
ACTION_HEADERS = {"x-refgenie-action": "1"}


#: The bridge's default allowlisted origin, and one that is never allowed.
PUBLIC_ORIGIN = "https://refgenie.org"


EVIL_ORIGIN = "https://evil.example"


#: Shared capability vocabulary -- the same seventeen keys /service-info emits.
CAPABILITY_KEY_SET = {
    "pull",
    "build",
    "delete",
    "aliases_write",
    "subscriptions",
    "data_channels",
    "recipes_write",
    "asset_classes_write",
    "remote_browse",
    "genome_init",
    "jobs",
    "jobs_cancel",
    "downloads",
    "archives",
    "seqcol",
    "drs",
    "local_files",
}


def act(client, method, path, json=None, headers=None):
    """Issue a state-changing request with the action header attached."""
    return client.request(method, path, json=json, headers={**ACTION_HEADERS, **(headers or {})})


class FakeRefgenie:
    """Just enough of a Refgenie for the jobs router's preflight/stage check."""

    def __init__(self, genome_stage_folder=None):
        self.genome_stage_folder = genome_stage_folder


def job_result_ok():
    """The typed ``JobResult`` every fake runner returns."""
    from refgenie.server.jobs.schemas import JobResult

    return JobResult(asset_digest="d" * 64, registry_path="rCRSd/fasta:default")


def instant_runner(ctx):
    """A job runner that succeeds immediately."""
    return job_result_ok()


def _stub_find_subscription(rgc, url):
    """``ServerManager.find_subscription`` over a stub's ``subscriptions`` return value."""
    from refgenie.managers.sources import normalize_server_url

    wanted = normalize_server_url(url)
    for stored in rgc.servers.subscriptions():
        if normalize_server_url(stored) == wanted:
            return stored
    return None


def web_stub_rgc():
    """A MagicMock Refgenie safe to drive every actions route end to end.

    Return values are plain strings/lists so job results and ActionResult
    payloads serialize; the mocked asset satisfies the runners' ``_asset_result``.
    """
    rgc = stub_rgc()
    asset = MagicMock()
    asset.digest = "assetdigest123"
    asset.name = "test"
    asset.registry_path = "rCRSd/fasta:test"
    asset.asset_group.name = "fasta"
    asset.asset_group.genome.digest = "genomedigest123"
    rgc.transfer.pull.return_value = asset
    rgc.build.run.return_value = asset
    rgc.build.initialize_and_build.return_value = ("genomedigest123", True)
    rgc.asset.remove_by_digest.return_value = "rCRSd/fasta:test"
    rgc.servers.subscriptions.return_value = []
    rgc.servers.find_subscription.side_effect = lambda url: _stub_find_subscription(rgc, url)
    rgc.build.preflight.return_value = {"ok": True, "errors": [], "resolved": {}}
    return rgc


def preflight(client, path, origin, method="GET", **extra_headers):
    """Issue a CORS preflight (``OPTIONS``) for ``origin``."""
    headers = {"Origin": origin, "Access-Control-Request-Method": method, **extra_headers}
    return client.options(path, headers=headers)


#: The phase names the web UI knows how to render, per job kind.
#:
#: Transcribed from ``frontend/src/components/jobs/phases.ts``, which turns each
#: into a label and a "step 4 of 8" counter, and gives three of them explicit
#: "this is silent and slow" copy. A phase the backend emits that is not here
#: still renders (humanized), but loses the step counter and the reassurance --
#: so keep the two lists in sync.
PHASE_VOCABULARY = {
    "pull": {
        "resolve",
        "query",
        "stage_lookup",
        "download",
        "verify",
        "unpack",
        "register",
        "symlink",
    },
    "build": {
        "resolve_recipe",
        "seek_keys",
        "validate",
        "run",
        "digest",
        "register",
        "stage",
    },
    "genome_init": {"resolve", "digest", "register", "build"},
}
