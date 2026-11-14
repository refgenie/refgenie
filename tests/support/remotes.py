"""Test doubles for remote sources, server clients and the puller."""

from contextlib import ExitStack, contextmanager
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

from tests.support.factories import GENOME, OMIT


def mock_server_client(
    *,
    server_url="http://test.example.com",
    asset_group_name="bowtie2_index",
    genome_digest="genome_digest_001",
    asset_group_id=1,
    asset_name="default",
    asset_digest="abc123digest",
    serving_modes=("archive",),
    seek_keys=OMIT,
    is_default=OMIT,
    alias_digest=OMIT,
    staged_items=(),
    file_list=(),
    relations=None,
    aliases=OMIT,
    download_side_effect=None,
) -> MagicMock:
    """A MagicMock RefgenieserverClient for the pull/seekr dispatch tests.

    Any argument left at :data:`OMIT` omits that key from the payload (or, for
    ``aliases``, leaves the attribute *unset* on the mock -- an unset attribute
    returns a truthy MagicMock, which is not the same as ``[]``).

    ``seek_keys`` populates the asset payload's ``seek_keys`` (seekr reads the
    default/named seek key straight off it). ``alias_digest`` makes
    ``client.alias.resolve`` return that digest,
    for remote resolution of a non-local genome.
    """
    client = MagicMock()
    client.server_url = server_url

    group = {"id": asset_group_id, "name": asset_group_name}
    if genome_digest is not OMIT:
        group["genome_digest"] = genome_digest
    client.asset_group.list.return_value = [group]

    asset = {
        "id": 1,
        "name": asset_name,
        "digest": asset_digest,
        "names": [{"name": asset_name}],
    }
    if asset_group_id is not OMIT:
        asset["asset_group_id"] = asset_group_id
    if serving_modes is not OMIT:
        asset["serving_modes"] = list(serving_modes)
    if seek_keys is not OMIT:
        asset["seek_keys"] = [dict(sk) for sk in seek_keys]
    if is_default is not OMIT:
        asset["is_default"] = is_default
    client.asset.list.return_value = [asset]

    client.staged_asset.list.return_value = list(staged_items)
    client.asset.files.return_value = list(file_list)
    _relations_payload = (
        dict(relations) if relations is not None else {"parents": [], "children": []}
    )
    client.asset.relationships.return_value = _relations_payload
    if alias_digest is not OMIT:
        client.alias.resolve.return_value = alias_digest
    if aliases is not OMIT:
        client.alias.list.return_value = list(aliases)
    if download_side_effect is not None:
        client.download_with_progress.side_effect = download_side_effect
    return client


def real_client_no_init(url: str = "http://test-server"):
    """A REAL RefgenieserverClient built via ``__new__`` (no network in __init__).

    Deliberately separate from :func:`mock_server_client`: the point is that
    ``get_paginated`` and friends run their real code.
    """
    from refgenie.managers.sources.client import RefgenieserverClient

    client = RefgenieserverClient.__new__(RefgenieserverClient)
    client.server_url = url
    client.openapi_endpoint = "/openapi.json"
    client.api_version = "v4"
    client._http_client = None
    return client


class MockRemoteSource:
    """A faithful in-memory RemoteGenomeSource for testing."""

    def __init__(
        self,
        collections=None,
        store_url="https://store.example.com",
        collection_aliases=None,
        fhr=None,
    ):
        self._collections = collections or {}
        self._store_url = store_url
        # digest -> [(namespace, alias), ...]
        self._collection_aliases = collection_aliases or {}
        # digest -> camelCase FHR dict
        self._fhr = fhr or {}

    def verify_collection(self, digest):
        return self._collections.get(digest)

    def list_collections(self, page=0, page_size=100):
        return {"results": [{"digest": d} for d in self._collections]}

    def get_collection_aliases(self, digest):
        return list(self._collection_aliases.get(digest, []))

    def get_collection_fhr(self, digest):
        return self._fhr.get(digest)

    @property
    def store_url(self):
        return self._store_url

    @property
    def url(self):
        return "https://mock.example.com"


def assert_mock_remote_source_conforms() -> None:
    """The test double must actually satisfy the protocol, or every test using
    it proves nothing about real sources."""
    from refgenie.managers.sources.genomes import RemoteGenomeSource

    assert isinstance(MockRemoteSource(), RemoteGenomeSource)


@contextmanager
def mocked_puller(
    rg,
    mock_client,
    *,
    alias=GENOME,
    subscriptions=("http://test.example.com",),
    mock_genome=True,
    mock_download_modes=True,
    mock_asset_writes=True,
):
    """Point ``rg``'s AssetPuller at ``mock_client`` with the usual patches.

    Yields a namespace with ``puller``, ``client`` and (when patched) the
    ``archive``, ``file`` and ``ensure`` mocks so assertions still work.

    * ``mock_genome=False`` leaves ``genome.ensure_from_remote`` (and the asset
      ``exists`` check) real -- required by the rollback tests, which are
      *about* genome creation.
    * ``mock_download_modes=False`` lets the real ``_pull_*_mode`` run (and
      raise).
    * ``mock_asset_writes=False`` leaves the cross-manager writes real.
    """
    puller = rg.transfer._puller
    ns = SimpleNamespace(puller=puller, client=mock_client, archive=None, file=None, ensure=None)

    with ExitStack() as stack:
        stack.enter_context(
            patch.object(puller._servers, "subscriptions", return_value=list(subscriptions))
        )
        stack.enter_context(patch.object(puller._servers, "client", return_value=mock_client))
        if mock_genome:
            ensure = stack.enter_context(patch.object(puller._genome_manager, "ensure_from_remote"))
            ensure.return_value = MagicMock(
                genome_digest=rg.alias.resolve(alias),
                created_genome=False,
                created_alias=False,
                alias_name=alias,
            )
            ns.ensure = ensure
            stack.enter_context(patch.object(puller._asset, "exists", return_value=False))
        # The mock client serves no asset classes; treat the asset's class as
        # registered so these tests reach the download they are about.
        stack.enter_context(
            patch.object(puller._asset_class, "exists_any_version", return_value=True)
        )
        if mock_download_modes:
            ns.archive = stack.enter_context(patch.object(puller, "_pull_archive_mode"))
            ns.file = stack.enter_context(patch.object(puller, "_pull_file_mode"))
        if mock_asset_writes:
            stack.enter_context(
                patch.object(puller._asset.content, "add", return_value=MagicMock())
            )
            stack.enter_context(patch.object(puller._asset.group, "get_default", return_value=None))
            stack.enter_context(patch.object(puller._asset.group, "set_default"))
        yield ns
