"""
RemoteGenomeSource protocol and implementations.

Provides a common interface for remote genome metadata sources. Every source is
ultimately a :class:`~gtars.refget.RefgetStore`: a store URL is opened directly,
and a *server* URL is resolved to its backing store exactly once, through a
GA4GH service-info fetch, and then never spoken to over seqcol REST again.

The one thing a store cannot answer is flat (non-namespaced) alias resolution,
which is not a seqcol capability at all -- refgenie's aliases are flat SQL rows,
while seqcol's alias endpoints are namespaced. That goes to refgenie's own v4
API via the OpenAPI/operationId-driven RefgenieserverClient, never through a
seqcol path.
"""

from pathlib import Path
from typing import Protocol, runtime_checkable

from refgenie.logger import logger
from refgenie.managers.sources.api_ids import API_ID_ALIAS_DIGEST


@runtime_checkable
class RemoteGenomeSource(Protocol):
    """Common interface for remote genome metadata sources."""

    def verify_collection(self, digest: str) -> dict | None:
        """Verify a collection exists and return level2 data, or None."""
        ...

    def list_collections(self, page: int = 0, page_size: int = 100) -> dict:
        """List available collections with pagination."""
        ...

    def resolve_alias(self, alias: str) -> str | None:
        """Resolve a human-friendly name to a digest. Returns None if not found."""
        ...

    def get_collection_aliases(self, digest: str) -> list[tuple[str, str]]:
        """All (namespace, alias) pairs the source holds for a collection."""
        ...

    def get_collection_fhr(self, digest: str) -> dict | None:
        """The collection's FHR record as a camelCase dict, or None if absent."""
        ...

    @property
    def store_url(self) -> str | None:
        """The RefgetStore URL for sequence fetching. Stored on Genome.remote_url."""
        ...

    @property
    def url(self) -> str:
        """The primary URL of this source (for display/logging)."""
        ...


class RefgetStoreSource:
    """Backend using RefgetStore.open_remote() (static files on S3)."""

    def __init__(self, store_url: str, cache_dir: Path):
        from gtars.refget import RefgetStore

        self._url = store_url
        self._store = RefgetStore.open_remote(cache_dir, store_url)
        # open_remote() does not pull the alias TSVs, so namespace/alias lookups
        # return empty (breaking `genome init --store --namespace`). Pull them
        # here; a store with no aliases is tolerated.
        try:
            self._store.pull_aliases()
        except Exception as e:
            logger.debug(f"pull_aliases failed for {store_url}: {e}")

    def verify_collection(self, digest: str) -> dict | None:
        meta = self._store.get_collection_metadata(digest)
        if meta is None:
            return None
        return self._store.get_collection_level2(digest)

    def list_collections(self, page: int = 0, page_size: int = 100) -> dict:
        result = self._store.list_collections(page=page, page_size=page_size)
        return result

    def resolve_alias(self, alias: str) -> str | None:
        return None  # Raw stores don't support alias resolution; use a server

    def resolve_alias_in_namespace(self, alias: str, namespace: str) -> str | None:
        """Look up an alias in a specific namespace. Not part of the protocol."""
        meta = self._store.get_collection_metadata_by_alias(namespace, alias)
        return meta.digest if meta else None

    def get_collection_aliases(self, digest: str) -> list[tuple[str, str]]:
        """All (namespace, alias) pairs for a collection, from the alias sidecars
        pulled at init. Answered from the already-loaded indexes -- no round trip."""
        return list(self._store.get_aliases_for_collection(digest))

    def get_collection_fhr(self, digest: str) -> dict | None:
        """The collection's FHR record as a camelCase dict, or None if absent.

        open_remote() does not fetch fhr/ sidecars, so pull this collection's on
        demand; a store (or collection) with no FHR is tolerated.
        """
        try:
            self._store.pull_fhr(digest)
        except Exception as e:
            logger.debug(f"pull_fhr({digest}) failed for {self._url}: {e}")
        fhr = self._store.get_fhr_metadata(digest)
        return fhr.to_dict() if fhr is not None else None

    @property
    def store_url(self) -> str | None:
        return self._url  # The store URL IS this source's URL

    @property
    def url(self) -> str:
        return self._url


# Well-known GA4GH discovery documents, tried in order. `/service-info` is the
# standalone seqcolapi location and refgenie's root document; `/seqcol/service-info`
# is refgenie's seqcol sub-service. This is a bounded bootstrap, not prefix
# discovery: nothing derives per-endpoint paths from which one answered.
SERVICE_INFO_PATHS = ("/service-info", "/seqcol/service-info")


def store_url_from_service_info(info: object) -> str | None:
    """Pull ``seqcol.refget_store.url`` out of a service-info document.

    Both refgenie's root document and any store-backed seqcol service-info use
    this nesting, so one reader covers every location in SERVICE_INFO_PATHS.
    Returns None when the document does not advertise an enabled store.
    """
    if not isinstance(info, dict):
        return None
    seqcol = info.get("seqcol")
    if not isinstance(seqcol, dict):
        return None
    store = seqcol.get("refget_store")
    if not isinstance(store, dict) or not store.get("enabled"):
        return None
    return store.get("url")


class RefgenieServerSource:
    """A refgenie / seqcolapi server URL, resolved once to its backing RefgetStore.

    Collection metadata comes from the store. Flat alias resolution goes to the
    server's own v4 API through RefgenieserverClient, which discovers the alias
    path from the server's published OpenAPI document -- never by prepending a
    seqcol prefix to it.

    Raises:
        ConnectionError: if no service-info responds, or if the one that does
            advertises no RefgetStore. A server with no store cannot supply
            sequences, so registering genomes from it would silently produce
            sequence-less genomes.
    """

    def __init__(self, server_url: str, cache_dir: Path | None = None):
        self._url = server_url.rstrip("/")
        self._client = None
        self._client_failed = False

        info = self._fetch_service_info()
        if info is None:
            raise ConnectionError(
                f"No service-info found at {self._url}. Tried: "
                + ", ".join(self._url + p for p in SERVICE_INFO_PATHS)
            )
        self._store_url = store_url_from_service_info(info)
        if self._store_url is None:
            raise ConnectionError(
                f"{self._url} does not advertise a RefgetStore "
                f"(service-info has no enabled seqcol.refget_store.url), so it "
                f"cannot be used as a genome source. A refgenie server needs at "
                f"least one enabled store in its registry to serve as one."
            )

        if cache_dir is None:
            import tempfile

            cache_dir = Path(tempfile.mkdtemp(prefix="refgenie_source_cache_"))
        logger.info(f"{self._url} resolves to RefgetStore {self._store_url}")
        self._store_source = RefgetStoreSource(self._store_url, cache_dir)

    def _fetch_service_info(self) -> dict | None:
        """Fetch the first service-info document that answers. One-time bootstrap."""
        from refgenie.utils.http import make_client

        with make_client(timeout=5.0) as client:
            for path in SERVICE_INFO_PATHS:
                url = self._url + path
                try:
                    response = client.get(url)
                except Exception as e:
                    logger.warning(f"service-info request to {url} failed: {e}")
                    continue
                if response.status_code != 200:
                    logger.debug(f"service-info at {url} returned {response.status_code}")
                    continue
                try:
                    return response.json()
                except Exception as e:
                    logger.warning(f"service-info at {url} returned unparseable JSON: {e}")
        return None

    def verify_collection(self, digest: str) -> dict | None:
        return self._store_source.verify_collection(digest)

    def list_collections(self, page: int = 0, page_size: int = 100) -> dict:
        return self._store_source.list_collections(page=page, page_size=page_size)

    def resolve_alias(self, alias: str) -> str | None:
        """Resolve a flat alias through the server's v4 alias endpoint.

        The path comes from the server's OpenAPI document (by operationId);
        prepending a seqcol prefix here silently breaks alias resolution.
        """
        from refgenie.managers.sources.client import RefgenieserverClient

        if self._client_failed:
            return None
        if self._client is None:
            try:
                self._client = RefgenieserverClient(server_url=self._url)
            except Exception as e:
                self._client_failed = True
                logger.warning(f"Could not read the OpenAPI spec from {self._url}: {e}")
                return None

        if not self._client.has_endpoint(API_ID_ALIAS_DIGEST):
            logger.warning(
                f"{self._url} publishes no '{API_ID_ALIAS_DIGEST}' operation; "
                f"cannot resolve alias '{alias}' there."
            )
            return None
        try:
            data = self._client.get(
                operation_id=API_ID_ALIAS_DIGEST,
                url_format_params={"name": alias},
            )
        except Exception as e:
            logger.warning(f"Alias '{alias}' could not be resolved at {self._url}: {e}")
            return None
        if isinstance(data, dict):
            return data.get("digest")
        logger.warning(f"Alias endpoint at {self._url} returned an unexpected payload: {data!r}")
        return None

    def get_collection_aliases(self, digest: str) -> list[tuple[str, str]]:
        return self._store_source.get_collection_aliases(digest)

    def get_collection_fhr(self, digest: str) -> dict | None:
        return self._store_source.get_collection_fhr(digest)

    @property
    def store_url(self) -> str | None:
        return self._store_url

    @property
    def url(self) -> str:
        return self._url


# Module-level cache: normalized URL -> RemoteGenomeSource
_source_cache: dict[str, RemoteGenomeSource] = {}


def clear_source_cache() -> None:
    """Clear the cached sources. Useful for testing."""
    _source_cache.clear()


# Fields every rgstore.json manifest carries. Require the generic `version`
# plus store-specific keys so a stray {"version": ...} JSON blob can't match.
_RGSTORE_REQUIRED_KEYS = ("version", "mode", "seqdata_path_template")


def _looks_like_rgstore_manifest(data: object) -> bool:
    """Check whether decoded JSON looks like a valid rgstore.json manifest."""
    return isinstance(data, dict) and all(k in data for k in _RGSTORE_REQUIRED_KEYS)


def make_source(url: str, cache_dir: Path | None = None) -> RemoteGenomeSource:
    """Auto-detect the kind of remote and create a RemoteGenomeSource.

    Two kinds exist, and both end up reading from a RefgetStore:

    1. The URL is a store root -- ``{url}/rgstore.json`` parses as a manifest.
       Opened directly as a :class:`RefgetStoreSource`.
    2. The URL is a refgenie / seqcolapi server -- a service-info document
       names its backing store. Wrapped in a :class:`RefgenieServerSource`.

    Results are cached by normalized URL, so detection happens at most once per
    URL per process.

    Args:
        url: URL of the remote source.
        cache_dir: Cache directory for the RefgetStore. A temporary directory is
            created if omitted.

    Returns:
        A RemoteGenomeSource implementation.

    Raises:
        ConnectionError: if the URL is neither a store nor a store-backed server.
    """
    from refgenie.utils.http import make_client

    normalized = url.rstrip("/")
    if normalized in _source_cache:
        return _source_cache[normalized]

    # Try RefgetStore: fetch {url}/rgstore.json and validate the body is a
    # real store manifest. A bare 200 status is not enough -- SPA hosts
    # (e.g. serving index.html for every path) return 200 for any URL.
    rgstore_url = normalized + "/rgstore.json"
    try:
        with make_client(timeout=5.0) as client:
            response = client.get(rgstore_url)
    except Exception as e:
        logger.warning(f"Store probe of {rgstore_url} failed: {e}")
        response = None

    if response is not None:
        if response.status_code == 200:
            content_type = response.headers.get("content-type", "")
            if "text/html" in content_type:
                # Hint only -- some static hosts (S3, etc.) serve rgstore.json
                # with a generic content-type, so this doesn't disqualify it
                # outright, but an HTML content-type is a strong signal of an
                # SPA fallback page. Log and let the JSON-parse check decide.
                logger.debug(
                    f"{rgstore_url} responded with text/html content-type; "
                    "likely not a real rgstore.json (SPA fallback?)."
                )
            try:
                data = response.json()
            except Exception as e:
                logger.debug(f"{rgstore_url} did not return JSON ({e}); not a store root.")
                data = None
            if _looks_like_rgstore_manifest(data):
                if cache_dir is None:
                    import tempfile

                    cache_dir = Path(tempfile.mkdtemp())
                logger.info(f"Detected RefgetStore at {url}")
                source: RemoteGenomeSource = RefgetStoreSource(normalized, cache_dir)
                _source_cache[normalized] = source
                return source
        else:
            logger.debug(f"{rgstore_url} returned {response.status_code}; not a store root.")

    # Not a store root -- treat it as a server and bootstrap to its store.
    logger.info(f"{url} is not a RefgetStore root; looking for a server service-info.")
    try:
        source = RefgenieServerSource(normalized, cache_dir)
    except ConnectionError as e:
        raise ConnectionError(
            f"Could not use {url} as a genome source. Probed {rgstore_url} for a "
            f"RefgetStore manifest and "
            + ", ".join(normalized + p for p in SERVICE_INFO_PATHS)
            + f" for a service-info. {e}"
        ) from e
    _source_cache[normalized] = source
    return source
