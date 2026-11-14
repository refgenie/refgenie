"""
RemoteGenomeSource protocol and implementations.

Every source reads from a :class:`~gtars.refget.RefgetStore`: a store URL is
opened directly, and a server URL is resolved to its backing store once,
through GA4GH service-info. Which genome an alias names is not a source's
question; ``ServerManager.locate`` asks the server's v4 API.
"""

import tempfile
from pathlib import Path
from typing import Protocol, runtime_checkable

from refgenie.logger import logger

_cache_root: tempfile.TemporaryDirectory | None = None


def _temp_cache_dir() -> Path:
    """A fresh store cache folder, inside one temp folder removed at interpreter exit."""
    global _cache_root
    if _cache_root is None:
        _cache_root = tempfile.TemporaryDirectory(prefix="refgenie_source_cache_")
    return Path(tempfile.mkdtemp(dir=_cache_root.name))


def normalize_server_url(url: str) -> str:
    """The one comparison form for server URLs: surrounding whitespace and trailing slashes removed."""
    return url.strip().rstrip("/")


@runtime_checkable
class RemoteGenomeSource(Protocol):
    """Common interface for remote genome metadata sources."""

    def verify_collection(self, digest: str) -> dict | None:
        """Verify a collection exists and return level2 data, or None."""
        ...

    def list_collections(self, page: int = 0, page_size: int = 100) -> dict:
        """List available collections with pagination."""
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
        self._store.set_quiet(True)
        # open_remote() does not pull the alias TSVs; a store with none is fine.
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

    def resolve_alias_in_namespace(self, alias: str, namespace: str) -> str | None:
        """Look up an alias in a specific namespace. Not part of the protocol."""
        meta = self._store.get_collection_metadata_by_alias(namespace, alias)
        return meta.digest if meta else None

    def get_collection_aliases(self, digest: str) -> list[tuple[str, str]]:
        """All (namespace, alias) pairs for a collection, from the aliases pulled at init."""
        return list(self._store.get_aliases_for_collection(digest))

    def get_collection_fhr(self, digest: str) -> dict | None:
        """The collection's FHR record as a camelCase dict, or None if absent.

        open_remote() does not fetch FHR sidecars, so this pulls them on demand.
        """
        try:
            self._store.pull_fhr(digest)
        except Exception as e:
            logger.debug(f"pull_fhr({digest}) failed for {self._url}: {e}")
        fhr = self._store.get_fhr_metadata(digest)
        return fhr.to_dict() if fhr is not None else None

    @property
    def store_url(self) -> str | None:
        return self._url

    @property
    def url(self) -> str:
        return self._url


# GA4GH service-info locations, tried in order: standalone seqcolapi or
# refgenie's root, then refgenie's seqcol sub-service.
SERVICE_INFO_PATHS = ("/service-info", "/seqcol/service-info")


def store_url_from_service_info(info: object) -> str | None:
    """Pull ``seqcol.refget_store.url`` out of a service-info document.

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

    Collection metadata comes from the store the server's service-info
    advertises first.

    Raises:
        ConnectionError: if no service-info responds, or it advertises no
            RefgetStore (genomes from it would have no sequences).
    """

    def __init__(self, server_url: str, cache_dir: Path | None = None):
        self._url = server_url.rstrip("/")

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

        logger.info(f"{self._url} resolves to RefgetStore {self._store_url}")
        self._store_source = RefgetStoreSource(self._store_url, cache_dir or _temp_cache_dir())

    def _fetch_service_info(self) -> dict | None:
        """Fetch the first service-info document that answers."""
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


# Normalized URL -> RemoteGenomeSource
_source_cache: dict[str, RemoteGenomeSource] = {}


def clear_source_cache() -> None:
    """Clear the cached sources. Useful for testing."""
    _source_cache.clear()


# Store-specific keys, so a stray {"version": ...} JSON blob does not match.
_RGSTORE_REQUIRED_KEYS = ("version", "mode", "seqdata_path_template")


def _looks_like_rgstore_manifest(data: object) -> bool:
    """Check whether decoded JSON looks like a valid rgstore.json manifest."""
    return isinstance(data, dict) and all(k in data for k in _RGSTORE_REQUIRED_KEYS)


def make_source(url: str, cache_dir: Path | None = None) -> RemoteGenomeSource:
    """Auto-detect the kind of remote and create a RemoteGenomeSource.

    A URL whose ``rgstore.json`` parses as a manifest is a
    :class:`RefgetStoreSource`; otherwise it is treated as a server and wrapped
    in a :class:`RefgenieServerSource`. Results are cached by normalized URL.

    Args:
        url: URL of the remote source.
        cache_dir: Cache directory for the RefgetStore. If omitted, a temporary
            one is used and removed at exit.

    Returns:
        A RemoteGenomeSource implementation.

    Raises:
        ConnectionError: if the URL is neither a store nor a store-backed server.
    """
    from refgenie.utils.http import make_client

    normalized = normalize_server_url(url)
    if normalized in _source_cache:
        return _source_cache[normalized]

    # Validate the body, not just the status: SPA hosts return 200 for any path.
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
                # A hint only; the JSON check below decides.
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
                logger.info(f"Detected RefgetStore at {url}")
                source: RemoteGenomeSource = RefgetStoreSource(
                    normalized, cache_dir or _temp_cache_dir()
                )
                _source_cache[normalized] = source
                return source
        else:
            logger.debug(f"{rgstore_url} returned {response.status_code}; not a store root.")

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


def collection_source(server_url: str, genome_attrs: dict) -> RemoteGenomeSource | None:
    """Where to read a genome's sequence collection.

    The store the server's genome record names (a server serves several stores
    but advertises only the first), else the store the server advertises.

    Args:
        server_url: The server that knows the genome.
        genome_attrs: The server's genome record (``ServerManager.locate``).

    Returns:
        The source, or None if neither URL is usable as one.
    """
    for url in (genome_attrs.get("refget_store_url"), server_url):
        if not url:
            continue
        try:
            return make_source(url)
        except Exception as exc:
            logger.debug(f"Could not use {url} as a genome source: {exc}")
    return None
