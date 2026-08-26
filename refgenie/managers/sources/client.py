import json
import logging
import os
import time
from functools import cached_property
from pathlib import Path
from typing import Protocol, runtime_checkable
from urllib.parse import urljoin

import httpx
from rich.progress import Progress

from refgenie import progress
from refgenie.config import config
from refgenie.const import API_VERSION, DEFAULT_PAGE_SIZE, MAX_PAGE_SIZE
from refgenie.logger import logger
from refgenie.utils.console import RICH_PROGRESS_COLUMNS
from refgenie.utils.http import make_client
from refgenie.managers.sources.api_ids import (
    API_ID_ASSET_FILE_DOWNLOAD,
    API_ID_ASSET_FILES,
    API_ID_STAGED_ASSETS,
)
from refgenie.managers.sources.openapi import HttpMethod, RefgenieserverOpenApiSpec

logging.getLogger("httpx").setLevel(config.log_level.value)


@runtime_checkable
class ServerClient(Protocol):
    """
    Protocol for the server client.
    """

    @property
    def server_url(self) -> str:
        """
        The URL of the server the client is connected to.
        """
        ...

    def get(
        self,
        operation_id: str,
        params: dict | None = None,
        url_format_params: dict[str, str | None] = None,
    ) -> dict:
        """
        Send a GET request to the specified operation ID.
        """
        ...

    def get_asset_groups(self, params: dict | None = None) -> list[dict]:
        """Get one page of asset groups, unwrapping the paginated response envelope."""
        ...

    def get_assets(self, params: dict | None = None) -> list[dict]:
        """Get one page of assets, unwrapping the paginated response envelope."""
        ...

    def get_staged_assets(self, params: dict | None = None) -> list[dict]:
        """Get one page of staged assets, unwrapping the paginated response envelope."""
        ...

    def download_with_progress(
        self,
        operation_id: str,
        output_path: Path,
        params: dict | None = None,
        url_format_params: dict[str, str | None] = None,
        name: str | None = None,
    ) -> Path:
        """
        Download the asset served by the given operation_id to output_path,
        showing progress along the way.

        Implementations should report byte progress to `refgenie.progress`
        whenever `refgenie.progress.active_sink()` is not None, and should not
        build a `rich` live display in that case -- something other than a
        terminal is consuming the progress.
        """
        ...

    def get_all_assets(self, params: dict | None = None) -> list[dict]:
        """Get all assets using pagination."""
        ...

    def get_all_asset_groups(self, params: dict | None = None) -> list[dict]:
        """Get all asset groups using pagination."""
        ...

    def get_all_genomes(self, params: dict | None = None) -> list[dict]:
        """Get all genomes using pagination."""
        ...

    def get_all_aliases(self, params: dict | None = None) -> list[dict]:
        """Get all aliases using pagination."""
        ...

    def get_all_staged_assets(self, params: dict | None = None) -> list[dict]:
        """Get all staged assets using pagination."""
        ...

    def get_asset_file_list(self, asset_digest: str) -> list[str]:
        """Get the list of files available for a file-mode asset."""
        ...

    def download_file(self, asset_digest: str, file_path: str, output_path: Path) -> Path:
        """Download a single file from a file-mode asset."""
        ...


#: Report to a progress sink at most this often, in seconds...
_SINK_MIN_INTERVAL_SECONDS = 0.2
#: ...or every this many bytes, whichever comes first.
_SINK_MIN_BYTES = 1024 * 1024

#: Timeout for archive/file downloads. httpx's own default read timeout is 5s,
#: which a genome-scale archive chunk can easily exceed under load; the rest of
#: this file uses a flat 10s for small metadata calls (`_send_request`), but a
#: multi-GB stream needs a longer read/write allowance. `pool=None` disables
#: the connection-pool-checkout timeout -- irrelevant for a one-off stream.
_DOWNLOAD_TIMEOUT = httpx.Timeout(connect=10, read=60, write=60, pool=None)


class _RichDownloadReporter:
    """Drive a `rich` progress bar. The terminal path; unchanged behavior."""

    def __init__(self, rich_progress: Progress, label: str | None):
        self._progress = rich_progress
        self._label = label
        self._task = None

    def start(self, total: int | None) -> None:
        self._task = self._progress.add_task("Download", n=self._label, total=total)

    def update(self, downloaded: int) -> None:
        self._progress.update(self._task, completed=downloaded)

    def finish(self, downloaded: int) -> None:
        self._progress.update(self._task, completed=downloaded)


class _SinkDownloadReporter:
    """Report byte counts to the installed `refgenie.progress` sink.

    Throttled: at most one event per `_SINK_MIN_INTERVAL_SECONDS` or per
    `_SINK_MIN_BYTES`, whichever comes first, plus one final event. A 40 GB
    archive downloaded in 64 KB chunks is ~650,000 chunks; an unthrottled
    event per chunk would drown the job manager's ring buffer in a second.

    `total` is None when the server sent no Content-Length -- the UI renders
    an indeterminate bar in that case rather than dividing by zero.
    """

    def __init__(self, label: str | None):
        self._label = label
        self._total: int | None = None
        self._last_time = 0.0
        self._last_bytes = 0

    def start(self, total: int | None) -> None:
        self._total = total
        self._last_time = time.monotonic()
        self._last_bytes = 0
        progress.emit(
            "progress",
            message=self._label,
            current=0,
            total=total,
            unit="bytes",
            phase="download",
        )

    def update(self, downloaded: int) -> None:
        now = time.monotonic()
        if (
            now - self._last_time < _SINK_MIN_INTERVAL_SECONDS
            and downloaded - self._last_bytes < _SINK_MIN_BYTES
        ):
            return
        self._last_time = now
        self._last_bytes = downloaded
        progress.emit(
            "progress",
            message=self._label,
            current=downloaded,
            total=self._total,
            unit="bytes",
            phase="download",
        )

    def finish(self, downloaded: int) -> None:
        progress.emit(
            "progress",
            message=self._label,
            current=downloaded,
            total=self._total if self._total is not None else downloaded,
            unit="bytes",
            phase="download",
        )


class RefgenieserverClient:
    def __init__(
        self,
        server_url: str,
        api_version: str = API_VERSION,
        openapi_endpoint: str = "/openapi.json",
        http_client: httpx.Client | None = None,
    ):
        self.server_url = server_url
        self.openapi_endpoint = openapi_endpoint
        self.api_version = api_version
        self._http_client = http_client  # For testing with custom transport

        logger.debug(f"Connecting to the server at {self.server_url}")

        # Fetch the spec eagerly so a connection failure surfaces here, not later
        # from an unusable client with an unrelated message.
        logger.info(f"Connected to server: {self.server_url} {self.openapi_spec.info}")

    def __repr__(self):
        return f"RefgenieserverClient(server_url={self.server_url})"

    @cached_property
    def openapi_spec(self) -> RefgenieserverOpenApiSpec:
        """
        Validate the openAPI JSON description of the server. Memoized to avoid
        multiple requests to the server.

        Returns:
            RefgenieserverOpenApiSpec: openAPI JSON description
        """
        logger.debug(f"Retrieving OpenAPI spec from {self.openapi_endpoint}")
        return RefgenieserverOpenApiSpec.model_validate(self._get(endpoint=self.openapi_endpoint))

    @cached_property
    def endpoints_mapping(self) -> dict:
        """
        Map of the server endpoints by their operationIds. Memoized to avoid
        unnecessary parsing of the openAPI spec.

        Returns:
            dict: endpoints mapped by their operationIds
        """
        mapping = {}
        for path, path_item in self.openapi_spec.paths.items():
            for operation in [
                path_item.get,
                path_item.post,
                path_item.put,
                path_item.delete,
                path_item.head,
            ]:
                if operation and operation.operationId:
                    mapping[operation.operationId] = path
        return mapping

    def has_endpoint(self, operation_id: str) -> bool:
        """Check if the server supports a given operation ID."""
        return operation_id in self.endpoints_mapping

    def _send_request(
        self,
        endpoint: str,
        params: dict | None = None,
        method: HttpMethod = HttpMethod.GET,
    ) -> httpx.Response:
        """
        Send a request to the specified endpoint.

        Args:
            endpoint (str): server endpoint
            params (dict, optional): query parameters
            method (HttpMethod): HTTP method

        Returns:
            httpx.Response: server response
        """
        url = urljoin(self.server_url, endpoint)
        if self._http_client:
            response = self._http_client.request(
                method=method.value,
                url=url,
                params=params,
            )
        else:
            with make_client() as client:
                response = client.request(
                    method=method.value,
                    url=url,
                    params=params,
                )
        response.raise_for_status()
        logger.debug(response.json())
        return response

    def _get(self, endpoint: str, params: dict | None = None) -> dict | str:
        """
        Send a GET request to the specified endpoint.

        Args:
            endpoint (str): server endpoint
            params (dict, optional): query parameters

        Returns:
            dict or str: served data
        """
        response = self._send_request(endpoint=endpoint, params=params, method=HttpMethod.GET)
        try:
            return response.json()
        except json.JSONDecodeError:
            return response.text

    @staticmethod
    def extract_items_from_response(response: dict | list | str) -> list:
        """
        Extract items from a potentially paginated response.

        Args:
            response: The response from the server, either a list or a dict
                     with pagination structure {"items": [...], "pagination": {...}}
                     or a string (which we'll treat as empty)

        Returns:
            List: The items from the response
        """
        if isinstance(response, dict) and "items" in response:
            # Paginated response format
            return response["items"]
        elif isinstance(response, list):
            # Direct list format (non-paginated)
            return response
        else:
            # Fallback for other formats (including strings)
            return []

    def get(
        self,
        operation_id: str,
        params: dict | None = None,
        url_format_params: dict[str, str | None] = None,
    ) -> dict | str:
        """
        Send a GET request to the specified operation ID.

        Args:
            operation_id (str): operation ID
            params (dict, optional): query parameters

        Returns:
            dict or str: served data
        """
        if operation_id not in self.endpoints_mapping:
            raise ValueError(
                f"Operation ID '{operation_id}' is not present in the server's endpoints. "
                f"Available operation IDs: {self.endpoints_mapping.keys()}"
            )
        url_template = self.endpoints_mapping[operation_id]
        url = url_template.format(**url_format_params) if url_format_params else url_template
        return self._get(url, params)

    def download_with_progress(
        self,
        operation_id: str,
        output_path: Path,
        params: dict | None = None,
        url_format_params: dict[str, str | None] = None,
        name: str | None = None,
    ) -> Path:
        """
        Download the asset served by the given operation_id to output_path,
        showing progress along the way.

        Progress goes to a `rich` bar on the terminal, unless a
        `refgenie.progress` sink is installed for this context -- then it is
        reported to the sink instead and no `rich` display is built. See the
        comment on the reporting branch below.

        Args:
            operation_id: operationId of the download endpoint.
            output_path: where to save the file.
            params: query parameters.
            url_format_params: values for URL path placeholders.
            name: label shown in the progress bar.

        Returns:
            Path: The path to the downloaded asset.
        """
        # TODO: dry
        if operation_id not in self.endpoints_mapping:
            raise ValueError(
                f"Operation ID '{operation_id}' is not present in the server's endpoints."
                f"Available operation IDs: {self.endpoints_mapping.keys()}"
            )
        endpoint_template = self.endpoints_mapping[operation_id]
        endpoint = (
            endpoint_template.format(**url_format_params)
            if url_format_params
            else endpoint_template
        )
        url = urljoin(self.server_url, endpoint)

        if not output_path.parent.exists():
            output_path.parent.mkdir(parents=True, exist_ok=True)

        label = name or output_path.name

        # Download to a temp file beside the final destination, then rename
        # into place only once the whole body has arrived. Writing straight to
        # output_path left a truncated file sitting where a complete one
        # belongs whenever a download was interrupted or failed partway.
        tmp_path = output_path.with_name(output_path.name + ".part")

        def _stream(response, reporter):
            """Write the response body to disk, reporting bytes as they arrive."""
            response.raise_for_status()
            content_length_header = response.headers.get("content-length")
            content_length = int(content_length_header) if content_length_header else None
            reporter.start(content_length)
            for chunk in response.iter_bytes():
                download_file.write(chunk)
                reporter.update(response.num_bytes_downloaded)
            reporter.finish(response.num_bytes_downloaded)

        def _open_stream(reporter):
            if self._http_client:
                # Use injected client (for testing)
                with self._http_client.stream(
                    "GET", url, follow_redirects=True, params=params
                ) as response:
                    _stream(response, reporter)
            else:
                # Use default httpx client. Explicit timeout: httpx's own
                # 5-second default read timeout applies per-chunk, which a
                # genome-scale archive trips under any load.
                with httpx.stream(
                    "GET",
                    url,
                    follow_redirects=True,
                    params=params,
                    timeout=_DOWNLOAD_TIMEOUT,
                ) as response:
                    _stream(response, reporter)

        try:
            if progress.active_sink() is None:
                # No consumer: the terminal is the consumer. Unchanged CLI behavior.
                with (
                    tmp_path.open("wb") as download_file,
                    Progress(*RICH_PROGRESS_COLUMNS) as rich_progress,
                ):
                    _open_stream(_RichDownloadReporter(rich_progress, label))
            else:
                # A sink is installed, so DO NOT build a rich Progress here.
                #
                # Two reasons, and the second is load-bearing. First, a live rich
                # display under uvicorn writes ANSI control sequences into the
                # server console. Second, rich permits only ONE live display per
                # console, so two overlapping downloads that each construct a
                # Progress raise LiveError -- which is exactly why the job manager
                # can run pulls concurrently: it always installs a sink, so this
                # branch is always the one taken in job context, and the
                # single-concurrent-pull constraint disappears with it.
                progress.emit(
                    "stage", message=f"Downloading {label}", phase="download"
                )
                with tmp_path.open("wb") as download_file:
                    _open_stream(_SinkDownloadReporter(label))
        except BaseException:
            # Interrupted or failed: the partial file must not be mistaken for
            # a complete one, whether by a later pull run or by a caller
            # inspecting output_path.
            tmp_path.unlink(missing_ok=True)
            raise

        os.replace(tmp_path, output_path)
        return output_path

    def get_paginated(
        self,
        operation_id: str,
        params: dict | None = None,
        url_format_params: dict[str, str | None] = None,
        page_size: int = DEFAULT_PAGE_SIZE,
        max_items: int | None = None,
        strict: bool = True,
    ) -> list[dict]:
        """
        Get all results from a paginated endpoint by automatically iterating through pages.

        Args:
            operation_id (str): operation ID
            params (dict, optional): query parameters (excluding offset/limit)
            url_format_params (dict, optional): parameters to format URL path with
            page_size (int): number of items per page (default: 100, max: 1000)
            max_items (int, optional): maximum total items to retrieve
            strict (bool): re-raise if a page fails. When False, the pages
                collected so far are returned and a warning is logged -- only
                appropriate for display-only callers that can tolerate an
                incomplete list.

        Returns:
            list[Dict]: all results from all pages combined

        Raises:
            Exception: whatever the underlying request raised, when strict.
        """
        # Ensure page_size is within server limits
        page_size = min(max(page_size, 1), MAX_PAGE_SIZE)

        all_results = []
        offset = 0
        total_retrieved = 0

        while True:
            page_params = dict(params or {})
            page_params["offset"] = offset
            page_params["limit"] = page_size

            # If max_items is set, adjust limit for final page
            if max_items is not None:
                remaining = max_items - total_retrieved
                if remaining <= 0:
                    break
                page_params["limit"] = min(page_size, remaining)

            logger.debug(f"Fetching page: offset={offset}, limit={page_params['limit']}")

            try:
                response = self.get(
                    operation_id=operation_id,
                    params=page_params,
                    url_format_params=url_format_params,
                )

                # Handle different response formats
                if isinstance(response, dict):
                    if "items" in response:
                        # PaginatedResponse format
                        items = response["items"]
                        supports_pagination = True
                    elif "data" in response:
                        # Alternative format
                        items = response["data"]
                        supports_pagination = True
                    else:
                        # Direct list format (legacy) - server doesn't support pagination
                        items = response if isinstance(response, list) else [response]
                        supports_pagination = False
                elif isinstance(response, list):
                    # Direct list response - server doesn't support pagination
                    items = response
                    supports_pagination = False
                else:
                    logger.warning(f"Unexpected response format: {type(response)}")
                    break

                if not items:
                    # No more results
                    break

                all_results.extend(items)
                total_retrieved += len(items)

                # If server doesn't support pagination, it returned all results at once
                if not supports_pagination:
                    logger.debug(
                        f"Server doesn't support pagination, retrieved all {total_retrieved} items"
                    )
                    break

                offset += len(items)

                logger.debug(f"Retrieved {len(items)} items, total so far: {total_retrieved}")

                # Check if we've got everything
                if len(items) < page_params["limit"]:
                    # This page had fewer items than requested, so we're at the end
                    break

                if max_items is not None and total_retrieved >= max_items:
                    break

            except Exception as e:
                if strict:
                    # A truncated collection is indistinguishable from a small
                    # one; callers must not mistake it for the whole thing.
                    logger.error(f"Error fetching page at offset {offset}: {e}")
                    raise
                logger.warning(
                    f"Error fetching page at offset {offset}: {e}. "
                    f"Returning {len(all_results)} item(s) collected so far -- "
                    f"this list is INCOMPLETE."
                )
                break

        return all_results

    def get_all_assets(self, params: dict | None = None, strict: bool = True) -> list[dict]:
        """Get all assets using pagination."""
        return self.get_paginated("list_assets_v4_assets_get", params, strict=strict)

    def get_all_asset_groups(self, params: dict | None = None, strict: bool = True) -> list[dict]:
        """Get all asset groups using pagination."""
        return self.get_paginated("list_asset_groups_v4_asset_groups_get", params, strict=strict)

    def get_all_genomes(self, params: dict | None = None, strict: bool = True) -> list[dict]:
        """Get all genomes using pagination."""
        return self.get_paginated("list_genomes_v4_genomes_get", params, strict=strict)

    def get_all_aliases(self, params: dict | None = None, strict: bool = True) -> list[dict]:
        """Get all aliases using pagination."""
        return self.get_paginated("list_aliases_v4_aliases_get", params, strict=strict)

    def get_all_staged_assets(self, params: dict | None = None, strict: bool = True) -> list[dict]:
        """Get all staged assets using pagination."""
        if not self.has_endpoint(API_ID_STAGED_ASSETS):
            logger.debug(
                f"Server {self.server_url} does not expose staged_assets endpoint; skipping"
            )
            return []
        return self.get_paginated(API_ID_STAGED_ASSETS, params, strict=strict)

    def get_asset_groups(self, params: dict | None = None) -> list[dict]:
        """Get one page of asset groups, unwrapping the paginated response envelope."""
        response = self.get("list_asset_groups_v4_asset_groups_get", params)
        return self.extract_items_from_response(response)

    def get_assets(self, params: dict | None = None) -> list[dict]:
        """Get one page of assets, unwrapping the paginated response envelope."""
        response = self.get("list_assets_v4_assets_get", params)
        return self.extract_items_from_response(response)

    def get_staged_assets(self, params: dict | None = None) -> list[dict]:
        """Get one page of staged assets, unwrapping the paginated response envelope."""
        if not self.has_endpoint(API_ID_STAGED_ASSETS):
            logger.debug(
                f"Server {self.server_url} does not expose staged_assets endpoint; skipping"
            )
            return []
        response = self.get(API_ID_STAGED_ASSETS, params)
        return self.extract_items_from_response(response)

    def get_asset_file_list(self, asset_digest: str) -> list[str]:
        """Get the list of files available for a file-mode asset.

        Args:
            asset_digest: The digest of the asset.

        Returns:
            list[str]: List of relative file paths within the asset.
        """
        response = self.get(
            operation_id=API_ID_ASSET_FILES,
            url_format_params={"asset_digest": asset_digest},
        )
        if isinstance(response, dict):
            return response.get("files", [])
        return []

    def download_file(self, asset_digest: str, file_path: str, output_path: Path) -> Path:
        """Download a single file from a file-mode asset.

        Args:
            asset_digest: The digest of the asset.
            file_path: The relative file path within the asset.
            output_path: The local path to save the file to.

        Returns:
            Path: The path to the downloaded file.
        """
        return self.download_with_progress(
            operation_id=API_ID_ASSET_FILE_DOWNLOAD,
            output_path=output_path,
            url_format_params={"asset_digest": asset_digest, "file_path": file_path},
            name=file_path,
        )
