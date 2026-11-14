from __future__ import annotations

import json
import logging
import os
import time
from functools import cached_property
from pathlib import Path
from typing import ClassVar
from urllib.parse import urljoin

import httpx2
from rich.progress import Progress

from refgenie import progress
from refgenie.config import config
from refgenie.const import API_VERSION, DEFAULT_PAGE_SIZE, MAX_PAGE_SIZE
from refgenie.logger import logger
from refgenie.managers.sources.api_ids import (
    API_ID_ALIAS_DIGEST,
    API_ID_ALIASES,
    API_ID_ASSET_CLASS,
    API_ID_ASSET_FILE_DOWNLOAD,
    API_ID_ASSET_FILES,
    API_ID_ASSET_GROUPS,
    API_ID_ASSET_RELATIONSHIPS,
    API_ID_ASSETS,
    API_ID_GENOME_ATTRS,
    API_ID_GENOME_COMPARE,
    API_ID_GENOMES,
    API_ID_STAGED_ASSETS,
)
from refgenie.managers.sources.openapi import RefgenieserverOpenApiSpec
from refgenie.utils.console import CONSOLE, RICH_PROGRESS_COLUMNS
from refgenie.utils.http import make_client

# One INFO line per HTTP request drowns refgenie's own INFO progress; show
# them only when debugging.
logging.getLogger("httpx2").setLevel(
    logging.DEBUG if config.log_level.value == "DEBUG" else logging.WARNING
)


#: Report to a progress sink at most this often, in seconds...
_SINK_MIN_INTERVAL_SECONDS = 0.2
#: ...or every this many bytes, whichever comes first.
_SINK_MIN_BYTES = 1024 * 1024

#: Timeout for archive/file downloads; httpx2's 5s default read timeout is too
#: short for large archive chunks under load.
_DOWNLOAD_TIMEOUT = httpx2.Timeout(connect=10, read=60, write=60, pool=None)


class _RichDownloadReporter:
    """Drive a `rich` progress bar on the terminal."""

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

    Throttled so per-chunk events do not flood the job manager's ring buffer.
    `total` is None when the server sent no Content-Length.
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
        http_client: httpx2.Client | None = None,
    ):
        self.server_url = server_url
        self.openapi_endpoint = openapi_endpoint
        self.api_version = api_version
        self._http_client = http_client

        self.asset = AssetClient(self)
        self.asset_group = AssetGroupClient(self)
        self.asset_class = AssetClassClient(self)
        self.genome = GenomeClient(self)
        self.alias = AliasClient(self)
        self.staged_asset = StagedAssetClient(self)

        logger.debug(f"Connecting to the server at {self.server_url}")

        # Fetch the spec eagerly so a connection failure surfaces here.
        logger.info(f"Connected to server: {self.server_url} {self.openapi_spec.info}")

    def __repr__(self):
        return f"RefgenieserverClient(server_url={self.server_url})"

    @cached_property
    def openapi_spec(self) -> RefgenieserverOpenApiSpec:
        """
        The server's validated OpenAPI description, fetched once.

        Returns:
            RefgenieserverOpenApiSpec: openAPI JSON description
        """
        logger.debug(f"Retrieving OpenAPI spec from {self.openapi_endpoint}")
        return RefgenieserverOpenApiSpec.model_validate(
            self._get(urljoin(self.server_url, self.openapi_endpoint))
        )

    @cached_property
    def endpoints_mapping(self) -> dict:
        """
        Server endpoint paths keyed by operationId.

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

    def _url(self, operation_id: str, url_format_params: dict[str, str | None] | None) -> str:
        """The absolute URL of ``operation_id``, with its path placeholders filled."""
        if operation_id not in self.endpoints_mapping:
            raise ValueError(
                f"Operation ID '{operation_id}' is not present in the server's endpoints. "
                f"Available operation IDs: {self.endpoints_mapping.keys()}"
            )
        endpoint = self.endpoints_mapping[operation_id]
        if url_format_params:
            endpoint = endpoint.format(**url_format_params)
        return urljoin(self.server_url, endpoint)

    def _send_request(self, url: str, params: dict | None = None) -> httpx2.Response:
        """Send a GET request to ``url`` and raise on an error status."""
        if self._http_client:
            response = self._http_client.get(url, params=params)
        else:
            with make_client() as client:
                response = client.get(url, params=params)
        response.raise_for_status()
        return response

    def _get(self, url: str, params: dict | None = None) -> dict | str:
        """GET ``url``: the JSON body, or the text if it is not JSON."""
        response = self._send_request(url, params=params)
        try:
            return response.json()
        except json.JSONDecodeError:
            return response.text

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
        return self._get(self._url(operation_id, url_format_params), params)

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

        Progress goes to a `rich` bar, or to the `refgenie.progress` sink when
        one is installed.

        Args:
            operation_id: operationId of the download endpoint.
            output_path: where to save the file.
            params: query parameters.
            url_format_params: values for URL path placeholders.
            name: label shown in the progress bar.

        Returns:
            Path: The path to the downloaded asset.
        """
        url = self._url(operation_id, url_format_params)

        if not output_path.parent.exists():
            output_path.parent.mkdir(parents=True, exist_ok=True)

        label = name or output_path.name

        # Download to a .part file and rename once complete, so an interrupted
        # download never leaves a truncated file at output_path.
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
                with self._http_client.stream(
                    "GET", url, follow_redirects=True, params=params
                ) as response:
                    _stream(response, reporter)
            else:
                with httpx2.stream(
                    "GET",
                    url,
                    follow_redirects=True,
                    params=params,
                    timeout=_DOWNLOAD_TIMEOUT,
                ) as response:
                    _stream(response, reporter)

        try:
            if progress.active_sink() is None:
                with (
                    tmp_path.open("wb") as download_file,
                    Progress(*RICH_PROGRESS_COLUMNS, console=CONSOLE) as rich_progress,
                ):
                    _open_stream(_RichDownloadReporter(rich_progress, label))
            else:
                # No rich display with a sink: rich allows one live display per
                # console, so concurrent job pulls would raise LiveError.
                progress.emit("stage", message=f"Downloading {label}", phase="download")
                with tmp_path.open("wb") as download_file:
                    _open_stream(_SinkDownloadReporter(label))
        except BaseException:
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
    ) -> list[dict]:
        """Fetch every page of a GA4GH-style paged endpoint and combine the results.

        Pages are requested from 0 with ``page_size`` clamped to 1..MAX_PAGE_SIZE,
        until ``pagination.total`` is reached or ``max_items`` results are collected.
        """
        page_size = min(max(page_size, 1), MAX_PAGE_SIZE)

        all_results = []
        page = 0

        while max_items is None or len(all_results) < max_items:
            page_params = dict(params or {}, page=page, page_size=page_size)
            logger.debug(f"Fetching page {page} (page_size={page_size})")

            response = self.get(
                operation_id=operation_id,
                params=page_params,
                url_format_params=url_format_params,
            )

            all_results.extend(response["results"])
            logger.debug(f"Retrieved page {page}, total so far: {len(all_results)}")
            if (page + 1) * page_size >= response["pagination"]["total"]:
                break
            page += 1

        if max_items is not None:
            all_results = all_results[:max_items]
        return all_results


class _ResourceClient:
    """One server resource, reached through the shared transport."""

    def __init__(self, transport: RefgenieserverClient):
        self._transport = transport


class _ListResourceClient(_ResourceClient):
    """A resource whose ``list()`` pages through one operation."""

    operation_id: ClassVar[str]

    def list(self, params: dict | None = None) -> list[dict]:
        """Every record matching ``params``, across all pages."""
        return self._transport.get_paginated(self.operation_id, params)


class AssetClient(_ListResourceClient):
    operation_id = API_ID_ASSETS

    def files(self, asset_digest: str) -> list[str]:
        """Relative paths of the files a file-mode asset serves."""
        response = self._transport.get(
            operation_id=API_ID_ASSET_FILES,
            url_format_params={"asset_digest": asset_digest},
        )
        if isinstance(response, dict):
            return response.get("files", [])
        return []

    def relationships(self, asset_digest: str) -> dict:
        """The asset's parents and children, expanded to full records."""
        return self._transport.get(
            operation_id=API_ID_ASSET_RELATIONSHIPS,
            url_format_params={"asset_digest": asset_digest},
            params={"expand": True},
        )

    def download_file(self, asset_digest: str, file_path: str, output_path: Path) -> Path:
        """Download one file of a file-mode asset to ``output_path``."""
        return self._transport.download_with_progress(
            operation_id=API_ID_ASSET_FILE_DOWNLOAD,
            output_path=output_path,
            url_format_params={"asset_digest": asset_digest, "file_path": file_path},
            name=file_path,
        )


class AssetGroupClient(_ListResourceClient):
    operation_id = API_ID_ASSET_GROUPS


class AssetClassClient(_ResourceClient):
    def get(self, asset_class_id: int) -> dict:
        """One asset class by the server's id for it.

        Older servers omit ``seek_keys``.

        Args:
            asset_class_id: The server's id of the asset class, as an asset
                group's ``asset_class_id`` gives it.

        Returns:
            dict: The asset class (``name``, ``version``, ``description``,
            ``serving_modes`` and, when sent, ``seek_keys``).
        """
        response = self._transport.get(
            operation_id=API_ID_ASSET_CLASS,
            url_format_params={"id": str(asset_class_id)},
        )
        if not isinstance(response, dict):
            raise ValueError(f"Unexpected asset class response from {self._transport.server_url}")
        return response


class GenomeClient(_ListResourceClient):
    operation_id = API_ID_GENOMES

    def get(self, digest: str) -> dict:
        """The server's record of one genome, or an empty dict if the body is not one."""
        response = self._transport.get(
            operation_id=API_ID_GENOME_ATTRS, url_format_params={"digest": digest}
        )
        return response if isinstance(response, dict) else {}

    def compare(self, digest: str, other: str) -> dict:
        """The server's GA4GH seqcol comparison of two genomes it serves."""
        response = self._transport.get(
            operation_id=API_ID_GENOME_COMPARE,
            url_format_params={"digest": digest, "other": other},
        )
        if not isinstance(response, dict):
            raise ValueError(f"Unexpected compare response from {self._transport.server_url}")
        return response


class AliasClient(_ListResourceClient):
    operation_id = API_ID_ALIASES

    def resolve(self, name: str) -> str | None:
        """The digest of the genome an alias names, or None if the body has none.

        Raises the HTTP error when the server does not resolve the alias.
        """
        response = self._transport.get(
            operation_id=API_ID_ALIAS_DIGEST, url_format_params={"name": name}
        )
        return response.get("digest") if isinstance(response, dict) else None


class StagedAssetClient(_ListResourceClient):
    operation_id = API_ID_STAGED_ASSETS
