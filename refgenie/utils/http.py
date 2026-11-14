"""
One way to build the httpx2 clients that talk to remote refgenie servers.

Build every remote-server client with `make_client`, so all send the package
User-Agent (Cloudflare blocks httpx2's default), follow redirects, and share one
channel auth rule. ``httpx2`` is imported at module level, so import this module
lazily from hot import paths.
"""

from typing import Any

import httpx2

from refgenie.const import HTTP_HEADERS

#: Long enough for a real download to get going, short enough that an
#: unreachable server fails the command rather than hanging it.
DEFAULT_TIMEOUT = 10.0


def make_client(
    *,
    timeout: float = DEFAULT_TIMEOUT,
    credentials: dict | None = None,
    **kwargs: Any,
) -> httpx2.Client:
    """
    Build an httpx2 client with refgenie's headers, timeout and redirect policy.

    Args:
        timeout: Seconds to wait. Lower it for a probe that must fail fast.
        credentials: Optional channel credentials. A ``token`` becomes a bearer
            ``Authorization`` header; a ``username``/``password`` pair becomes
            basic auth. Anything else is ignored.
        **kwargs: Passed to `httpx2.Client`.

    Returns:
        httpx2.Client: The client, unopened. Use it as a context manager.
    """
    headers = dict(HTTP_HEADERS)
    auth = None
    if credentials:
        if credentials.get("token"):
            headers["Authorization"] = f"Bearer {credentials['token']}"
        elif credentials.get("username") and credentials.get("password"):
            auth = httpx2.BasicAuth(credentials["username"], credentials["password"])
    return httpx2.Client(
        headers=headers,
        auth=auth,
        timeout=timeout,
        follow_redirects=True,
        **kwargs,
    )
