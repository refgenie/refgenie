"""
One way to build the httpx clients that talk to remote refgenie servers.

Five call sites used to construct their own, and they disagreed: some sent the
package User-Agent and some sent httpx's default (which Cloudflare, in front of
api.refgenie.org, answers with a 403), some followed redirects and some did
not, and the two channel handlers each rebuilt the same token/basic-auth
selection. `make_client` settles all of it in one place.

``httpx`` is imported here at module level, so import this module from inside a
function if the caller is on a hot import path.
"""

from typing import Any

import httpx

from refgenie.const import HTTP_HEADERS

#: Long enough for a real download to get going, short enough that an
#: unreachable server fails the command rather than hanging it.
DEFAULT_TIMEOUT = 10.0


def make_client(
    *,
    timeout: float = DEFAULT_TIMEOUT,
    credentials: dict | None = None,
    **kwargs: Any,
) -> httpx.Client:
    """
    Build an httpx client with refgenie's headers, timeout and redirect policy.

    Args:
        timeout: Seconds to wait. Lower it for a probe that must fail fast.
        credentials: Optional channel credentials. A ``token`` becomes a bearer
            ``Authorization`` header; a ``username``/``password`` pair becomes
            basic auth. Anything else is ignored.
        **kwargs: Passed to `httpx.Client`.

    Returns:
        httpx.Client: The client, unopened. Use it as a context manager.
    """
    headers = dict(HTTP_HEADERS)
    auth = None
    if credentials:
        if credentials.get("token"):
            headers["Authorization"] = f"Bearer {credentials['token']}"
        elif credentials.get("username") and credentials.get("password"):
            auth = httpx.BasicAuth(credentials["username"], credentials["password"])
    return httpx.Client(
        headers=headers,
        auth=auth,
        timeout=timeout,
        follow_redirects=True,
        **kwargs,
    )
