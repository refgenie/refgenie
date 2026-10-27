"""
Where refgenie gets things from outside this node.

- ``ServerManager`` (``servers.py``, ``rgc.servers``): the refgenieservers this
  node pulls assets from -- subscriptions, clients, their catalogs, remote seek.
- ``SourceManager`` (``manager.py``, ``rgc.sources``): data channels, which
  publish recipes and asset classes, and ``sync_channel``, which registers
  them locally.
- ``genomes.py``: ``RemoteGenomeSource``, the store-backed view of a remote used
  to create genomes; ``client.py``: the HTTP client for a refgenieserver.
"""

from refgenie.managers.sources.client import RefgenieserverClient, ServerClient
from refgenie.managers.sources.genomes import (
    RemoteGenomeSource,
    make_source,
    normalize_server_url,
)
from refgenie.managers.sources.manager import IndexFile, SourceManager
from refgenie.managers.sources.servers import ServerManager, estimate_pull_size

__all__ = [
    "SourceManager",
    "ServerManager",
    "RefgenieserverClient",
    "ServerClient",
    "IndexFile",
    "RemoteGenomeSource",
    "estimate_pull_size",
    "make_source",
    "normalize_server_url",
]
