"""
Pulling from subscribed servers, reached as ``rgc.transfer``.

- ``manager.py`` (``TransferManager``): the public entry points -- one asset,
  many genomes, or a full mirror -- plus the bulk-pull confirmation.
- ``puller.py`` (``AssetPuller``): the pull engine for one asset, with its
  transaction, download modes, and large-archive and SIGINT policy.
"""

from refgenie.managers.transfer.manager import TransferManager

__all__ = ["TransferManager"]
