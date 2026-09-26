"""
Refgenie managers for domain-specific operations.

Managers handle CRUD operations for their respective domains:
- AliasManager: genome alias operations
- GenomeManager: genome operations
- AssetManager: asset operations
- StageManager: staging operations (archive or folder serving modes)
- AssetClassManager: asset class operations
- ConfigurationManager: configuration operations
- RemoteManager: push targets (remotes) and the staged assets queued for each
- SourceManager: data channels, and syncing their recipes and asset classes
- RecipeManager: recipe operations
- DatabaseManager: schema creation, backend init, migrations, and purge
- SequenceManager: sequence retrieval through the RefgetStore router
- TransferManager: pulls from subscribed servers (one asset, many genomes, mirror)
- BuildManager: building assets from recipes, build preflight, and Snakemake targets
"""

from refgenie.managers.base import ResourceManager
from refgenie.managers.alias import (
    AliasBackend,
    AliasManager,
    FederatedAliasManager,
    StoreAliasManager,
)
from refgenie.managers.asset import AssetManager
from refgenie.managers.genome import GenomeManager
from refgenie.managers.stage import StageManager
from refgenie.managers.asset_class import AssetClassManager
from refgenie.managers.configuration import ConfigurationManager
from refgenie.managers.remote import RemoteManager
from refgenie.managers.sources import SourceManager
from refgenie.managers.recipe import RecipeManager
from refgenie.managers.store import StoreManager
from refgenie.managers.database import DatabaseManager
from refgenie.managers.sequence import SequenceManager
from refgenie.managers.transfer import TransferManager
from refgenie.managers.build import BuildManager

__all__ = [
    "ResourceManager",
    "AliasBackend",
    "AliasManager",
    "FederatedAliasManager",
    "StoreAliasManager",
    "AssetManager",
    "GenomeManager",
    "StageManager",
    "AssetClassManager",
    "ConfigurationManager",
    "RemoteManager",
    "SourceManager",
    "RecipeManager",
    "StoreManager",
    "DatabaseManager",
    "SequenceManager",
    "TransferManager",
    "BuildManager",
]
