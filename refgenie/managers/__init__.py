"""
Refgenie managers for domain-specific operations.

Managers handle CRUD operations for their respective domains:
- AliasManager: genome alias operations
- GenomeManager: genome operations
- AssetManager: asset operations
- StageManager: staging operations (archive or folder serving modes)
- AssetClassManager: asset class operations
- ConfigurationManager: configuration operations
- SourceManager: external data sources (servers and data channels)
- RecipeManager: recipe operations
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
from refgenie.managers.sources import SourceManager
from refgenie.managers.recipe import RecipeManager
from refgenie.managers.store import StoreManager

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
    "SourceManager",
    "RecipeManager",
    "StoreManager",
]
