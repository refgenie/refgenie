"""
One manager per kind of thing refgenie keeps: genomes, aliases, assets, recipes,
stages, sources, transfers, builds and more. The ``Refgenie`` root creates them
and exposes each as an attribute, such as ``rgc.genome`` or ``rgc.asset``.

Managers sit between the root and the database. They call the query functions
in ``refgenie.db.queries`` and never write SQL themselves. Each module and
subpackage describes its own manager.
"""

from refgenie.managers.base import ResourceManager
from refgenie.managers.alias import (
    AliasBackend,
    AliasManager,
    FederatedAliasManager,
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
from refgenie.managers.species import SpeciesManager
from refgenie.managers.database import DatabaseManager
from refgenie.managers.sequence import SequenceManager
from refgenie.managers.transfer import TransferManager
from refgenie.managers.build import BuildManager

__all__ = [
    "ResourceManager",
    "AliasBackend",
    "AliasManager",
    "FederatedAliasManager",
    "AssetManager",
    "GenomeManager",
    "StageManager",
    "AssetClassManager",
    "ConfigurationManager",
    "RemoteManager",
    "SourceManager",
    "RecipeManager",
    "StoreManager",
    "SpeciesManager",
    "DatabaseManager",
    "SequenceManager",
    "TransferManager",
    "BuildManager",
]
