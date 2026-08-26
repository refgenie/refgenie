"""
Asset management submodule.

Exports only the public interface - AssetManager.
Internal classes (AssetBuilder, AssetPuller, AssetRelations) are hidden.
"""

from refgenie.managers.asset.manager import AssetManager

__all__ = ["AssetManager"]
