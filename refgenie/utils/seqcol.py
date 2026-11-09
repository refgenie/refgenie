"""Seqcol level 2 and comparison, shared by the server, the CLI and MCP.

Every compare path (two genomes in one store, in two stores, or on a remote
store) builds each side's level 2 with :func:`seqcol_level2` and compares with
:func:`compare_level2`, so all of them return the same GA4GH comparison.
"""

from typing import Any

from refget.utils import compare_seqcols

from refgenie.logger import logger

__all__ = ["compare_level2", "compare_in_stores", "seqcol_level2"]


def seqcol_level2(store: Any, digest: str) -> dict:
    """A collection's level 2 with every attribute the store holds.

    The store's level 2 has only names, lengths and sequences; the derived
    attributes (such as ``name_length_pairs``) are resolved from level 1. An
    attribute that cannot be resolved is left out with a warning, so the
    comparison still covers the core attributes.

    Raises:
        ValueError: If the store has no level 2 for ``digest``.
    """
    level2 = store.get_collection_level2(digest)
    if level2 is None:
        raise ValueError(f"No collection metadata found for {digest}")
    level2 = dict(level2)
    try:
        level1 = dict(store.get_collection_level1(digest))
    except Exception as e:  # noqa: BLE001 - compare the core attributes only
        logger.warning(f"Could not read level 1 of {digest}; comparing core attributes: {e}")
        return level2
    for attr, attr_digest in level1.items():
        if attr in level2:
            continue
        try:
            resolved = store.get_attribute(attr, attr_digest)
        except Exception as e:  # noqa: BLE001 - leave this attribute out
            logger.warning(f"Could not resolve attribute '{attr}' of {digest}: {e}")
            continue
        if resolved is not None:
            level2[attr] = resolved
    return level2


def compare_level2(level2_a: dict, level2_b: dict, digest_a: str, digest_b: str) -> dict:
    """The GA4GH seqcol comparison of two level 2 collections."""
    result = compare_seqcols(level2_a, level2_b)
    result.setdefault("digests", {"a": digest_a, "b": digest_b})
    return result


def compare_in_stores(store_a: Any, digest_a: str, store_b: Any, digest_b: str) -> dict:
    """Compare two collections, each read from its own store (the same one is fine)."""
    return compare_level2(
        seqcol_level2(store_a, digest_a), seqcol_level2(store_b, digest_b), digest_a, digest_b
    )
