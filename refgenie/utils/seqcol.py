"""Seqcol level 2 and comparison, shared by the server, the CLI and MCP.

Every compare path (two genomes in one store, in two stores, or on a remote
store) builds each side's level 2 with :func:`seqcol_level2` and compares with
:func:`compare_level2`, so all of them return the same GA4GH comparison.
"""

from typing import Any

from refget.utils import compare_seqcols

from refgenie.logger import logger

__all__ = ["compare_level2", "seqcol_level2", "summarize_comparison"]


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


def _names_by_digest(level2: dict) -> dict[str, list[str]]:
    by_digest: dict[str, list[str]] = {}
    for name, digest in zip(level2["names"], level2["sequences"]):
        by_digest.setdefault(digest, []).append(name)
    return by_digest


def _pairs(renames: dict[str, str]) -> str:
    """``a->b`` pairs, at most 4; more shows the first 3, ``...`` and the last."""
    pairs = [f"{a}->{b}" for a, b in renames.items()]
    if len(pairs) > 4:
        pairs = pairs[:3] + ["...", pairs[-1]]
    return ", ".join(pairs)


def summarize_comparison(
    level2_a: dict, level2_b: dict, comparison: dict, label_a: str, label_b: str
) -> tuple[str, dict[str, str]]:
    """One sentence on how two genomes differ, and a map of renamed sequences.

    The map takes A's name to B's name for each sequence the two share by
    content but not by name, in A's order.
    """
    counts = comparison["array_elements"]
    n_a, n_b = counts["a_count"]["sequences"], counts["b_count"]["sequences"]
    shared = counts["a_and_b_count"]["sequences"]
    if "name_length_pairs" in counts["a_and_b_count"]:
        same_nl = counts["a_and_b_count"]["name_length_pairs"]
    else:
        pairs_a = set(zip(level2_a["names"], level2_a["lengths"]))
        same_nl = len(pairs_a & set(zip(level2_b["names"], level2_b["lengths"])))

    names_a, names_b = _names_by_digest(level2_a), _names_by_digest(level2_b)
    name_map: dict[str, str] = {}
    for digest, names in names_a.items():
        for name_a, name_b in zip(names, names_b.get(digest, [])):
            if name_a != name_b:
                name_map[name_a] = name_b

    digests = comparison.get("digests") or {}
    if level2_a == level2_b or (digests.get("a") and digests.get("a") == digests.get("b")):
        return "Identical.", {}

    if name_map:
        which = f"all {shared}" if len(name_map) == shared else str(len(name_map))
        naming = f"; {which} names differ: {_pairs(name_map)}."
    else:
        naming = "; same names."
    smaller, other = (label_a, label_b) if n_a <= n_b else (label_b, label_a)
    low = min(n_a, n_b)

    if shared == n_a == n_b:
        same_order = counts["a_and_b_same_order"]["sequences"] is not False
        order = "in the same order" if same_order else "in a different order"
        return f"Same {n_a} sequences {order}{naming}", name_map
    if same_nl == low and shared < same_nl:
        differ = same_nl - shared
        if n_a == n_b:
            return (
                f"Same {n_a} names and lengths, but {differ} of {n_a} sequences differ in content.",
                name_map,
            )
        return (
            f"Every name and length of {smaller} is also in {other}, but {differ} of those "
            f"{same_nl} sequences differ in content.",
            name_map,
        )
    if shared == low:
        return (
            f"All {shared} sequences of {smaller} are in {other}, which has "
            f"{max(n_a, n_b) - shared} more{naming}",
            name_map,
        )
    if shared > 0:
        return (
            f"{shared} sequences are shared ({label_a} has {n_a}, {label_b} has {n_b}){naming}",
            name_map,
        )
    return "No sequences in common.", {}
