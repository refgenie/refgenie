"""Seqcol level 2 and comparison, shared by the server, the CLI and MCP.

Every compare path (two genomes in one store, in two stores, or on a remote
store) builds each side's level 2 with :func:`seqcol_level2` and compares with
:func:`compare_level2`, so all of them return the same GA4GH comparison.
"""

from collections import Counter
from typing import Any

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


def _compare_seqcols(A: dict, B: dict) -> dict:
    """The GA4GH seqcol comparison function (copied from refget 0.13.0).

    Copied so a plain install does not need refget. Its output is the spec's
    comparison, which the server, the CLI and MCP all return, so it must not
    change.

    Args:
        A: Sequence collection A
        B: Sequence collection B

    Returns:
        dict: Following formal seqcol specification comparison function return value
    """
    a_keys = list(A.keys())
    b_keys = list(B.keys())
    a_keys.sort()
    b_keys.sort()

    all_keys = a_keys + list(set(b_keys) - set(a_keys))
    all_keys.sort()

    # Compute lengths of each array; only do this for array attributes
    a_lengths = {}
    b_lengths = {}
    for k in a_keys:
        a_lengths[k] = len(A[k])
    for k in b_keys:
        b_lengths[k] = len(B[k])

    return_obj = {
        "attributes": {"a_only": [], "b_only": [], "a_and_b": []},
        "array_elements": {
            "a_count": a_lengths,
            "b_count": b_lengths,
            "a_and_b_count": {},
            "a_and_b_same_order": {},
        },
    }

    for k in all_keys:
        if k not in A:
            return_obj["attributes"]["b_only"].append(k)
        elif k not in B:
            return_obj["attributes"]["a_only"].append(k)
        else:
            return_obj["attributes"]["a_and_b"].append(k)
            res = _compare_elements(A[k], B[k])
            return_obj["array_elements"]["a_and_b_count"][k] = res["a_and_b"]
            return_obj["array_elements"]["a_and_b_same_order"][k] = res["a_and_b_same_order"]
    return return_obj


def _compare_elements(A: list, B: list) -> dict[str, int | bool | None]:
    """Compare two arrays' elements for :func:`_compare_seqcols` (copied from refget 0.13.0)."""
    A_filtered = list(filter(lambda x: x in B, A))
    B_filtered = list(filter(lambda x: x in A, B))
    A_count = len(A_filtered)
    B_count = len(B_filtered)
    overlap = min(len(A_filtered), len(B_filtered))  # counts duplicates

    if A_count + B_count < 1:
        # order match requires at least 2 matching elements
        order = None
    elif not (A_count == B_count == overlap):
        # duplicated matches means order match is undefined
        order = None
    else:
        order = A_filtered == B_filtered
    return {"a_and_b": overlap, "a_and_b_same_order": order}


def compare_level2(level2_a: dict, level2_b: dict, digest_a: str, digest_b: str) -> dict:
    """The GA4GH seqcol comparison of two level 2 collections."""
    result = _compare_seqcols(level2_a, level2_b)
    result.setdefault("digests", {"a": digest_a, "b": digest_b})
    return result


def _names_by_digest(level2: dict) -> dict[str, list[str]]:
    by_digest: dict[str, list[str]] = {}
    for name, digest in zip(level2["names"], level2["sequences"]):
        by_digest.setdefault(digest, []).append(name)
    return by_digest


def _n(count: int, noun: str) -> str:
    """``1 sequence``, ``2 sequences``."""
    return f"{count} {noun}" if count == 1 else f"{count} {noun}s"


def _abbreviate(items: list[str]) -> str:
    """At most 4 items; more shows the first 3, ``...`` and the last."""
    if len(items) > 4:
        items = items[:3] + ["...", items[-1]]
    return ", ".join(items)


def _prefix_change(a: str, b: str) -> str | None:
    """The prefix dropped (``-chr``) or added (``+chr``) to make ``b`` from ``a``, if that is all."""
    if a.endswith(b):
        return "-" + a[: len(a) - len(b)]
    if b.endswith(a):
        return "+" + b[: len(b) - len(a)]
    return None


def _pairs(renames: dict[str, str]) -> str:
    """``a->b`` pairs: those with the most common prefix change (such as dropping
    ``chr``), abbreviated, then the rest (such as ``chrM->Mito``), so no odd one hides."""
    changes = {a: _prefix_change(a, b) for a, b in renames.items()}
    common = Counter(c for c in changes.values() if c).most_common(1)
    pattern = common[0][0] if common else None
    regular = [f"{a}->{b}" for a, b in renames.items() if changes[a] == pattern]
    irregular = [f"{a}->{b}" for a, b in renames.items() if changes[a] != pattern]
    return ", ".join(part for part in (_abbreviate(regular), _abbreviate(irregular)) if part)


def summarize_comparison(
    level2_a: dict, level2_b: dict, comparison: dict, label_a: str, label_b: str
) -> tuple[str, dict[str, str]]:
    """One sentence on how two genomes differ, and a map of renamed sequences.

    The map takes A's name to B's name for each sequence the two share by
    content but not by name, in A's order. A name both genomes use for
    different sequences adds a warning.
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
        k = len(name_map)
        which = "1 name differs" if k == 1 else f"{'all ' if k == shared else ''}{k} names differ"
        naming = f"; {which}: {_pairs(name_map)}."
    else:
        naming = "; same names."
    seq_b = dict(zip(level2_b["names"], level2_b["sequences"]))
    clash = [
        name
        for name, digest in zip(level2_a["names"], level2_a["sequences"])
        if seq_b.get(name, digest) != digest
    ]
    warning = (
        f" Warning: {_n(len(clash), 'name')} in both {'holds' if len(clash) == 1 else 'hold'} "
        f"different sequences: {_abbreviate(clash)}."
        if clash
        else ""
    )
    smaller, other = (label_a, label_b) if n_a <= n_b else (label_b, label_a)
    low = min(n_a, n_b)

    if shared == n_a == n_b:
        same_order = counts["a_and_b_same_order"]["sequences"] is not False
        order = "in the same order" if same_order else "in a different order"
        return f"Same {_n(n_a, 'sequence')} {order}{naming}{warning}", name_map
    if same_nl == low and shared < same_nl:
        differ = same_nl - shared
        verb = "differs" if differ == 1 else "differ"
        listed = f": {_abbreviate(clash)}." if clash else "."
        if n_a == n_b:
            return (
                f"Same names and lengths, but {differ} of {_n(n_a, 'sequence')} {verb} "
                f"in content{listed}",
                name_map,
            )
        return (
            f"Every name and length of {smaller} is also in {other}, but {differ} of those "
            f"{_n(same_nl, 'sequence')} {verb} in content{listed}",
            name_map,
        )
    if shared == low:
        if shared > 1:
            whole = f"All {shared} sequences of {smaller} are"
        else:
            whole = f"The one sequence of {smaller} is"
        return (
            f"{whole} in {other}, which has {max(n_a, n_b) - shared} more{naming}{warning}",
            name_map,
        )
    if shared > 0:
        return (
            f"{_n(shared, 'sequence')} {'is' if shared == 1 else 'are'} shared "
            f"({label_a} has {n_a}, {label_b} has {n_b}){naming}{warning}",
            name_map,
        )
    return f"No sequences in common.{warning}", {}
