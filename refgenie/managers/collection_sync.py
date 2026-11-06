"""Register every collection a remote source or a federated store lists.

``sync_collections`` pages through a source and hands each collection to a
register function. ``sync_source`` (``refgenie genome sync``) and
``sync_store`` (``refgenie store sync``) supply the two registration policies.
"""

from collections.abc import Callable
from dataclasses import dataclass

from refgenie.db.tables import Store, StoreType
from refgenie.exceptions import MissingGenomeError
from refgenie.logger import logger
from refgenie.managers.alias import AliasManager

__all__ = [
    "SYNC_ALIAS_NAMESPACE",
    "SyncReport",
    "collection_field",
    "sync_collections",
    "sync_source",
    "sync_store",
]

#: The only store alias namespace sync copies into refgenie's flat aliases: it
#: is 1:1 with collections. Others (genome_assembly, refseq, ...) hold
#: duplicated names, like 'hg38', that would collide with built genomes.
SYNC_ALIAS_NAMESPACE = "name"


@dataclass
class SyncReport:
    """New genomes registered, and collections or listings that failed."""

    registered: int = 0
    failures: int = 0


Register = Callable[[str, str, SyncReport], None]


def sync_collections(source, label: str, page_size: int, register: Register) -> SyncReport:
    """Page through ``source.list_collections`` and register each collection.

    ``register(digest, description, report)`` records its own outcome on the
    report; an exception it raises counts as one failure. A failed listing
    counts as one failure and ends the sync.
    """
    report = SyncReport()
    page = 0
    while True:
        try:
            result = source.list_collections(page=page, page_size=page_size)
        except Exception as e:  # noqa: BLE001 - reported, not raised
            logger.error(f"Failed to list collections from {label}: {e}")
            report.failures += 1
            break
        items = (result or {}).get("results", [])
        for item in items:
            digest = collection_field(item, "digest")
            if not digest:
                continue
            try:
                register(digest, collection_field(item, "description") or "", report)
            except Exception as e:  # noqa: BLE001 - one bad collection never stops a sync
                logger.error(f"Failed to sync collection {digest} from {label}: {e}")
                report.failures += 1
        if len(items) < page_size:
            break
        page += 1
    return report


def collection_field(item, name: str):
    """A listing entry's field: refget lists dicts, stores list objects."""
    return item.get(name) if isinstance(item, dict) else getattr(item, name, None)


def sync_source(rgc, source, label: str, page_size: int) -> SyncReport:
    """Register each collection on a refgenie server or RefgetStore as a genome.

    Idempotent: an existing genome is kept, its missing aliases are added, and
    its FHR metadata is applied again so source-side edits propagate.
    """

    def register(digest: str, description: str, report: SyncReport) -> None:
        alias_names = _syncable_alias_names(rgc, digest, source.get_collection_aliases(digest))
        digest, created = rgc.genome.initialize_genome(
            source=source,
            digest=digest,
            description=description,
            alias_names=alias_names,
            use_existing=True,
        )
        if created:
            report.registered += 1
            logger.info(f"Registered genome: {digest}")
        # The genome stays registered when its FHR fails; that still counts.
        try:
            fhr = source.get_collection_fhr(digest)
            if fhr:
                rgc.genome.apply_fhr(digest, fhr)
        except Exception as e:  # noqa: BLE001
            logger.error(f"Genome {digest} registered but FHR apply failed: {e}")
            report.failures += 1

    return sync_collections(source, label, page_size, register)


def _syncable_alias_names(rgc, digest: str, alias_pairs: list[tuple[str, str]]) -> list[str]:
    """A collection's aliases in SYNC_ALIAS_NAMESPACE that no other genome holds.

    Sync never repoints an alias; one owned by a different genome is skipped.
    """
    names = []
    for namespace, name in alias_pairs:
        if namespace != SYNC_ALIAS_NAMESPACE or name in names:
            continue
        try:
            existing = rgc.alias.resolve(name)
        except Exception:  # noqa: BLE001 - unresolvable reads as unclaimed
            existing = None
        if existing is not None and existing != digest:
            logger.warning(
                f"Sync skipping alias '{name}' for {digest}: it already points "
                f"to {existing} and sync never repoints an alias."
            )
            continue
        names.append(name)
    return names


def sync_store(
    rgc, store, store_row: Store, priorities: dict[str, int], page_size: int
) -> SyncReport:
    """Ingest one federated store's collections, aliases and FHR columns.

    ``store`` is the opened refget store for ``store_row``. A genome in several
    stores is owned by the highest-priority one (lowest number); alias ties are
    settled the same way by ``AliasManager.federated_sync``.
    """
    # The SQL alias backend whatever the instance's mode: the federation
    # registry is a SQL concept the server reads.
    sql_alias = AliasManager(rgc.database_engine)

    def register(digest: str, description: str, report: SyncReport) -> None:
        _claim_genome(rgc, digest, description, store_row, priorities)
        report.registered += _ingest_aliases(sql_alias, store, digest, store_row, priorities)
        _ingest_fhr(rgc, store, digest, store_row)

    return sync_collections(store, f"store '{store_row.name}'", page_size, register)


def _claim_genome(rgc, digest, description, store_row: Store, priorities) -> None:
    """Add the genome owned by this store, or take ownership if this store outranks the owner."""
    try:
        genome = rgc.genome.get(digest)
    except MissingGenomeError:
        rgc.genome.add(
            digest=digest, description=description, alias_names=[], store_name=store_row.name
        )
        return
    current = genome.store_name
    if current == store_row.name:
        return
    if current is None or store_row.priority < priorities.get(current, float("inf")):
        if current is not None:
            logger.info(
                f"Genome {digest}: store '{store_row.name}' (priority {store_row.priority}) "
                f"outranks owner '{current}'; repointing ownership."
            )
        rgc.genome.set_store_name(digest, store_row.name)
    else:
        logger.debug(
            f"Genome {digest}: keeping owner '{current}'; '{store_row.name}' has lower priority."
        )


def _ingest_aliases(sql_alias, store, digest, store_row: Store, priorities) -> int:
    """Apply the collision policy to each of a collection's 'name' aliases.

    Returns 1 if any alias was newly inserted, else 0: a coarse count for the
    summary line.
    """
    try:
        alias_pairs = list(store.get_aliases_for_collection(digest))
    except Exception as e:  # noqa: BLE001 - a collection with no aliases is normal
        logger.debug(f"No aliases for {digest} in '{store_row.name}': {e}")
        return 0
    inserted = 0
    for namespace, name in alias_pairs:
        if namespace != SYNC_ALIAS_NAMESPACE:
            continue
        if sql_alias.federated_sync(name, digest, store_row.name, priorities) == "inserted":
            inserted = 1
    return inserted


def _ingest_fhr(rgc, store, digest, store_row: Store) -> None:
    """Write the queryable FHR columns from the store's sidecar, if any.

    A remote store opens without its FHR sidecars, so this genome's is pulled first.
    """
    if store_row.type != StoreType.on_disk:
        try:
            store.pull_fhr(digest)
        except Exception as e:  # noqa: BLE001 - a genome with no FHR is normal
            logger.debug(f"pull_fhr({digest}) failed on store '{store_row.name}': {e}")
    fhr = store.get_fhr_metadata(digest)
    if fhr is not None:
        rgc.genome.apply_fhr_columns(digest, fhr.to_dict())
