"""The `store` command group: the federation registry.

``refgenie store`` manages the ``store`` table -- the single source of truth for
which refget stores this service federates over. ``store sync`` ingests each
store's collections and curated aliases into the SQL catalog so the server can
serve them as one, applying the priority-based alias collision policy.
"""

from collections.abc import Callable
from typing import Literal

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.errors import EXIT_NOT_FOUND, fail
from refgenie.logger import logger

# Only the curated 1:1 "name" namespace maps cleanly onto refgenie's flat alias
# space; the others (accession, refseq, genome_assembly, ...) are coarse or
# duplicated. Shared with genome sync's SYNC_ALIAS_NAMESPACE rationale.
SYNC_ALIAS_NAMESPACE = "name"


class StoreAddModel(BaseModel):
    """store add: register a store in the federation registry."""

    name: str = Field(description="Unique name for the store (e.g. 'jungle').")
    url: str = Field(
        description="Store root URL (remote) or local .refget_store path (on_disk).",
        validation_alias=AliasChoices("u", "url"),
    )
    type: Literal["remote", "on_disk"] = Field(
        default="remote",
        description="How the store is opened.",
        validation_alias=AliasChoices("t", "type"),
    )
    priority: int = Field(
        default=100,
        description="Lower integer = higher priority; wins alias ties.",
        validation_alias=AliasChoices("p", "priority"),
    )
    description: str | None = Field(default=None, description="Optional description.")


class StoreSyncModel(BaseModel):
    """store sync: ingest a store's collections and aliases into the catalog."""

    name: str | None = Field(
        default=None,
        description="Store to sync. If omitted, sync every enabled store.",
    )
    page_size: int = Field(
        default=100,
        description="Collections per page when listing the store.",
        validation_alias=AliasChoices("page-size", "page_size"),
    )


class StoreListModel(BaseModel):
    """store list: list all registered stores."""

    pass


class StoreRemoveModel(BaseModel):
    """store remove: remove a store from the registry."""

    name: str = Field(description="Store name to remove.")


class StoreConflictsModel(BaseModel):
    """store conflicts: show alias collisions across stores."""

    pass


class StoreParser(BaseModel):
    """Intermediate parser for store subcommands."""

    add: CliSubCommand[StoreAddModel] = Field(description="Add a store.")
    sync: CliSubCommand[StoreSyncModel] = Field(description="Ingest a store's genomes.")
    list: CliSubCommand[StoreListModel] = Field(description="List stores.")
    remove: CliSubCommand[StoreRemoveModel] = Field(description="Remove a store.")
    conflicts: CliSubCommand[StoreConflictsModel] = Field(description="Show alias collisions.")


def handle_store_add(cmd, refgenie) -> None:
    from refgenie.db.tables import StoreType
    from refgenie.exceptions import StoreExistsError

    store_type = StoreType(cmd.type)
    # Best-effort validation of a remote root: probe for its rgstore.json.
    if store_type == StoreType.remote:
        from refgenie.managers.sources import make_source

        try:
            make_source(cmd.url)
        except Exception as e:  # noqa: BLE001 - a temporarily-unreachable store may still be added
            logger.warning(
                f"Could not validate store root {cmd.url} ({e}); adding it anyway. "
                f"'store sync' will fail until it is reachable."
            )
    try:
        refgenie.store.add(
            name=cmd.name,
            url=cmd.url,
            store_type=store_type,
            priority=cmd.priority,
            description=cmd.description,
        )
    except StoreExistsError as e:
        fail(str(e))


def handle_store_remove(cmd, refgenie) -> None:
    from refgenie.exceptions import MissingStoreError

    try:
        refgenie.store.remove(cmd.name)
    except MissingStoreError as e:
        fail(str(e), EXIT_NOT_FOUND)


def handle_store_list(cmd, refgenie) -> None:
    stores = list(refgenie.store.list_all())
    if not stores:
        logger.info("No stores registered.")
        return
    rprint(refgenie.store.table())


def handle_store_conflicts(cmd, refgenie) -> None:
    """List alias collisions, i.e. losers stored under qualified ``store::alias`` names."""
    from sqlmodel import Session, select

    from refgenie.db.tables import Alias

    rows = []
    with Session(refgenie.database_engine) as session:
        for alias in session.exec(select(Alias)).all():
            if "::" in alias.name:
                store_name, _, bare = alias.name.partition("::")
                rows.append((bare, store_name, alias.genome_digest))
    if not rows:
        logger.info("No alias conflicts.")
        return
    from refgenie.utils.tables import build_table

    rprint(
        build_table(
            "Alias conflicts (loser mappings)",
            ["Alias", "Qualified store", "Genome digest"],
            rows,
        )
    )


def _sync_one_store(refgenie, store_row, priorities, page_size) -> tuple[int, int]:
    """Ingest one store's collections + aliases. Returns (genomes, failures)."""
    from refgenie.core.store_router import make_store
    from refgenie.managers.alias import AliasManager
    import tempfile
    from pathlib import Path

    # SQL alias backend regardless of the instance's mode: the federation
    # registry is a SQL concept the server reads.
    sql_alias = AliasManager(refgenie.database_engine)

    cache_dir = Path(tempfile.mkdtemp(prefix=f"refgenie_sync_{store_row.name}_"))
    store = make_store(store_row, cache_dir)

    registered = 0
    failures = 0
    page = 0
    while True:
        try:
            result = store.list_collections(page=page, page_size=page_size)
        except Exception as e:  # noqa: BLE001
            logger.error(f"Failed to list collections from store '{store_row.name}': {e}")
            failures += 1
            break
        if not result:
            break
        items = result.get("results", [])
        if not items:
            break

        for item in items:
            digest = item.get("digest") if isinstance(item, dict) else getattr(item, "digest", None)
            if not digest:
                continue
            description = (
                item.get("description", "") if isinstance(item, dict)
                else getattr(item, "description", "")
            )
            try:
                _upsert_genome_owner(refgenie, digest, description or "", store_row, priorities)
                registered += _ingest_aliases(sql_alias, store, digest, store_row, priorities)
                _ingest_fhr(refgenie, store, digest)
            except Exception as e:  # noqa: BLE001
                logger.error(f"Failed to sync collection {digest} from '{store_row.name}': {e}")
                failures += 1

        if len(items) < page_size:
            break
        page += 1

    return registered, failures


def _upsert_genome_owner(refgenie, digest, description, store_row, priorities) -> None:
    """Insert the genome (owned by this store) or resolve ownership by priority."""
    from refgenie.exceptions import MissingGenomeError

    try:
        genome = refgenie.genome.get(digest)
    except MissingGenomeError:
        refgenie.genome.add(
            digest=digest, description=description, alias_names=[], store_name=store_row.name
        )
        return

    current = genome.store_name
    if current is None or current == store_row.name:
        if current != store_row.name:
            refgenie.genome.set_store_name(digest, store_row.name)
        return

    cur_priority = priorities.get(current, float("inf"))
    if store_row.priority < cur_priority:
        logger.info(
            f"Genome {digest}: store '{store_row.name}' (priority {store_row.priority}) "
            f"outranks owner '{current}'; repointing ownership."
        )
        refgenie.genome.set_store_name(digest, store_row.name)
    else:
        logger.debug(
            f"Genome {digest}: keeping owner '{current}'; '{store_row.name}' has "
            f"lower priority."
        )


def _ingest_aliases(sql_alias, store, digest, store_row, priorities) -> int:
    """Apply the collision policy for each of a collection's 'name' aliases.

    Returns 1 if this was a newly registered genome-with-aliases (an insert),
    else 0 -- a coarse count for the summary line.
    """
    try:
        alias_pairs = list(store.get_aliases_for_collection(digest))
    except Exception as e:  # noqa: BLE001
        logger.debug(f"No aliases for {digest} in '{store_row.name}': {e}")
        return 0
    inserted = 0
    for namespace, name in alias_pairs:
        if namespace != SYNC_ALIAS_NAMESPACE:
            continue
        outcome = sql_alias.federated_sync(name, digest, store_row.name, priorities)
        if outcome == "inserted":
            inserted = 1
    return inserted


def _ingest_fhr(refgenie, store, digest) -> None:
    """Write the queryable FHR columns from the store's sidecar, if any."""
    fhr = store.get_fhr_metadata(digest)
    if fhr is not None:
        refgenie.genome.apply_fhr_columns(digest, fhr.to_dict())


def handle_store_sync(cmd, refgenie) -> None:
    if cmd.name:
        from refgenie.exceptions import MissingStoreError

        try:
            targets = [refgenie.store.get(cmd.name)]
        except MissingStoreError as e:
            fail(str(e), EXIT_NOT_FOUND)
    else:
        targets = refgenie.store.enabled_stores()
        if not targets:
            fail("No enabled stores to sync. Add one with 'refgenie store add'.")

    priorities = {s.name: s.priority for s in refgenie.store.list_all()}

    total = 0
    failures = 0
    for store_row in targets:
        logger.info(f"Syncing store '{store_row.name}' ({store_row.url})")
        registered, fails = _sync_one_store(refgenie, store_row, priorities, cmd.page_size)
        total += registered
        failures += fails
        logger.info(f"  {registered} new genome(s) from '{store_row.name}'")

    # Rebuild the router so a long-lived instance sees the freshly synced set.
    refgenie.reload_store_router()

    logger.info(f"Synced {total} new genome(s) total")
    if failures:
        fail(f"{failures} collection(s)/store(s) failed to sync")


STORE_DISPATCH: dict[type, Callable] = {
    StoreAddModel: handle_store_add,
    StoreSyncModel: handle_store_sync,
    StoreListModel: handle_store_list,
    StoreRemoveModel: handle_store_remove,
    StoreConflictsModel: handle_store_conflicts,
}


def handle_store_group(cmd, refgenie) -> None:
    leaf = get_subcommand(cmd, is_required=True)
    handler = STORE_DISPATCH.get(type(leaf))
    if handler is None:
        fail(f"Unknown store subcommand: {type(leaf).__name__}")
    handler(leaf, refgenie)
