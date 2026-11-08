"""The `store` command group: the federation registry.

``store sync`` ingests each store's collections and curated aliases into the
SQL catalog, resolving alias collisions by store priority, and loads species
names from each enabled store's optional ``species.json``.
"""

from collections.abc import Callable
from pathlib import Path
from typing import Literal

from pydantic import AliasChoices, BaseModel, Field
from pydantic_settings import CliSubCommand, get_subcommand
from rich import print as rprint

from refgenie.cli.errors import EXIT_NOT_FOUND, fail
from refgenie.logger import logger


class StoreAddModel(BaseModel):
    """store add: register a store in the federation registry."""

    name: str = Field(description="Unique name for the store (e.g. 'jungle'). 'local' is reserved.")
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


class StoreApplyModel(BaseModel):
    """store apply: make the registry match a JSON declaration file."""

    file: Path = Field(
        description="JSON array of stores, each {name, url, type?, priority?, "
        "description?, enabled?}. Stores not listed are disabled, never deleted.",
        validation_alias=AliasChoices("f", "file"),
    )


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
    apply: CliSubCommand[StoreApplyModel] = Field(
        description="Make the store list match a declaration file."
    )
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
    except (StoreExistsError, ValueError) as e:
        fail(str(e))


def handle_store_apply(cmd, refgenie) -> None:
    import json

    try:
        declared = json.loads(Path(cmd.file).read_text())
    except (OSError, json.JSONDecodeError) as e:
        fail(f"Could not read store declaration {cmd.file}: {e}")
    if not isinstance(declared, list):
        fail(f"Store declaration {cmd.file} must be a JSON array of store objects")
    try:
        changes = refgenie.store.apply_declaration(declared, disable_missing=True)
    except ValueError as e:
        fail(f"Invalid store declaration: {e}")
    rows = [(kind, name) for kind, names in changes.items() for name in names]
    if not rows:
        logger.info("Stores already match the declaration; nothing changed.")
        return
    from refgenie.utils.tables import build_table

    rprint(build_table("Store changes", ["Change", "Store"], rows))


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
    rows = refgenie.store.conflicts()
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
    """Open one store and ingest it. Returns (new genomes, failures)."""
    import tempfile
    from pathlib import Path

    from refgenie.core.store_router import make_store
    from refgenie.db.tables import StoreType
    from refgenie.managers.collection_sync import sync_store

    with tempfile.TemporaryDirectory(prefix=f"refgenie_sync_{store_row.name}_") as cache_dir:
        store = make_store(store_row, Path(cache_dir), readonly=False)
        if store_row.type != StoreType.on_disk:
            # make_store opens a remote store with collection stubs only; sync
            # is the one reader of its alias sidecars.
            try:
                store.pull_aliases()
            except Exception as e:  # noqa: BLE001 - a store with no aliases is tolerated
                logger.debug(f"pull_aliases failed for store '{store_row.name}': {e}")
        report = sync_store(refgenie, store, store_row, priorities, page_size)
    return report.registered, report.failures


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

    # Species names feed search only, so they never fail a sync. All enabled
    # stores are read, not only the ones synced here.
    species = refgenie.species.load_from_stores(refgenie.store.enabled_stores())
    if species is not None:
        logger.info(f"Loaded {species} species")

    # Rebuild the router so a long-lived instance sees the freshly synced set.
    refgenie.reload_store_router()

    logger.info(f"Synced {total} new genome(s) total")
    if failures:
        fail(f"{failures} collection(s)/store(s) failed to sync")


STORE_DISPATCH: dict[type, Callable] = {
    StoreAddModel: handle_store_add,
    StoreApplyModel: handle_store_apply,
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
