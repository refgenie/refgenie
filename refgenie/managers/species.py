"""SpeciesManager - the ``species`` and ``speciesname`` tables.

Species names come from the stores' optional ``species.json`` files (see
``refgenie.utils.species_file``). One write rule: per-taxon replace.
"""

from collections.abc import Iterable

from refgenie.db.queries import species as queries
from refgenie.db.tables import Species, SpeciesName, Store, StoreType
from refgenie.exceptions import SpeciesFileError
from refgenie.logger import logger
from refgenie.managers.base import ResourceManager
from refgenie.utils.search import SearchOperator
from refgenie.utils.species_file import (
    SPECIES_FILE,
    SpeciesEntry,
    merge_species,
    read_species_file,
)


class SpeciesManager(ResourceManager):
    """Manager for the ``species`` and ``speciesname`` tables (species names from the stores)."""

    def load_from_stores(self, stores: Iterable[Store]) -> int | None:
        """Write the merge of these stores' species.json files, per taxon.

        Stores are read in priority order (lowest integer first, then name); see
        ``merge_species`` for the rule. A store with no file adds nothing. If any
        store's file cannot be read, nothing is written and the old rows stay.
        Each merged taxon is rewritten (``replace``); taxa no file lists, such as
        rows a pull wrote, are kept.

        Returns:
            The number of species written, or None if the old rows were kept.
        """
        per_store = []
        for store in sorted(stores, key=lambda s: (s.priority, s.name)):
            try:
                entries = read_species_file(store.url, on_disk=store.type == StoreType.on_disk)
            except SpeciesFileError as e:
                logger.warning(f"Kept the old species rows; store '{store.name}': {e}")
                return None
            if entries is None:
                logger.debug(f"Store '{store.name}' has no {SPECIES_FILE}")
                continue
            logger.debug(f"Read {len(entries)} species from store '{store.name}'")
            per_store.append(entries)
        return self.replace(merge_species(per_store))

    def replace(self, entries: Iterable[SpeciesEntry]) -> int:
        """Write ``entries`` per taxon, in one transaction: upsert each taxon's row and rewrite its names.

        Taxa not in ``entries`` are not touched. ``entries`` must have unique
        taxon ids (``merge_species`` guarantees it; a pull passes one entry).
        """
        rows = [
            Species(
                taxon_id=e.taxon_id,
                scientific_name=e.scientific_name,
                common_name=e.common_name,
                names=[SpeciesName(name=n) for n in e.names()],
            )
            for e in entries
        ]
        with self._database_session as session:
            queries.replace_taxa(session, rows)
            session.commit()
        return len(rows)

    def get(self, taxon_id: int) -> Species | None:
        """The species with this taxon id, or None."""
        return self.get_many([taxon_id]).get(taxon_id)

    def get_many(self, taxon_ids: Iterable[int]) -> dict[int, Species]:
        """The species rows for these taxon ids, by taxon id. Ids with no row are left out."""
        with self._database_session as session:
            return {s.taxon_id: s for s in queries.list_by_taxon_ids(session, taxon_ids)}

    def matching_taxon_ids(
        self, term: str, operator: SearchOperator = SearchOperator.CONTAINS
    ) -> list[int]:
        """Taxon ids with any name matching ``term``, ascending."""
        with self._database_session as session:
            return queries.matching_taxon_ids(session, term, operator)
