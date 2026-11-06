"""Reading, validating and merging the stores' optional ``species.json`` files.

A store may carry a ``species.json`` next to its ``rgstore.json``::

    {"species": [{"taxon_id": 10116, "scientific_name": "Rattus norvegicus",
                  "common_name": "Norway rat", "synonyms": ["brown rat", "rat"]}]}

This is a leaf module: it never touches the database.
"""

import json
from collections.abc import Iterable
from pathlib import Path

from pydantic import BaseModel, ConfigDict, Field, ValidationError

from refgenie.exceptions import SpeciesFileError
from refgenie.logger import logger

SPECIES_FILE = "species.json"


def _unique(names: Iterable[str | None], exclude: Iterable[str] = ()) -> list[str]:
    """Non-blank names, stripped, each once ignoring case (first spelling kept).

    Names whose casefold is in ``exclude`` are dropped.
    """
    seen = {n.casefold() for n in exclude}
    out = []
    for name in names:
        if not isinstance(name, str):
            continue
        name = name.strip()
        if name and name.casefold() not in seen:
            seen.add(name.casefold())
            out.append(name)
    return out


class SpeciesEntry(BaseModel):
    """One entry of a store's species.json, as the registry writes it."""

    model_config = ConfigDict(strict=True, str_strip_whitespace=True, extra="ignore", frozen=True)

    taxon_id: int = Field(gt=0)
    scientific_name: str = Field(min_length=1)
    common_name: str | None = None
    synonyms: list[str] | None = None

    def names(self) -> list[str]:
        """Scientific name, common name, then synonyms: non-blank, each once ignoring case."""
        return _unique([self.scientific_name, self.common_name, *(self.synonyms or [])])


def parse_species(data: object, source: str) -> list[SpeciesEntry]:
    """The valid entries of a parsed species.json, in file order.

    Invalid entries and repeated taxon ids (the first wins) are skipped with a
    debug log.

    Raises:
        SpeciesFileError: If ``data`` is not an object with a ``species`` list.
    """
    if not isinstance(data, dict) or not isinstance(data.get("species"), list):
        raise SpeciesFileError(f"{source}: expected an object with a 'species' list")
    kept: dict[int, SpeciesEntry] = {}
    for i, item in enumerate(data["species"]):
        try:
            entry = SpeciesEntry.model_validate(item)
        except ValidationError as e:
            logger.debug(f"{source}: skipped species entry {i}: {e.errors()[0]['msg']}")
            continue
        if entry.taxon_id in kept:
            logger.debug(f"{source}: skipped species entry {i}: repeated taxon_id {entry.taxon_id}")
            continue
        kept[entry.taxon_id] = entry
    return list(kept.values())


def read_species_file(store_url: str, *, on_disk: bool) -> list[SpeciesEntry] | None:
    """Read and parse a store root's species.json.

    Returns:
        The entries, or None if the store has no file (or its root is a scheme
        that is not read).

    Raises:
        SpeciesFileError: If the file exists but cannot be read or has a bad top level.
    """
    if on_disk:
        path = Path(store_url) / SPECIES_FILE
        if not path.is_file():
            return None
        try:
            data = json.loads(path.read_text())
        except (OSError, ValueError) as e:
            raise SpeciesFileError(f"{path}: {e}") from e
        return parse_species(data, source=str(path))

    if not store_url.startswith(("http://", "https://")):
        logger.info(
            f"species.json not read from {store_url}: only on-disk and http(s) store roots are supported"
        )
        return None

    from refgenie.utils.http import make_client

    url = store_url.rstrip("/") + "/" + SPECIES_FILE
    try:
        with make_client() as client:
            response = client.get(url)
    except Exception as e:
        raise SpeciesFileError(f"{url}: {e}") from e
    # S3 answers 403 for a missing key.
    if response.status_code in (403, 404):
        return None
    if response.status_code != 200:
        raise SpeciesFileError(f"{url}: HTTP {response.status_code}")
    try:
        data = response.json()
    except ValueError as e:
        raise SpeciesFileError(f"{url}: {e}") from e
    return parse_species(data, source=url)


def merge_species(per_store: Iterable[list[SpeciesEntry]]) -> list[SpeciesEntry]:
    """Merge the stores' entries into one list, one entry per taxon id.

    ``per_store`` is in store-priority order, highest first. For each taxon:
    the scientific name is the first (highest-priority) store's; the common name
    is the first non-blank one in store order; the synonyms are every name any
    store gives the taxon, deduplicated ignoring case, minus the chosen
    scientific and common names. The result is sorted by taxon id.
    """
    grouped: dict[int, list[SpeciesEntry]] = {}
    for entries in per_store:
        for entry in entries:
            grouped.setdefault(entry.taxon_id, []).append(entry)
    merged = []
    for taxon_id in sorted(grouped):
        group = grouped[taxon_id]
        scientific = group[0].scientific_name
        common = next((c for e in group if (c := (e.common_name or "").strip())), None)
        synonyms = _unique(
            (n for e in group for n in e.names()),
            exclude=[scientific, common or ""],
        )
        merged.append(
            SpeciesEntry(
                taxon_id=taxon_id,
                scientific_name=scientific,
                common_name=common,
                synonyms=synonyms,
            )
        )
    return merged
