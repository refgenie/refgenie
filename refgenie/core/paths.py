"""A lazy, read-only view of local asset paths: ``paths[alias][group][seek_key]``.

:meth:`refgenie.core.root.Refgenie.paths` returns a fresh :class:`AssetPaths`.
It is keyed by genome alias, because pipeline configs name genomes by alias;
each alias resolves to its digest once. Nothing is looked up until it is read.
Each lookup goes through ``rg.asset.seek_components`` (the group's default
asset, then the named seek key), and
every answer, a miss included, is remembered for the life of the view. A new
view starts empty, so the caller decides how long answers stay valid (one
looper run, one request). Nothing is cached on disk: the database is the saved
lookup, and a disk copy would go stale after a pull or a remove.

A missing genome, group, asset or seek key raises :class:`KeyError`, which
Jinja treats as undefined, so ``{% if refgenie[g].blacklist is defined %}``
works. A group with no default asset counts as missing, as it does for
``refgenie seek``.

An asset group named like a ``Mapping`` method (``keys``, ``items``,
``values``, ``get``) must be read with brackets in Jinja
(``refgenie[g]["items"]``), because attribute access finds the method first.
"""

from collections.abc import Iterator, Mapping
from typing import TYPE_CHECKING

from refgenie.exceptions import (
    MissingAliasError,
    MissingAssetError,
    MissingAssetGroupError,
    MissingGenomeError,
    MissingSeekKeyError,
)
from refgenie.models import AssetRegistryPathComponents, GenomeAlias, GenomeDigest

if TYPE_CHECKING:  # pragma: no cover - annotation only; importing would cycle
    from refgenie.core.root import Refgenie

#: The lookups that mean "not here". Anything else (a DB error, a seek key with
#: no path) is a real failure and propagates.
_MISSING = (
    MissingAliasError,
    MissingGenomeError,
    MissingAssetGroupError,
    MissingAssetError,
    MissingSeekKeyError,
)


class AssetPaths(Mapping[str, "GenomePaths"]):
    """Read-only view: ``paths[alias][asset_group][seek_key] -> str``.

    Resolves on access through ``rg.asset.seek_components`` (default asset, then the
    named seek key) and remembers every answer for the life of the view.
    Missing genome, group, asset, or seek key raises KeyError, which Jinja
    treats as undefined, so ``{% if x is defined %}`` guards work.
    """

    def __init__(self, rg: "Refgenie"):
        self._rg = rg
        self._memo: dict[tuple, object] = {}

    def _remember(self, key: tuple, compute):
        if key not in self._memo:
            try:
                self._memo[key] = compute()
            except _MISSING:
                self._memo[key] = None
        return self._memo[key]

    def __getitem__(self, genome_alias: GenomeAlias) -> "GenomePaths":
        if not isinstance(genome_alias, str):
            raise KeyError(genome_alias)
        try:
            alias = GenomeAlias(genome_alias)
        except ValueError:
            raise KeyError(genome_alias) from None
        genome_digest = self._remember((alias,), lambda: self._rg.alias.resolve(alias))
        if genome_digest is None:
            raise KeyError(genome_alias)
        return GenomePaths(self, alias, genome_digest)

    def __iter__(self) -> Iterator[str]:
        return iter([alias.name for alias in self._rg.alias.list_all() if alias.name])

    def __len__(self) -> int:
        return sum(1 for _ in self)

    def to_dict(self) -> dict[str, dict[str, dict[str, str]]]:
        """Walk every genome, group and seek key into plain nested dicts.

        Leaves that do not resolve are left out, and so are groups and genomes
        left empty by that.
        """
        out: dict[str, dict[str, dict[str, str]]] = {}
        for genome in self:
            per_genome = {}
            for group in self[genome]:
                per_group = {}
                for seek_key in self[genome][group]:
                    try:
                        per_group[seek_key] = self[genome][group][seek_key]
                    except KeyError:
                        continue
                if per_group:
                    per_genome[group] = per_group
            if per_genome:
                out[genome] = per_genome
        return out


class GenomePaths(Mapping[str, "GroupPaths"]):
    """The asset groups of one genome. Built by :class:`AssetPaths`."""

    def __init__(self, root: AssetPaths, genome_alias: GenomeAlias, genome_digest: GenomeDigest):
        self._root = root
        self._alias = genome_alias
        self._digest = genome_digest

    def __getitem__(self, group: str) -> "GroupPaths":
        if not isinstance(group, str):
            raise KeyError(group)
        rg = self._root._rg
        # get_default raises MissingAssetGroupError for a missing group and
        # returns None for a group with no default; both are "not defined".
        default = self._root._remember(
            (self._alias, group),
            lambda: rg.asset.group.get_default(group, genome_digest=self._digest),
        )
        if default is None:
            raise KeyError(group)
        return GroupPaths(self._root, self._alias, self._digest, group)

    def __iter__(self) -> Iterator[str]:
        groups = self._root._rg.asset.group.list_all(genome_digests=[self._digest])
        return iter([g.name for g in groups if g.name in self])

    def __len__(self) -> int:
        return sum(1 for _ in self)


class GroupPaths(Mapping[str, str]):
    """The seek-key paths of one group's default asset. Built by :class:`GenomePaths`."""

    def __init__(
        self,
        root: AssetPaths,
        genome_alias: GenomeAlias,
        genome_digest: GenomeDigest,
        group: str,
    ):
        self._root = root
        self._alias = genome_alias
        self._digest = genome_digest
        self._group = group

    def __getitem__(self, seek_key: str) -> str:
        if not isinstance(seek_key, str):
            raise KeyError(seek_key)
        rg = self._root._rg
        value = self._root._remember(
            (self._alias, self._group, seek_key),
            lambda: str(
                rg.asset.seek_components(
                    AssetRegistryPathComponents(
                        genome=self._alias, asset_group=self._group, seek_key=seek_key
                    )
                )
            ),
        )
        if value is None:
            raise KeyError(seek_key)
        return value

    def __iter__(self) -> Iterator[str]:
        return iter(self._root._rg.asset.seek_key.list_all(self._digest, self._group))

    def __len__(self) -> int:
        return sum(1 for _ in self)
