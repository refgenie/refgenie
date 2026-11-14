"""
Wire models for the MCP tools.

List tools return pages of slim rows, since an agent reads every token; the
``get_*`` tools return one full record.
"""

from typing import Any

from pydantic import BaseModel

from refgenie.mcp.paging import ResultPage
from refgenie.db.tables import Asset, AssetClass, Recipe, Species
from refgenie.models import GenomeListing


def _summary(description: str | None) -> str | None:
    """The first line of a description, cut to 200 characters."""
    if not description:
        return None
    line = description.strip().splitlines()[0]
    return line if len(line) <= 200 else line[:197] + "..."


class SpeciesNames(BaseModel):
    """A genome's species from the species list.

    Separate from the HTTP ``SpeciesResponse``: the stdio install may lack the
    server extras that ``refgenie.server.schemas`` imports.
    """

    taxon_id: int
    scientific_name: str
    common_name: str | None
    #: The species' other common names, e.g. "brown rat".
    synonyms: list[str]

    @classmethod
    def of(cls, species: Species | None) -> "SpeciesNames | None":
        if species is None:
            return None
        return cls(
            taxon_id=species.taxon_id,
            scientific_name=species.scientific_name,
            common_name=species.common_name,
            synonyms=species.synonyms,
        )


class GenomeRow(BaseModel):
    digest: str
    aliases: list[str]
    species_name: str | None
    #: The species list's common name for the genome's taxon.
    common_name: str | None
    asset_count: int

    @classmethod
    def of(cls, listing: GenomeListing) -> "GenomeRow":
        genome = listing.genome
        return cls(
            digest=genome.digest,
            aliases=listing.aliases,
            species_name=genome.species_name,
            common_name=listing.common_name,
            asset_count=listing.asset_count,
        )


class GenomeDetail(GenomeRow):
    description: str | None
    #: Other names of the genome itself, from FHR ``genomeSynonym`` (e.g. GRCh38).
    synonyms: list[str] | None
    taxon_id: int | None
    assembly_source: str | None
    assembly_accession: str | None
    assembly_level: str | None
    #: The FHR metadata record, if the genome's store has one.
    fhr: dict[str, Any] | None = None
    #: The genome's species and all its names; None if the species list has no row.
    species: SpeciesNames | None = None


class AssetRow(BaseModel):
    digest: str
    asset_class: str | None
    #: The asset's canonical name in its group (``asset_class`` names the group).
    name: str
    genome_digest: str | None
    size: int | None
    serving_modes: list[str]

    @classmethod
    def of(cls, asset: Asset) -> "AssetRow":
        return cls(
            digest=asset.digest,
            asset_class=asset.asset_class_name,
            name=asset.name,
            genome_digest=asset.genome_digest,
            size=asset.size,
            serving_modes=asset.serving_modes,
        )


class AssetDetail(AssetRow):
    asset_group: str | None
    #: Every name this content has in its group.
    names: list[str]
    description: str | None
    #: Seek key name -> raw value, relative to the asset's folder.
    seek_keys: dict[str, str]
    parents: list[AssetRow]
    children: list[AssetRow]
    #: Server mode, archive-served assets: the URL of the .tgz archive.
    archive_url: str | None = None
    #: Path seek key -> where its file is: a local path, or a URL on a server
    #: (file-served assets only).
    files: dict[str, str]
    #: The command that pulls this asset.
    pull_command: str


class RecipeRow(BaseModel):
    name: str
    version: str
    output_asset_class: str | None
    summary: str | None
    #: Input asset classes, then input files marked ``(file)``.
    needs: list[str]

    @classmethod
    def of(cls, recipe: Recipe) -> "RecipeRow":
        return cls(
            name=recipe.name,
            version=recipe.version,
            output_asset_class=recipe.output_asset_class.name
            if recipe.output_asset_class
            else None,
            summary=_summary(recipe.description),
            needs=[link.asset_class.name for link in recipe.input_asset_classes or []]
            + [f"{name} (file)" for name in recipe.input_files or {}],
        )


class RecipeDetail(BaseModel):
    name: str
    version: str
    description: str | None
    output_asset_class: str | None
    input_assets: dict[str, Any] | None
    input_files: dict[str, Any] | None
    input_params: dict[str, Any] | None
    command_templates: list[str]
    docker_image: str | None
    custom_seek_keys: dict[str, str] | None
    default_asset: str
    #: The command that builds an asset with this recipe.
    next: str

    @classmethod
    def of(cls, recipe: Recipe) -> "RecipeDetail":
        output = recipe.output_asset_class.name if recipe.output_asset_class else None
        return cls(
            name=recipe.name,
            version=recipe.version,
            description=recipe.description,
            output_asset_class=output,
            input_assets=recipe.input_assets,
            input_files=recipe.input_files,
            input_params=recipe.input_params,
            command_templates=recipe.command_templates,
            docker_image=recipe.docker_image,
            custom_seek_keys=recipe.custom_seek_keys,
            default_asset=recipe.default_asset,
            next=f"refgenie build <genome>/{output or recipe.name}",
        )


class AssetClassRow(BaseModel):
    name: str
    version: str
    summary: str | None
    seek_keys: list[str]

    @classmethod
    def of(cls, asset_class: AssetClass) -> "AssetClassRow":
        return cls(
            name=asset_class.name,
            version=asset_class.version,
            summary=_summary(asset_class.description),
            seek_keys=[seek_key.name for seek_key in asset_class.seek_keys],
        )


class DigestLookup(BaseModel):
    #: ``genome``, ``asset`` or ``not_found``.
    type: str
    digest: str
    genome: GenomeRow | None = None
    asset: AssetRow | None = None


class SpeciesCount(BaseModel):
    species: str
    genomes: int
    asset_classes: int
    assets: int


class CatalogSummary(BaseModel):
    #: Which server this is: the local database, or the server's URL.
    server: str
    genomes: int
    asset_groups: int
    assets: int
    #: Species with at least one asset, most assets first.
    species_with_assets: list[SpeciesCount]


class SequenceRow(BaseModel):
    name: str
    length: int
    #: The refget sequence id (``SQ.``...).
    digest: str


class SequencePage(ResultPage[SequenceRow]):
    """A page of a genome's sequences, with totals for the whole genome."""

    n_sequences: int
    total_length: int


class GenomeComparison(BaseModel):
    #: One sentence: whether the genomes are the same and how they differ.
    summary: str
    #: A's sequence name -> B's name, for sequences that share content but not a name.
    name_map: dict[str, str]
    #: The GA4GH seqcol comparison object.
    comparison: dict[str, Any]
