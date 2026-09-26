"""The `pull` and `mirror` commands: models and handlers.

`LargeArchiveFlags` lives here too -- these are the only two commands that
take the large-archive prompt flags.
"""

from pydantic import AliasChoices, BaseModel, ConfigDict, Field, model_validator

from refgenie.cli.framework import CliList
from refgenie.cli.commands.helpers import (
    _ASSET_REGISTRY_PATHS_DESCRIPTION,
    _GENOME_DESCRIPTION,
    _GENOME_DIGEST_DESCRIPTION,
    _validate_registry_paths,
    data_channel_hint,
    genome_ref,
    no_data_channel_synced,
)
from refgenie.cli.errors import (
    EXIT_GENERAL_ERROR,
    EXIT_INVALID_INPUT,
    EXIT_NOT_FOUND,
    EXIT_OK,
    fail,
    run_over_paths,
)
from refgenie.const import DEFAULT_PULL_SIZE_CUTOFF_GB
from refgenie.logger import logger


class LargeArchiveFlags(BaseModel):
    """Shared prompt-handling flags for pull and mirror commands.

    Note: The flag is --skip-large, not --no-large: cli_implicit_flags=True
    auto-generates --no-* negations for every boolean, so --no-large would be
    ambiguous (it could mean "set no_large=True" or "negate large=False").
    """

    model_config = ConfigDict(populate_by_name=True)

    skip_large: bool = Field(
        default=False,
        description="Do not pull archives over the size cutoff.",
        validation_alias=AliasChoices("skip-large"),
    )
    pull_large: bool = Field(
        default=False,
        description="Pull all archives regardless of size.",
        validation_alias=AliasChoices("pull-large"),
    )
    size_cutoff: float = Field(
        # Single source of truth: refgenie.const.DEFAULT_PULL_SIZE_CUTOFF_GB.
        # This class is the only definition of the CLI-facing default -- do not
        # restate the number anywhere else.
        default=DEFAULT_PULL_SIZE_CUTOFF_GB,
        description="Maximum archive size to pull without confirmation, in GB.",
        validation_alias=AliasChoices("size-cutoff"),
    )
    batch: bool = Field(
        default=False,
        description="Batch mode: pull all archives regardless of size.",
        validation_alias=AliasChoices("batch"),
    )

    @model_validator(mode="after")
    def check_large_exclusion(self):
        if self.skip_large and self.pull_large:
            raise ValueError("--skip-large and --pull-large are mutually exclusive")
        return self

    @model_validator(mode="after")
    def apply_batch_mode(self):
        """--batch is a convenience flag that sets pull_large=True.

        If --batch is combined with --skip-large, raise an error
        since the user's intent is contradictory.
        """
        if self.batch:
            if self.skip_large:
                raise ValueError("--batch (implies pull-large) conflicts with --skip-large")
            self.pull_large = True
        return self

    def resolve_force_large(self) -> bool | None:
        """Resolve the large flag to a tri-state value for the library layer.

        Returns False if --skip-large, True if --pull-large, None if neither.
        """
        if self.skip_large:
            return False
        if self.pull_large:
            return True
        return None


class PullModel(LargeArchiveFlags):
    # `asset_registry_paths` is a named (non-positional) field because
    # pydantic-settings rejects positional args with defaults; argv
    # preprocessing in `refgenie.cli.main` rewrites bare positional paths to
    # `--asset-registry-paths <paths>`.
    asset_registry_paths: CliList | None = Field(
        default=None,
        description=_ASSET_REGISTRY_PATHS_DESCRIPTION,
        validation_alias=AliasChoices("asset-registry-paths"),
    )
    genome: str | None = Field(
        default=None,
        description=_GENOME_DESCRIPTION + " Comma-separate several.",
        validation_alias=AliasChoices("g", "genome"),
    )
    genome_digest: str | None = Field(
        default=None,
        description=_GENOME_DIGEST_DESCRIPTION + " Comma-separate several for --all or --asset.",
        validation_alias=AliasChoices("genome-digest"),
    )
    all: bool = Field(
        default=False,
        description="Pull all assets for the specified genome(s).",
        validation_alias=AliasChoices("all"),
    )
    all_genomes: bool = Field(
        default=False,
        description="Apply the operation to all genomes available on subscribed servers.",
        validation_alias=AliasChoices("all-genomes"),
    )
    asset: str | None = Field(
        default=None,
        description="Pull a specific asset type across the specified genomes (e.g. --asset fasta).",
        validation_alias=AliasChoices("asset"),
    )
    init: bool = Field(
        default=False,
        description="Register genome(s) in the local store (aliases, metadata) "
        "without downloading asset files.",
        validation_alias=AliasChoices("init"),
    )
    force: bool = Field(
        default=False,
        description="Skip confirmation prompts for multi-asset operations.",
        validation_alias=AliasChoices("f", "force"),
    )


class MirrorModel(LargeArchiveFlags):
    force: bool = Field(
        default=False,
        description="Skip confirmation prompt.",
        validation_alias=AliasChoices("f", "force"),
    )


def _split(value: str | None) -> list[str]:
    return [v.strip() for v in (value or "").split(",") if v.strip()]


def handle_pull(cmd, refgenie) -> None:
    from refgenie.models import GenomeAlias

    genome_names = _split(cmd.genome)
    genome_digests = _split(cmd.genome_digest)
    genomes = [genome_ref(g, None) for g in genome_names] + [
        genome_ref(None, d) for d in genome_digests
    ]

    asset_registry_paths = cmd.asset_registry_paths
    force_large = cmd.resolve_force_large()

    if getattr(cmd, "init", False):
        if genome_digests:
            fail("--init registers genomes by alias; use -g, not --genome-digest.")
        registered = refgenie.transfer.init_genomes(
            aliases=[GenomeAlias(g) for g in genome_names] or None,
            all_genomes=getattr(cmd, "all_genomes", False),
        )
        if not registered:
            fail("Init failed: no genomes were registered.", EXIT_NOT_FOUND)
    elif getattr(cmd, "all", False):
        if not genomes and not getattr(cmd, "all_genomes", False):
            fail("Specify genome(s) with -g or --genome-digest, or use --all-genomes")
        if getattr(cmd, "all_genomes", False):
            fail("Use 'refgenie mirror' to pull all assets for all genomes.")
        pulled = refgenie.transfer.pull_genomes(
            genomes=genomes,
            force=cmd.force or None,
            force_large=force_large,
            size_cutoff=cmd.size_cutoff,
        )
        if not pulled:
            fail(f"Pull failed: no assets were pulled for {', '.join(genomes)}.", EXIT_NOT_FOUND)
    elif getattr(cmd, "asset", None):
        if not genomes and not getattr(cmd, "all_genomes", False):
            fail("Specify genome(s) with -g or --genome-digest, or use --all-genomes")
        pulled = refgenie.transfer.pull_genomes(
            genomes=genomes or None,
            asset_group_name=cmd.asset,
            all_genomes=getattr(cmd, "all_genomes", False),
            force=cmd.force or None,
            force_large=force_large,
            size_cutoff=cmd.size_cutoff,
        )
        if not pulled:
            fail(f"Pull failed: no '{cmd.asset}' assets were pulled.", EXIT_NOT_FOUND)
    elif asset_registry_paths:
        from refgenie.exceptions import (
            AssetExistsError,
            MissingAssetClassError,
            MissingRecipeError,
            PullFailedError,
        )

        parsed_asset_registry_paths = [
            refgenie.parse_asset_registry_path(path) for path in asset_registry_paths
        ]
        if len(genome_digests) > 1:
            fail("Pulling registry paths takes one --genome-digest.", EXIT_INVALID_INPUT)
        genome_digest = genome_digests[0] if genome_digests else None
        _validate_registry_paths(parsed_asset_registry_paths, asset_registry_paths, genome_digest)

        def pull_one(parsed) -> int:
            genome = genome_ref(parsed.genome, genome_digest)
            try:
                pulled = refgenie.transfer.pull(
                    genome=genome,
                    asset_group_name=parsed.asset_group,
                    asset_name=parsed.asset,
                    force=cmd.force or None,
                    force_large=force_large,
                    size_cutoff=cmd.size_cutoff,
                )
            except AssetExistsError as e:
                logger.error(f"{e}. Use --force to overwrite.")
                return EXIT_GENERAL_ERROR
            except PullFailedError as e:
                logger.error(str(e))
                return EXIT_NOT_FOUND
            except (MissingAssetClassError, MissingRecipeError) as e:
                logger.error(str(e))
                if no_data_channel_synced(refgenie):
                    retry = (
                        f"refgenie pull {parsed.asset_group} --genome-digest {genome}"
                        if genome_digest
                        else f"refgenie pull {genome}/{parsed.asset_group}"
                    )
                    logger.info(data_channel_hint(retry))
                return EXIT_NOT_FOUND
            if pulled is None:
                # transfer.pull signals failure by RETURNING None (no server could
                # serve the asset, or the pull was skipped), not only by raising.
                # Checking only for exceptions reports a failed pull as a success.
                logger.error(
                    f"Pull failed: no asset was retrieved for '{genome}/{parsed.asset_group}'."
                )
                return EXIT_NOT_FOUND
            return EXIT_OK

        run_over_paths(parsed_asset_registry_paths, pull_one)
    else:
        fail("Specify assets to pull, or use --all, --asset, or --init")


def handle_mirror(cmd, refgenie) -> None:
    refgenie.transfer.mirror(
        force=cmd.force or None,
        force_large=cmd.resolve_force_large(),
        size_cutoff=cmd.size_cutoff,
    )
