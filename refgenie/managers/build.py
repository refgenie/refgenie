"""
BuildManager - building assets from recipes with pypiper (``rgc.build``).
"""

from collections.abc import Callable
from pathlib import Path
from typing import TYPE_CHECKING, Any

from pypiper import PipelineManager
from pypiper.exceptions import SubprocessError
from pypiper.manager import COMPLETE_FLAG, FAIL_FLAG
from refgenie.db.engine import Engine

from refgenie.db.tables import Asset
from refgenie.exceptions import (
    BuildFailedError,
    ContainerEngineMissingError,
    MissingAliasError,
    MissingAssetError,
    MissingAssetGroupError,
    MissingBuildInputError,
    MissingRecipeError,
)
from refgenie.logger import logger
from refgenie.managers.asset.colocation import (
    create_colocation_symlinks,
    get_colocation_metadata,
)
from refgenie.managers import build_inputs, containers
from refgenie.managers.base import ResourceManager
from refgenie.models import (
    AssetRegistryPathComponents,
    BuildCommandValues,
    BuildParams,
    GenomeAlias,
    GenomeDigest,
)
from refgenie.plugins.events import NULL_EVENTS, EventSink, update_scope
from refgenie.plugins.hooks import POST_BUILD, PRE_BUILD, HookEvent
from refgenie.utils.io import coerce_cli_kwargs
from refgenie.utils.paths import get_build_dir
from refgenie.utils.signals import sigint_guard

if TYPE_CHECKING:
    from refgenie.managers.alias import AliasBackend
    from refgenie.managers.asset.manager import AssetManager
    from refgenie.managers.genome import GenomeManager
    from refgenie.managers.sources import RemoteGenomeSource
    from refgenie.managers.recipe import RecipeManager
    from refgenie.managers.sequence import SequenceManager
    from refgenie.managers.stage import StageManager


class BuildManager(ResourceManager):
    """
    The public build manager, at ``rgc.build``: building assets from recipes.

    Public surface:

    - ``run``: build one asset (optionally pulling the parent genome and
      staging the result).
    - ``preflight``: validate a prospective build without starting it.
    - ``initialize_and_build``: initialize a genome (from FASTA or a remote
      source) and build its ``fasta`` asset.
    - ``target_template`` / ``genome_init_target_template``: the Snakemake
      target paths for a built asset and for genome init.
    - ``resolve_custom_seek_keys``: a recipe's custom seek keys, cached, for
      naming an asset outside a build.
    - ``default_asset_name``: a recipe's default asset name, outside a build.
    - ``resolve_default_asset``: render a ``default_asset`` template.
    """

    def __init__(
        self,
        database_engine: Engine,
        recipe_manager: "RecipeManager",
        asset: "AssetManager",
        alias_manager: "AliasBackend",
        genome_manager: "GenomeManager",
        sequence_manager: "SequenceManager",
        stage_manager: "StageManager",
        stage_folder_getter: Callable[[], Path | None],
        pull_parent: Callable[[GenomeAlias], object],
        events: EventSink | None = None,
    ):
        """
        Initialize the BuildManager.

        Args:
            database_engine: The database engine.
            recipe_manager: The RecipeManager for getting recipes.
            asset: The AssetManager the built asset is recorded in. Its
                folders, ``content``, ``group``, ``tree`` and ``links`` are
                what a build writes through.
            alias_manager: The alias manager, for the genome being built and
                for input assets named by registry path.
            genome_manager: The GenomeManager, for ``initialize_and_build``.
            sequence_manager: Fetches a remote genome's sequences before a build reads them.
            stage_manager: The StageManager, for ``run(stage=True)``.
            stage_folder_getter: Returns the genome stage folder (or None) at
                call time.
            pull_parent: Pulls a missing parent genome, for
                ``run(pull_parents=True)``.
            events: Where ``pre_build`` / ``post_build`` are emitted for plugins.
        """
        super().__init__(database_engine)
        self._events = events or NULL_EVENTS
        self._recipe_manager = recipe_manager
        self._asset = asset
        self._alias_manager = alias_manager
        self._genome_manager = genome_manager
        self._sequence = sequence_manager
        self._stage_manager = stage_manager
        self._stage_folder_getter = stage_folder_getter
        self._pull_parent = pull_parent
        #: (recipe name, recipe version, container image) -> resolved custom seek keys.
        self._custom_seek_key_cache: dict[tuple[str, str | None, str | None], dict[str, Any]] = {}

    def resolve_custom_seek_keys(self, recipe: Any) -> dict[str, Any]:
        """Resolve a recipe's custom seek keys, to name an asset outside a build.

        Probes run in the recipe's ``docker_image`` when it declares one, else
        on the host. Results are cached per (recipe, version, image) because the
        build form preflights on every keystroke; failures are not cached.
        ``run`` does not use this: a build probes fresh at build time.

        Args:
            recipe: The recipe whose `custom_seek_keys` to resolve.

        Returns:
            dict[str, Any]: seek key name -> resolved value. Empty, with
            nothing run, for a recipe that declares no custom seek keys.

        Raises:
            CustomSeekKeyError: A probe failed, could not start, or timed out.
        """
        commands = recipe.custom_seek_keys or {}
        if not commands:
            return {}
        image = recipe.docker_image
        cache_key = (recipe.name, recipe.version, image)
        if cache_key not in self._custom_seek_key_cache:
            self._custom_seek_key_cache[cache_key] = build_inputs.probe_custom_seek_keys(
                commands, image, containers.detect_engine() if image else None
            )
        return dict(self._custom_seek_key_cache[cache_key])

    def default_asset_name(
        self, recipe: Any, command_values: BuildCommandValues | None = None
    ) -> str:
        """A recipe's default asset name, outside a build.

        Resolves the recipe's custom seek keys (cached) and renders
        ``default_asset``. Shared by ``preflight`` and the generated Snakefile.

        Args:
            recipe: The recipe to name an asset for.
            command_values: The namespace to render in. A minimal one is made
                when omitted. Its ``custom_seek_keys`` are overwritten.

        Returns:
            str: The resolved default asset name.

        Raises:
            CustomSeekKeyError: A probe failed, could not start, or timed out.
            ValueError: The template resolved to an empty name.
        """
        if command_values is None:
            command_values = BuildCommandValues(
                custom_seek_keys={},
                asset_group_name=recipe.name,
                genome_digest="",
                genome_folder=Path(""),
            )
        command_values.custom_seek_keys = self.resolve_custom_seek_keys(recipe) or {}
        return build_inputs.resolve_default_asset(recipe.default_asset, command_values)

    def _get_build_dir(
        self,
        genome_name: str,
        asset_group_name: str,
        asset_name: str,
    ) -> Path:
        """
        Get the build bookkeeping directory for a build invocation.

        Args:
            genome_name: The genome alias name (or a ``{genome_name}`` placeholder).
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.

        Returns:
            Path: The build directory.
        """
        return get_build_dir(
            genome_folder=self._asset.genome_folder,
            genome_name=genome_name,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )

    def _get_build_flag(
        self,
        genome_name: str,
        asset_group_name: str,
        asset_name: str,
    ) -> Path:
        """
        Get the path to the build-completion flag for a build invocation.

        Args:
            genome_name: The genome alias name (or a ``{genome_name}`` placeholder).
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.

        Returns:
            Path: The completion flag path.
        """
        build_dir = self._get_build_dir(genome_name, asset_group_name, asset_name)
        return build_dir / f"{genome_name}_{asset_group_name}__{asset_name}.flag"

    def genome_init_target_template(self) -> Path:
        """Return a path template for genome init sentinel files.

        The returned path contains ``{genome_name}`` for Snakemake wildcard
        substitution.  After ``refgenie1 genome init`` succeeds the sentinel
        file is touched so that downstream Snakemake rules can depend on it.
        """
        return self._asset.alias_folder / "{genome_name}" / ".genome_init_complete"

    def target_template(
        self,
        asset_group_name: str,
        asset_name: str,
    ) -> Path:
        """
        Get the build target template path.

        The ``{genome_name}`` placeholder is left unresolved for snakemake
        wildcard substitution. Substituted, it is the completion flag path.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.

        Returns:
            Path: The build target template path.
        """
        return self._get_build_flag("{genome_name}", asset_group_name, asset_name)

    def _resolve_input_assets(
        self,
        genome_digest: GenomeDigest,
        recipe_name: str,
        recipe_version: str | None = None,
        build_params: BuildParams | None = None,
    ) -> dict[str, Asset | None]:
        """
        Resolve input assets for a recipe based on user-provided input assets and defaults.

        Args:
            genome_digest: The digest of the genome.
            recipe_name: The name of the recipe.
            recipe_version: The version of the recipe.
            build_params: The build parameters.

        Returns:
            dict[str, Asset] | None: The resolved input assets.
        """
        if (
            inputs := self._recipe_manager.get(recipe_name, recipe_version).input_asset_classes
        ) is None:
            return None
        asset_by_input_asset_name: dict[str, Asset | None] = {i.name: None for i in inputs}
        if build_params is not None and build_params.assets is not None:
            # resolve user-provided input assets
            for asset_input_name, asset_registry_path in build_params.assets.items():
                if asset_input_name not in asset_by_input_asset_name:
                    raise ValueError(
                        f"Provided input asset is not required: {asset_input_name}. "
                        f"Required input assets are: {list(asset_by_input_asset_name.keys())}"
                    )
                parsed = AssetRegistryPathComponents.parse_registry_path(asset_registry_path)
                # The genome in a registry path is an alias; none means this genome.
                input_genome = (
                    self._alias_manager.resolve(GenomeAlias(parsed.genome))
                    if parsed.genome
                    else genome_digest
                )
                input_asset = self._asset.get(
                    genome_digest=input_genome,
                    asset_group_name=parsed.asset_group,
                    asset_name=parsed.asset
                    or self._asset.group.get_default(
                        parsed.asset_group, genome_digest=input_genome
                    ),
                )
                asset_by_input_asset_name[asset_input_name] = input_asset
        # fill the inputs the user did not provide with the recipe defaults
        for input_asset_name, default in asset_by_input_asset_name.items():
            if default is not None:
                continue
            for input in inputs:
                if input.name != input_asset_name:
                    continue
                input_asset = self._asset.get(
                    genome_digest=genome_digest,
                    asset_group_name=input.default,
                    asset_name=self._asset.group.get_default(
                        asset_group_name=input.default, genome_digest=genome_digest
                    )
                    or "default",
                )
                asset_by_input_asset_name[input_asset_name] = input_asset

        logger.debug(f"Resolved input assets: {asset_by_input_asset_name}")

        return asset_by_input_asset_name  # type: ignore

    @update_scope
    def run(
        self,
        recipe_name: str,
        genome_alias: GenomeAlias,
        asset_group_name: str,
        asset_name: str | None = None,
        recipe_version: str | None = None,
        params: BuildParams | None = None,
        stage: bool = False,
        container: bool = False,
        mounts: list[str] | None = None,
        asset_description: str | None = None,
        pull_parents: bool = False,
        pipeline_kwargs: dict[str, Any | None] = None,
        push_to: list[str] | None = None,
    ) -> Asset | None:
        """
        Build an asset for a specified genome based on provided recipe.

        Resolves the genome digest, runs the build (firing ``pre_build`` and
        ``post_build``), and stages the result if asked. The genome must already
        be initialized (via genome.initialize_genome or genome init CLI).

        ``post_build`` fires on every outcome, with ``succeeded`` True when an
        asset came back (built, or already present) and False when the build
        failed or raised. See :meth:`_build` for the build itself.

        Args:
            recipe_name: The name of the recipe to use
            genome_alias: The alias of the genome to build the asset for. The
                build folder is named after it.
            asset_group_name: The name of the asset group to use
            asset_name: The name of the asset to build
            recipe_version: The version of the recipe to use. If not provided, the latest version is used.
            params: The build parameters. Some recipes require input assets to be built.
            stage: Whether to stage the asset after building. Only possible if the stage folder is set.
            container: Run the recipe's commands inside its container image, if
                it declares one, with Docker or Apptainer.
            mounts: More host folders to make visible inside the container.
            asset_description: The description of the asset, used only if the asset does not exist.
            pull_parents: Whether to pull the parent genome if not found.
            pipeline_kwargs: Additional kwargs for PipelineManager.
            push_to: Optional remote names or ids to create push intent records for.
        Returns:
            Asset: The asset that was built (or the existing one, if skipped), or None if the build failed
        """
        # Check staging up front so a long build is not lost to a config error.
        if stage and self._stage_folder_getter() is None:
            raise ValueError(
                "Can't stage the asset: genome_stage_folder is not set. "
                "Set it in the refgenie config or build without --stage."
            )

        try:
            genome_digest = self._alias_manager.resolve(genome_alias)
        except MissingAliasError:
            msg = f"Genome '{genome_alias}' not found."
            if pull_parents:
                msg += " Attempting to pull from remote server."
                logger.warning(msg)
                self._pull_parent(GenomeAlias(genome_alias))
                genome_digest = self._alias_manager.resolve(genome_alias)
            else:
                msg += " Initialize it first with 'refgenie genome init'."
                logger.error(msg)
                raise

        self._events.emit(
            HookEvent(
                hook=PRE_BUILD, genome=genome_alias, asset_group=asset_group_name, asset=asset_name
            )
        )
        asset = None
        try:
            asset = self._build(
                recipe_name,
                genome_name=genome_alias,
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
                recipe_version=recipe_version,
                params=params,
                container=container,
                mounts=mounts,
                asset_description=asset_description,
                pipeline_kwargs=pipeline_kwargs,
            )
        finally:
            self._events.emit(
                HookEvent(
                    hook=POST_BUILD,
                    genome=genome_alias,
                    asset_group=asset_group_name,
                    asset=getattr(asset, "name", None) or asset_name,
                    succeeded=asset is not None,
                )
            )

        if stage and asset is not None:
            self._stage_manager.create(
                asset=asset,
                genome_folder=self._asset.genome_folder,
                genome_stage_folder=self._stage_folder_getter(),
                push_to=push_to,
                build_dir=self._get_build_dir(genome_alias, asset_group_name, asset.name),
            )

        return asset

    def preflight(
        self,
        recipe_name: str,
        genome_alias: GenomeAlias,
        asset_group_name: str,
        asset_name: str | None = None,
        recipe_version: str | None = None,
        params: BuildParams | None = None,
        stage: bool = False,
    ) -> dict[str, Any]:
        """Validate a prospective build without starting it.

        Runs the same resolution and validation as `run` but starts no
        pipeline and writes nothing. Used by the web build form.

        Args:
            recipe_name: The name of the recipe to use.
            genome_alias: The alias of the genome to build for.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset to build. If omitted, the
                recipe's ``default_asset`` template is resolved.
            recipe_version: The version of the recipe. Latest if omitted.
            params: The build parameters (assets/params/files).
            stage: Whether the build would stage the asset afterwards.

        Returns:
            ``{"ok": bool, "errors": [{"field", "code", "message"}, ...],
            "resolved": {...}}`` -- ``errors`` is field-scoped so a form can
            attach each problem to its input; ``resolved`` reports what the
            build would use (genome digest, recipe version, input assets,
            asset name).
        """
        errors: list[dict[str, str]] = []
        resolved: dict[str, Any] = {}

        def problem(field: str, code: str, message: str) -> None:
            errors.append({"field": field, "code": code, "message": message})

        if stage and self._stage_folder_getter() is None:
            problem(
                "stage",
                "missing_build_input",
                "Staging is not configured: genome_stage_folder is not set. "
                "Set it in the refgenie config or build without staging.",
            )

        genome_digest: GenomeDigest | None = None
        try:
            genome_digest = self._alias_manager.resolve(genome_alias)
            resolved["genome_digest"] = genome_digest
        except MissingAliasError:
            problem(
                "genome",
                "genome_not_found",
                f"Genome '{genome_alias}' not found. Initialize it first with "
                "'refgenie genome init', or build with pull_parents.",
            )

        recipe = None
        try:
            recipe = self._recipe_manager.get(recipe_name, recipe_version)
            resolved["recipe"] = recipe_name
            resolved["recipe_version"] = recipe.version
        except MissingRecipeError as exc:
            problem("recipe", "recipe_not_found", str(exc))

        if genome_digest is None or recipe is None:
            return {"ok": False, "errors": errors, "resolved": resolved}

        engine_missing = bool(recipe.docker_image) and containers.detect_engine() is None
        if engine_missing:
            problem(
                "recipe",
                "container_engine_missing",
                str(ContainerEngineMissingError(recipe.docker_image)),
            )

        # The steps of `_prepare`, one at a time: a build stops at the first
        # problem, a preflight reports each under its own field.
        if params is not None:
            params.populate_with_defaults_from_recipe(recipe)

        input_assets = None
        try:
            input_assets = self._resolve_input_assets(
                genome_digest=genome_digest,
                recipe_name=recipe.name,
                recipe_version=recipe.version,
                build_params=params,
            )
            resolved["input_assets"] = {
                name: (asset.registry_path if asset is not None else None)
                for name, asset in (input_assets or {}).items()
            }
        except (MissingAssetError, MissingAssetGroupError, MissingAliasError) as exc:
            problem("params.assets", "asset_not_found", str(exc))
        except ValueError as exc:
            problem("params.assets", "conflict", str(exc))

        command_values = self._build_command_values(
            genome_digest, asset_group_name, params, input_assets, {}
        )
        try:
            build_inputs.validate_build_inputs(recipe, command_values)
        except MissingBuildInputError as exc:
            problem("params", "missing_build_input", str(exc))

        if asset_name:
            resolved["asset_name"] = asset_name
        elif engine_missing and recipe.custom_seek_keys:
            # Naming needs the image too; the recipe problem above says why.
            pass
        else:
            # Its own code, so a form asks for a name instead of showing an error.
            try:
                resolved["asset_name"] = self.default_asset_name(recipe, command_values)
            except Exception as exc:  # noqa: BLE001 - report, don't crash a preflight
                problem("asset", "asset_name_required", str(exc))

        return {"ok": not errors, "errors": errors, "resolved": resolved}

    @update_scope
    def initialize_and_build(
        self,
        genome_names: list[GenomeAlias],
        fasta_file_path: Path | None = None,
        source: "RemoteGenomeSource | None" = None,
        digest: GenomeDigest | None = None,
        description: str = "",
        species_name: str | None = None,
        use_existing: bool = False,
        build_fasta: bool = True,
    ) -> tuple[GenomeDigest, bool]:
        """Initialize a genome and build its fasta asset.

        This is the recommended way to add a new genome. Give either a local
        FASTA or a remote ``source`` (with ``digest``); see
        ``GenomeManager.initialize_genome``. The fasta asset (fa, fai,
        chrom.sizes) is built only for a newly created genome, and only when a
        ``fasta`` recipe is registered.

        Args:
            genome_names: Alias names for the genome.
            fasta_file_path: Path to a local FASTA file.
            source: A remote genome source, in place of a FASTA.
            digest: The genome's digest on ``source``.
            description: Genome description.
            species_name: Species name.
            use_existing: Allow re-init of existing genome.
            build_fasta: Whether to build the fasta asset (default True).

        Returns:
            Tuple of (genome_digest, was_created).

        Raises:
            BuildFailedError: The genome was registered but its fasta build failed.
        """
        digest, created = self._genome_manager.initialize_genome(
            fasta_file_path=fasta_file_path,
            source=source,
            digest=digest,
            description=description,
            alias_names=genome_names,
            use_existing=use_existing,
            species_name=species_name,
        )
        if not (build_fasta and created):
            return digest, created
        if not self._recipe_manager.exists("fasta"):
            logger.info("Skipping automatic fasta build: no 'fasta' recipe is registered.")
            return digest, created

        alias = genome_names[0]
        if self.run(recipe_name="fasta", genome_alias=alias, asset_group_name="fasta") is None:
            raise BuildFailedError(
                f"Genome '{alias}' ({digest}) was registered, but its fasta asset failed to build."
            )
        return digest, created

    def _build_command_values(
        self,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        params: BuildParams | None,
        input_assets: dict[str, Asset | None] | None,
        custom_seek_keys: dict[str, Any],
    ) -> BuildCommandValues:
        """
        The values a recipe's command templates are rendered with."""
        return BuildCommandValues(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            params=params.params if params is not None else {},
            files=params.files if params is not None else {},
            assets=input_assets,
            custom_seek_keys=custom_seek_keys,
            genome_folder=self._asset.genome_folder,
            refget_store_path=str(self._asset.genome_folder / ".refget_store"),
        )

    def _prepare(
        self,
        recipe: Any,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        params: BuildParams | None,
        custom_seek_keys: dict[str, Any],
    ) -> tuple[dict[str, Asset | None] | None, BuildCommandValues]:
        """
        Everything a build needs before it runs anything, checked.

        Fills ``params`` with the recipe's defaults, resolves the input assets
        and validates the required inputs. ``preflight`` runs the same steps
        one by one so it can report every problem, not just the first.

        Args:
            recipe: The recipe to build with.
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            params: The build parameters. Filled with recipe defaults in place.
            custom_seek_keys: The resolved custom seek keys.

        Returns:
            (input assets by input name, the command values).

        Raises:
            MissingAssetError, MissingAssetGroupError, MissingAliasError,
                ValueError: From resolving the input assets.
            MissingBuildInputError: If a required file or parameter is missing.
        """
        if params is not None:
            params.populate_with_defaults_from_recipe(recipe)
        input_assets = self._resolve_input_assets(
            genome_digest=genome_digest,
            recipe_name=recipe.name,
            recipe_version=recipe.version,
            build_params=params,
        )
        command_values = self._build_command_values(
            genome_digest, asset_group_name, params, input_assets, custom_seek_keys
        )
        build_inputs.validate_build_inputs(recipe, command_values)
        return input_assets, command_values

    def _build(
        self,
        recipe_name: str,
        *,
        genome_name: GenomeAlias,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str | None = None,
        recipe_version: str | None = None,
        params: BuildParams | None = None,
        container: bool = False,
        mounts: list[str] | None = None,
        asset_description: str | None = None,
        pipeline_kwargs: dict[str, Any | None] = None,
    ) -> Asset | None:
        """
        Build an asset using a recipe: prepare, run the pipeline, record.

        Args:
            recipe_name: The name of the recipe to use.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset to build.
            recipe_version: The version of the recipe.
            params: The build parameters.
            container: Run in the recipe's container image, if it has one.
            mounts: More host folders to make visible inside the container.
            asset_description: Description of the asset.
            pipeline_kwargs: Additional kwargs for PipelineManager.

        Returns:
            Asset | None: The built asset (or the existing asset when the
                build is skipped), or None if the build failed.

        Raises:
            ContainerEngineMissingError: The build is to run in the recipe's
                image, and this host has neither Docker nor Apptainer.
        """
        recipe = self._recipe_manager.get(recipe_name, recipe_version)
        logger.info(f"Building '{genome_name}/{asset_group_name}' using recipe '{recipe}'")

        image = recipe.docker_image if container else None
        engine = containers.detect_engine() if image else None
        if image and engine is None:
            raise ContainerEngineMissingError(image)
        custom_seek_keys = (
            build_inputs.probe_custom_seek_keys(recipe.custom_seek_keys or {}, image, engine) or {}
        )
        logger.debug(f"Build parameters: {params}")
        logger.debug(f"Custom seek keys: {custom_seek_keys}")
        input_assets, command_values = self._prepare(
            recipe, genome_digest, asset_group_name, params, custom_seek_keys
        )
        # Not `default_asset_name`, which uses the cached probe.
        asset_name = asset_name or build_inputs.resolve_default_asset(
            recipe.default_asset, command_values
        )

        if self._asset.exists(asset_group_name, asset_name, genome_digest=genome_digest):
            logger.warning(
                f"Asset '{genome_name}/{asset_group_name}:{asset_name}' already exists. "
                f"Skipping build"
            )
            # Snakemake still expects the flag; it may have been pruned from builds/.
            flag = self._get_build_flag(genome_name, asset_group_name, asset_name)
            if not flag.exists():
                logger.info(f"Restoring missing build flag for skipped build: {flag}")
                flag.parent.mkdir(parents=True, exist_ok=True)
                flag.touch()
            return self._asset.get(asset_group_name, asset_name, genome_digest=genome_digest)

        if any("refget_store_path" in t for t in recipe.command_templates):
            self._sequence.fetch(genome_digest)

        build_output_folder = (
            self._asset.data_folder / genome_digest / asset_group_name / asset_name
        )
        command_values.output_folder = build_output_folder

        if not self._run_pipeline(
            recipe,
            command_values,
            input_assets,
            genome_name=genome_name,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            image=image,
            engine=engine,
            mounts=mounts,
            pipeline_kwargs=pipeline_kwargs,
        ):
            return None

        return self._record(
            recipe,
            command_values,
            input_assets,
            genome_name=genome_name,
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            description=asset_description or recipe.description,
            image=image,
            engine=engine,
        )

    def _run_pipeline(
        self,
        recipe: Any,
        command_values: BuildCommandValues,
        input_assets: dict[str, Asset | None] | None,
        *,
        genome_name: GenomeAlias,
        asset_group_name: str,
        asset_name: str,
        image: str | None,
        engine: containers.ContainerEngine | None,
        mounts: list[str] | None,
        pipeline_kwargs: dict[str, Any | None] | None,
    ) -> bool:
        """
        Run the recipe's commands through pypiper, into ``command_values.output_folder``.

        With an ``image``, each command is wrapped to run inside it under
        ``engine`` (`containers.wrap`), with the output folder, the genome
        folder and ``mounts`` visible at the same paths inside.

        Returns:
            bool: Whether the pipeline succeeded.
        """
        build_target_string = f"{genome_name}/{asset_group_name}:{asset_name}"
        build_output_folder = command_values.output_folder

        # Before the commands run: tools like BWA need the parent file (.fa)
        # present in the output folder.
        create_colocation_symlinks(
            output_folder=build_output_folder,
            genome_folder=self._asset.genome_folder,
            input_assets=recipe.input_assets,
            resolved_assets=input_assets,
        )

        commands = build_inputs.populate_commands(
            command_templates=recipe.command_templates,
            command_values=command_values,
        )
        build_stats_output_folder = self._get_build_dir(genome_name, asset_group_name, asset_name)
        build_flag = self._get_build_flag(genome_name, asset_group_name, asset_name)
        target = build_flag.as_posix()
        # The flag means "this asset is in the catalog", so _record writes it
        # after content.add commits, never a recipe command.
        logger.debug(f"Commands to be executed: {commands}")

        # The asset is not in the catalog, so a flag here is from an interrupted
        # build; left in place, pypiper would skip every command.
        if build_flag.exists():
            logger.warning(
                f"Removing a stale build flag for an asset that is not in the catalog: {build_flag}"
            )
            build_flag.unlink()

        coerced_kwargs = coerce_cli_kwargs(pipeline_kwargs) if pipeline_kwargs else {}
        pm = PipelineManager(
            name=f"refgenie_{genome_name}_{asset_group_name}_{asset_name}",
            outfolder=build_stats_output_folder,
            recover=True,
            **coerced_kwargs,
        )
        # None lets pypiper split pipes itself.
        shell = None
        if image:
            logger.info(f"Running the recipe in {image} with {engine}")
            # The output folder must exist to be mounted: Apptainer refuses a
            # missing bind source, and Docker would create it owned by root.
            build_output_folder.mkdir(parents=True, exist_ok=True)
            # Mount the genome folder too: colocation symlinks are relative and
            # point into sibling asset folders.
            all_mounts = [
                *(mounts or []),
                build_output_folder.as_posix(),
                self._asset.genome_folder.as_posix(),
            ]
            commands = [containers.wrap(engine, image, c, mounts=all_mounts) for c in commands]
            # Each wrapped command is one engine call whose quoted argument may
            # hold pipes; pypiper splitting it on `|` would cut it apart.
            shell = True

        return_code = None
        try:
            failed = False
            with sigint_guard(
                f"The build was interrupted for {genome_name}/{asset_group_name}:{asset_name}"
            ):
                if commands:
                    return_code = pm.run(commands, target, shell=shell)
                else:
                    # e.g. output is colocation symlinks alone
                    return_code = 0
        except SubprocessError:
            failed = True

        if return_code is None or return_code > 0:
            logger.error(f"Asset '{build_target_string}' build failed. Exit code: {return_code}")
            failed = True
        else:
            logger.info(f"Asset '{build_target_string}' build succeeded")
        # stop_pipeline() defaults to 'completed'; pass the real status or a
        # failed build records success in pipestat and the on-disk flags.
        pm.stop_pipeline(status=FAIL_FLAG if failed else COMPLETE_FLAG)
        return not failed

    def _record(
        self,
        recipe: Any,
        command_values: BuildCommandValues,
        input_assets: dict[str, Asset | None] | None,
        *,
        genome_name: GenomeAlias,
        genome_digest: GenomeDigest,
        asset_group_name: str,
        asset_name: str,
        description: str | None,
        image: str | None,
        engine: containers.ContainerEngine | None,
    ) -> Asset | None:
        """
        Record a finished build: the asset, its completion flag, parents and
        provenance. The alias tree is rendered inside ``content.add``.

        Returns:
            Asset | None: The added asset.
        """
        input_asset_digests = {
            name: asset.digest
            for name, asset in (input_assets or {}).items()
            if asset is not None and asset.digest
        }
        build_provenance = build_inputs.build_provenance(
            command_values=command_values,
            recipe=recipe,
            genome_digest=genome_digest,
            input_asset_digests=input_asset_digests,
            input_file_digests=build_inputs.input_file_digests(command_values),
            image=image,
            engine=engine,
        )

        added_asset = self._asset.content.add(
            asset_class_name=recipe.output_asset_class.name,
            path=command_values.output_folder.relative_to(self._asset.genome_folder),
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            description=description,
            recipe=recipe,
            custom_seek_keys=command_values.custom_seek_keys,
            colocate=get_colocation_metadata(recipe.input_assets),
            set_default=True,  # Deliberate builds become the group default
            build_provenance=build_provenance,
        )
        logger.info(f"Added asset: '{added_asset}'")
        if added_asset is not None:
            build_flag = self._get_build_flag(genome_name, asset_group_name, asset_name)
            build_flag.parent.mkdir(parents=True, exist_ok=True)
            build_flag.touch()

        if input_asset_digests:
            self._asset.links.set_parents(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
                parent_asset_digests=list(input_asset_digests.values()),
            )

        return added_asset
