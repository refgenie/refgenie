"""
AssetBuilder - internal service for asset construction using recipes and pypiper.
"""

import signal
import threading
from datetime import datetime, timezone
from pathlib import Path
from subprocess import CalledProcessError, check_output
from typing import Any, TYPE_CHECKING

from pypiper import PipelineManager
from pypiper.exceptions import SubprocessError
from pypiper.manager import COMPLETE_FLAG, FAIL_FLAG
from sqlalchemy.engine import Engine

from refgenie.db.tables import Asset
from refgenie.logger import logger
from refgenie.models import AssetRegistryPathComponents, BuildCommandValues, BuildParams
from refgenie.managers.base import ResourceManager
from refgenie.utils.build import (
    BUILD_DIGEST_SCHEME,
    BuildProvenance,
    build_digest,
    build_level2,
    checksum,
    get_build_dir,
    handle_build_sigint,
)
from refgenie.managers.asset.colocation import (
    create_colocation_symlinks,
    get_colocation_metadata,
)
from refgenie.utils.io import coerce_cli_kwargs
from refgenie.utils.templating import jinja_render_template_strictly

if TYPE_CHECKING:
    from refgenie.managers.recipe import RecipeManager
    from refgenie.managers.asset_class import AssetClassManager
    from refgenie.managers.asset.manager import AssetManager
    from refgenie.managers.alias import AliasManager
    from refgenie.managers.asset.relations import AssetRelations
    from refgenie.managers.genome import GenomeManager


class AssetBuilder(ResourceManager):
    """
    Internal service for building assets using recipes and pypiper.

    This is NOT user-facing - it's used by Refgenie.build_asset() internally.
    """

    def __init__(
        self,
        database_engine: Engine,
        genome_folder: Path,
        alias_folder: Path,
        data_folder: Path,
        recipe_manager: "RecipeManager",
        asset_class_manager: "AssetClassManager",
        asset_manager: "AssetManager",
        alias_manager: "AliasManager",
        asset_relations: "AssetRelations",
        genome_manager: "GenomeManager",
    ):
        """
        Initialize the AssetBuilder.

        Args:
            database_engine: The database engine.
            genome_folder: Path to genome data folder.
            alias_folder: Path to alias folder.
            data_folder: Path to data folder.
            recipe_manager: The RecipeManager for getting recipes.
            asset_class_manager: The AssetClassManager.
            asset_manager: The AssetManager for asset operations.
            alias_manager: The AliasManager for alias operations.
            asset_relations: The AssetRelations for parent/child relationships.
            genome_manager: The GenomeManager for genome operations.
        """
        super().__init__(database_engine)
        self._genome_folder = genome_folder
        self._alias_folder = alias_folder
        self._data_folder = data_folder
        self._recipe_manager = recipe_manager
        self._asset_class_manager = asset_class_manager
        self._asset_manager = asset_manager
        self._alias_manager = alias_manager
        self._asset_relations = asset_relations
        self._genome_manager = genome_manager

    @staticmethod
    def _populate_commands(
        command_templates: list[str], command_values: BuildCommandValues
    ) -> list[str]:
        """
        Populate the command templates with values.

        Args:
            command_templates: The command templates.
            command_values: The values to populate the templates with.

        Returns:
            list[str]: The populated command templates.
        """
        return [jinja_render_template_strictly(cmd, command_values) for cmd in command_templates]

    @staticmethod
    def _validate_build_inputs(recipe, command_values: BuildCommandValues) -> None:
        """Validate that all required build inputs are available before template rendering.

        Raises MissingBuildInputError with a clear message if any required input is missing.
        """
        errors = []

        # Check required input files
        if recipe.input_files:
            for file_name, file_spec in recipe.input_files.items():
                if command_values.files is None or file_name not in command_values.files:
                    if file_spec.get("default") is None:
                        errors.append(
                            f"Missing required file '{file_name}': "
                            f"{file_spec.get('description', '')}. "
                            f"Provide it with: --files {file_name}=/path/to/file"
                        )

        # Check required input params (without defaults)
        if recipe.input_params:
            for param_name, param_spec in recipe.input_params.items():
                if (
                    command_values.params is None or param_name not in command_values.params
                ) and param_spec.get("default") is None:
                    errors.append(
                        f"Missing required parameter '{param_name}': "
                        f"{param_spec.get('description', '')}. "
                        f"Provide it with: --params {param_name}=value"
                    )

        # Check that refget_store_path exists when set (common for fasta recipe)
        if command_values.refget_store_path:
            store_path = Path(command_values.refget_store_path)
            if not store_path.exists():
                errors.append(
                    f"RefgetStore not found at '{store_path}'. "
                    f"Initialize the genome first with: refgenie genome init --fasta /path/to/file.fa"
                )

        if errors:
            error_list = "\n  - ".join(errors)
            from refgenie.exceptions import MissingBuildInputError

            raise MissingBuildInputError(
                f"Cannot build '{command_values.asset_group_name}': missing required inputs:\n"
                f"  - {error_list}\n"
                f"Use 'refgenie build --requirements {command_values.asset_group_name}' "
                f"to see all required inputs."
            )

    @staticmethod
    def _run_in_docker(cmd: str, docker_image: str) -> bytes:
        """
        Run a command in a docker container.

        Args:
            cmd: The command to run.
            docker_image: The docker image.

        Returns:
            bytes: The result of the command.
        """
        get_id_cmd = f"docker run -itd --rm {docker_image}"
        container_id = check_output(get_id_cmd, shell=True).decode("utf-8").strip()
        get_result_cmd = f"docker exec -it {container_id} {cmd}"
        return check_output(get_result_cmd, shell=True)

    def _resolve_custom_seek_keys(
        self, custom_seek_keys: dict[str, str], docker_image: str | None = None
    ) -> dict[str, Any]:
        """
        Resolve custom seek keys by executing shell commands.

        Args:
            custom_seek_keys: Dict mapping seek key names to shell commands.
            docker_image: If provided, run commands in this docker container.

        Returns:
            dict[str, Any]: A dictionary of resolved custom seek key values.
        """
        return (
            {
                key: (
                    self._run_in_docker(commands, docker_image)
                    if docker_image
                    else check_output(commands, shell=True)
                )
                .decode("utf-8")
                .strip()
                for key, commands in custom_seek_keys.items()
            }
            if custom_seek_keys
            else {}
        )

    @staticmethod
    def _input_file_digests(command_values: BuildCommandValues) -> dict[str, str]:
        """sha256 per input file that exists on disk, for the build digest."""
        return {
            name: checksum(str(path))
            for name, path in (command_values.files or {}).items()
            if Path(path).is_file()
        }

    @staticmethod
    def _build_provenance(
        command_values: BuildCommandValues,
        recipe: Any,
        genome_digest: str,
        input_asset_digests: dict[str, str],
        input_file_digests: dict[str, str],
        docker: bool = False,
        docker_image: str | None = None,
    ) -> BuildProvenance:
        """
        Record what this build was, for its ``AssetName`` row.

        Provenance belongs to the build, not to the bytes it produced, so it
        goes on the name row. It used to be injected as seek keys on the asset
        row, where a second build producing identical content silently lost it.
        """
        from importlib.metadata import version

        try:
            refgenie_version = version("refgenie")
        except Exception:
            refgenie_version = "unknown"

        resolved_image_digest = None
        if docker and docker_image:
            try:
                resolved_image_digest = (
                    check_output(
                        ["docker", "inspect", "--format", "{{index .RepoDigests 0}}", docker_image],
                    )
                    .decode("utf-8")
                    .strip()
                ) or None
            except (CalledProcessError, FileNotFoundError, OSError):
                pass

        build_timestamp = datetime.now(timezone.utc)
        level2 = build_level2(
            genome_digest=genome_digest,
            recipe_name=recipe.name,
            recipe_version=recipe.version,
            input_asset_digests=input_asset_digests,
            input_file_digests=input_file_digests,
            params=command_values.params,
            build_timestamp=build_timestamp,
            refgenie_version=refgenie_version,
            docker_image=docker_image if docker else None,
            docker_image_digest=resolved_image_digest,
        )
        digest, level1 = build_digest(level2)
        return BuildProvenance(
            build_digest=digest,
            build_level1=level1,
            build_digest_scheme=BUILD_DIGEST_SCHEME,
            build_timestamp=build_timestamp,
            refgenie_version=refgenie_version,
            inputs=level2,
            docker_image=docker_image if docker else None,
            docker_image_digest=resolved_image_digest,
            recipe_id=recipe.id,
        )

    @staticmethod
    def _resolve_default_asset(default_asset: str, namespaces: BuildCommandValues) -> str:
        """
        Resolve the default asset name using template rendering.

        Args:
            default_asset: The default asset template string.
            namespaces: The build command values.

        Returns:
            str: The resolved default asset name.

        Raises:
            ValueError: If the template resolves to an empty name. See
                ``Refgenie.resolve_default_asset`` for why falling back to
                "default" here silently publishes mis-named assets.
        """
        resolved = jinja_render_template_strictly(default_asset, namespaces)
        if not resolved or not resolved.strip():
            raise ValueError(
                f"default_asset template {default_asset!r} resolved to an empty asset "
                f"name. Refusing to fall back to 'default': that would build and "
                f"publish a mis-named asset. Check that the recipe's custom_seek_keys "
                f"commands run in this environment and actually print a version."
            )
        return resolved

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
            genome_folder=self._genome_folder,
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

        The single formula for the flag path. The advertised template, the
        write in ``build()``, and the skip-build guard all derive from here, so
        they cannot drift apart.

        Args:
            genome_name: The genome alias name (or a ``{genome_name}`` placeholder).
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.

        Returns:
            Path: The completion flag path.
        """
        build_dir = self._get_build_dir(genome_name, asset_group_name, asset_name)
        return build_dir / f"{genome_name}_{asset_group_name}__{asset_name}.flag"

    def _get_build_target_template(
        self,
        asset_group_name: str,
        asset_name: str,
    ) -> Path:
        """
        Get the build target template path.

        The ``{genome_name}`` placeholder is left unresolved for snakemake
        wildcard substitution. Substituting it yields exactly the path
        ``build()`` writes the completion flag to — the advertised path and the
        real path are the same file, with no symlink reconciliation between them.

        Args:
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.

        Returns:
            Path: The build target template path.
        """
        return self._get_build_flag("{genome_name}", asset_group_name, asset_name)

    def _create_symlinks(
        self,
        alias_name: str,
        asset_group_name: str,
        asset_name: str,
    ) -> None:
        """
        Create the name-addressed alias symlinks for the built asset.

        The content lives under a digest-named directory; the alias tree is the
        name-addressed view of it (``alias/<alias>/<group>/<name>/``). A build
        produces exactly one name, so this renders that one name across every
        genome alias.

        Args:
            alias_name: The alias name for the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset.
        """
        genome_digest = self._alias_manager.resolve(alias_name)
        asset_name = asset_name or self._asset_manager.get_default(
            asset_group_name, genome_digest=genome_digest
        )
        self._asset_manager._render_asset_alias_symlinks(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        )

    def _resolve_input_assets(
        self,
        genome_digest: str,
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
        # check if recipe input assets are available
        if (
            inputs := self._recipe_manager.get(recipe_name, recipe_version).input_asset_classes
        ) is None:
            return None
        # resolve user-provided input assets
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
                input_asset = self._asset_manager.get(
                    genome_digest=(
                        self._alias_manager.resolve(parsed.genome)
                        if parsed.genome
                        else genome_digest
                    ),
                    asset_group_name=parsed.asset_group,
                    asset_name=parsed.asset
                    or self._asset_manager.get_default(
                        parsed.asset_group,
                        genome_name=parsed.genome,
                    ),
                )
                asset_by_input_asset_name[asset_input_name] = input_asset
        # resolve remaining input assets. Ones that were not provided by the user (= None)
        for input_asset_name, default in asset_by_input_asset_name.items():
            if default is not None:
                # skip if the input asset class has been provided by the user
                continue
            for input in inputs:
                if input.name != input_asset_name:
                    continue
                input_asset = self._asset_manager.get(
                    genome_digest=genome_digest,
                    asset_group_name=input.default,
                    asset_name=self._asset_manager.get_default(
                        asset_group_name=input.default, genome_digest=genome_digest
                    )
                    or "default",
                )
                asset_by_input_asset_name[input_asset_name] = input_asset

        logger.debug(f"Resolved input assets: {asset_by_input_asset_name}")

        return asset_by_input_asset_name  # type: ignore

    def build(
        self,
        recipe_name: str,
        *,
        genome_name: str,
        genome_digest: str,
        asset_group_name: str,
        asset_name: str | None = None,
        recipe_version: str | None = None,
        params: BuildParams | None = None,
        docker: bool = False,
        docker_volumes: list[str] | None = None,
        asset_description: str | None = None,
        pipeline_kwargs: dict[str, Any | None] = None,
    ) -> Asset | None:
        """
        Build an asset using a recipe.

        Args:
            recipe_name: The name of the recipe to use.
            genome_name: The name of the genome.
            genome_digest: The digest of the genome.
            asset_group_name: The name of the asset group.
            asset_name: The name of the asset to build.
            recipe_version: The version of the recipe.
            params: The build parameters.
            docker: Whether to use docker.
            docker_volumes: Docker volumes to mount.
            asset_description: Description of the asset.
            pipeline_kwargs: Additional kwargs for PipelineManager.

        Returns:
            Asset | None: The built asset (or the existing asset when the
                build is skipped), or None if the build failed.
        """
        recipe = self._recipe_manager.get(recipe_name, recipe_version)
        # Extract asset class name while recipe is still in session
        output_asset_class_name = recipe.output_asset_class.name

        build_target_string = f"{genome_name}/{asset_group_name}"
        logger.info(f"Building '{build_target_string}' using recipe '{recipe}'")

        custom_seek_keys = (
            self._resolve_custom_seek_keys(
                recipe.custom_seek_keys or {}, recipe.docker_image if docker else None
            )
            or {}
        )
        logger.debug(f"Build parameters: {params}")
        logger.debug(f"Custom seek keys: {custom_seek_keys}")
        if params is not None:
            params.populate_with_defaults_from_recipe(recipe)
        asset_by_input_asset_name = self._resolve_input_assets(
            genome_digest=genome_digest,
            build_params=params,
            recipe_name=recipe_name,
            recipe_version=recipe_version,
        )
        command_values = BuildCommandValues(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            params=params.params if params is not None else {},
            files=params.files if params is not None else {},
            assets=asset_by_input_asset_name,
            custom_seek_keys=custom_seek_keys,
            genome_folder=self._genome_folder,
            refget_store_path=str(self._genome_folder / ".refget_store"),
        )
        # Validate required inputs before template rendering
        self._validate_build_inputs(recipe, command_values)
        # resolve default asset if not provided using command build values
        asset_name = asset_name or self._resolve_default_asset(recipe.default_asset, command_values)

        build_target_string += f":{asset_name}"
        # check if the asset already exists, and if so, skip the build
        if self._asset_manager.exists(
            genome_digest=genome_digest,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
        ):
            logger.warning(f"Asset '{build_target_string}' already exists. Skipping build")
            # The build is skipped, but snakemake still checks the completion
            # flag declared by _get_build_target_template(). It normally exists
            # already, having been written when the asset was built. The builds/
            # tree has no sweeper though, so it can be pruned independently of
            # the asset it describes — and a skipped build would then leave the
            # declared output missing and raise MissingOutputException. Touching
            # the one canonical path is idempotent.
            flag = self._get_build_flag(genome_name, asset_group_name, asset_name)
            if not flag.exists():
                logger.info(f"Restoring missing build flag for skipped build: {flag}")
                flag.parent.mkdir(parents=True, exist_ok=True)
                flag.touch()
            # Return existing asset so caller can archive if needed
            return self._asset_manager.get(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
            )

        # compose build output folder and update command values with the result
        build_output_folder = self._data_folder / genome_digest / asset_group_name / asset_name
        command_values.output_folder = build_output_folder

        # Create colocation symlinks BEFORE recipe commands run.
        # Tools like BWA need the parent file (e.g., .fa) present in the output
        # directory when running commands like `bwa index output/genome.fa`.
        create_colocation_symlinks(
            output_folder=build_output_folder,
            genome_folder=self._genome_folder,
            input_assets=recipe.input_assets,
            resolved_assets=asset_by_input_asset_name,
        )

        commands = self._populate_commands(
            command_templates=recipe.command_templates,
            command_values=command_values,
        )
        # Build bookkeeping lives outside the asset directory, in the builds/
        # tree. This path is identical to _get_build_target_template() with
        # {genome_name} substituted, so snakemake's declared output is the
        # real flag rather than a symlink to one written elsewhere.
        build_stats_output_folder = self._get_build_dir(genome_name, asset_group_name, asset_name)
        build_flag = self._get_build_flag(genome_name, asset_group_name, asset_name)
        target = build_flag.as_posix()
        # The flag is written below, after add_from_path commits -- never as a
        # recipe command: a flag written before the catalog row means a build
        # killed in the gap leaves snakemake and pypiper skipping the rebuild
        # forever. Its meaning is "this asset is in the catalog".
        logger.debug(f"Commands to be executed: {commands}")

        # We only get here when the catalog has no such asset, so any flag left
        # on disk is debris from an interrupted build. Left in place it would
        # make pypiper skip every command and hand an empty output folder to
        # add_from_path.
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
        if docker:
            logger.info(f"Getting docker container for image: {recipe.docker_image}")
            docker_volumes = docker_volumes or []
            docker_volumes.append(build_output_folder.as_posix())
            # Also mount the genome folder so colocation symlinks resolve in-container.
            # Colocation symlinks in the output folder are RELATIVE and point into
            # sibling asset folders (e.g. ../../fasta/<name>/<digest>.fa). Mounting only
            # the output folder leaves those targets unmounted and dangling, causing
            # tools like `bwa index` to fail with "No such file or directory". The
            # genome folder contains both the child output and the parent asset, so
            # mounting it makes the relative symlink targets resolve inside the container.
            genome_folder_mount = self._genome_folder.as_posix()
            if genome_folder_mount not in docker_volumes:
                docker_volumes.append(genome_folder_mount)
            pm.get_container(image=recipe.docker_image, mounts=docker_volumes)

        return_code = None
        try:
            failed = False
            # run build command
            #
            # signal.signal() only works in the main thread of the main
            # interpreter; off the main thread it raises ValueError. A build
            # driven from a worker thread -- the local web UI's job manager
            # runs every build on a dedicated pool thread -- must not die
            # here. AssetPuller guards its own SIGINT registration the same
            # way. In the CLI this branch is always taken.
            if threading.current_thread() is threading.main_thread():
                signal.signal(
                    signal.SIGINT,
                    handle_build_sigint(genome_name, asset_group_name, asset_name),
                )
            if commands:
                return_code = pm.run(commands, target, container=pm.container)
            else:
                # A recipe with no commands produces its output some other way
                # (colocation symlinks alone, say). Nothing ran, nothing failed.
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

        if failed:
            return None

        # The same parent digests collected for set_parents below; resolved once.
        input_asset_digests = {
            name: asset.digest
            for name, asset in asset_by_input_asset_name.items()
            if asset is not None and asset.digest
        }
        build_provenance = self._build_provenance(
            command_values=command_values,
            recipe=recipe,
            genome_digest=genome_digest,
            input_asset_digests=input_asset_digests,
            input_file_digests=self._input_file_digests(command_values),
            docker=docker,
            docker_image=recipe.docker_image if docker else None,
        )

        # Add the built asset using AssetManager
        colocate_metadata = get_colocation_metadata(recipe.input_assets)
        added_asset = self._asset_manager.add_from_path(
            asset_class_name=output_asset_class_name,
            path=build_output_folder.relative_to(self._genome_folder),
            genome_name=genome_name,
            asset_group_name=asset_group_name,
            asset_name=asset_name,
            description=asset_description or recipe.description,
            recipe=recipe,
            custom_seek_keys=custom_seek_keys,
            colocate=colocate_metadata,
            set_default=True,  # Deliberate builds become the group default
            build_provenance=build_provenance,
        )
        logger.info(f"Added asset: '{added_asset}'")
        # Post-operation: write the completion flag, then create symlinks. The
        # flag follows the commit that made the asset real; anything that reads
        # it to decide whether to rebuild is therefore reading the catalog's
        # answer, not the recipe's.
        if added_asset is not None:
            build_flag.parent.mkdir(parents=True, exist_ok=True)
            build_flag.touch()
            self._create_symlinks(
                alias_name=genome_name,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
            )

        # Set parent relationships if there are input assets
        if input_asset_digests:
            self._asset_relations.set_parents(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
                parent_asset_digests=list(input_asset_digests.values()),
            )

        return added_asset
