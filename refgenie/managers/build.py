"""
BuildManager - building assets from recipes with pypiper (``rgc.build``).
"""

import signal
import sys
import threading
from collections.abc import Callable
from datetime import datetime, timezone
from pathlib import Path
from subprocess import DEVNULL, CalledProcessError, TimeoutExpired, check_output
from typing import TYPE_CHECKING, Any

from pypiper import PipelineManager
from pypiper.exceptions import SubprocessError
from pypiper.manager import COMPLETE_FLAG, FAIL_FLAG
from sqlalchemy.engine import Engine

from refgenie.db.tables import Asset
from refgenie.exceptions import (
    CustomSeekKeyError,
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
from refgenie.utils.digest import (
    BUILD_DIGEST_SCHEME,
    BuildProvenance,
    build_digest,
    build_level2,
    checksum,
)
from refgenie.utils.io import coerce_cli_kwargs
from refgenie.utils.paths import get_build_dir
from refgenie.utils.templating import jinja_render_template_strictly

if TYPE_CHECKING:
    from refgenie.managers.alias import AliasBackend
    from refgenie.managers.asset.manager import AssetManager
    from refgenie.managers.genome import GenomeManager
    from refgenie.managers.recipe import RecipeManager
    from refgenie.managers.stage import StageManager


#: Seconds one custom-seek-key probe may run before it is killed. The probes
#: themselves are `tool --version` one-liners, but a docker probe may have to
#: pull the image first, which is why the bound is generous rather than tight.
#: It exists at all because the web build form resolves an asset name on a
#: request thread: `check_output` with no timeout on a `docker run` pins that
#: thread for the life of the process.
CUSTOM_SEEK_KEY_TIMEOUT = 120


def handle_build_sigint(genome_name: str, asset_group_name: str, asset_name: str):
    """
    Build a SIGINT handler that reports what was interrupted, then exits.

    Args:
        genome_name: The name of the genome being built.
        asset_group_name: The name of the asset group being built.
        asset_name: The name of the asset being built.

    Returns:
        The SIGINT handling function.
    """

    def handle(sig, frame):
        logger.warning(
            f"\nThe build was interrupted for {genome_name}/{asset_group_name}:{asset_name}"
        )
        sys.exit(0)

    return handle


class BuildManager(ResourceManager):
    """
    The public build manager, at ``rgc.build``: building assets from recipes.

    Public surface:

    - ``run``: build one asset (optionally pulling the parent genome and
      staging the result).
    - ``preflight``: validate a prospective build without starting it.
    - ``initialize_and_build``: initialize a genome from FASTA and build its
      ``fasta`` asset.
    - ``target_template`` / ``genome_init_target_template``: the Snakemake
      target paths for a built asset and for genome init.
    - ``resolve_custom_seek_keys``: a recipe's custom seek keys, cached, for
      naming an asset outside a build.
    - ``default_asset_name``: a recipe's default asset name, outside a build.
    - ``resolve_default_asset``: render a ``default_asset`` template.

    Everything else is a private step. It depends on ``rgc.asset`` one way;
    nothing in the asset package holds a reference to it.
    """

    def __init__(
        self,
        database_engine: Engine,
        recipe_manager: "RecipeManager",
        asset: "AssetManager",
        alias_manager: "AliasBackend",
        genome_manager: "GenomeManager",
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
            stage_manager: The StageManager, for ``run(stage=True)``.
            stage_folder_getter: Returns the genome stage folder (or None) at
                call time.
            pull_parent: Pulls a missing parent genome, for
                ``run(pull_parents=True)``. The only link from build to pull.
            events: Where ``pre_build`` / ``post_build`` are emitted for plugins.
        """
        super().__init__(database_engine)
        self._events = events or NULL_EVENTS
        self._recipe_manager = recipe_manager
        self._asset = asset
        self._alias_manager = alias_manager
        self._genome_manager = genome_manager
        self._stage_manager = stage_manager
        self._stage_folder_getter = stage_folder_getter
        self._pull_parent = pull_parent
        #: (recipe name, recipe version, docker image) -> resolved custom seek
        #: keys, for `resolve_custom_seek_keys`. See its docstring for why the
        #: naming path caches and the build path does not.
        self._custom_seek_key_cache: dict[tuple[str, str | None, str | None], dict[str, Any]] = {}

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
            raise MissingBuildInputError(
                f"Cannot build '{command_values.asset_group_name}': missing required inputs:\n"
                f"  - {error_list}\n"
                f"Use 'refgenie build --requirements {command_values.asset_group_name}' "
                f"to see all required inputs."
            )

    @staticmethod
    def _run_in_docker(
        cmd: str, docker_image: str, timeout: float = CUSTOM_SEEK_KEY_TIMEOUT
    ) -> bytes:
        """
        Run a one-shot command in a docker container and return its stdout.

        A single `docker run --rm`, never a `run -itd` + `exec -it` pair. That
        pair breaks two ways: `--rm` fires only when the container's own
        entrypoint exits, so every detached probe leaks a container that runs
        until the daemon stops; and `docker exec -it` refuses outright
        ("cannot attach stdin to a TTY-enabled container") whenever stdin is
        not a terminal -- which is every server thread, snakemake run and CI
        job. The refusal is silent: it goes to stderr, and the empty stdout
        becomes an empty tool version.

        Args:
            cmd: The command to run, interpreted by `sh` inside the container.
            docker_image: The docker image.
            timeout: Seconds before the container is killed.

        Returns:
            bytes: The stdout of the command.
        """
        return check_output(
            ["docker", "run", "--rm", "--entrypoint", "sh", docker_image, "-c", cmd],
            stdin=DEVNULL,
            timeout=timeout,
        )

    def resolve_custom_seek_keys(self, recipe: Any) -> dict[str, Any]:
        """Resolve a recipe's custom seek keys, to NAME an asset rather than build one.

        Used by the build form's preflight and by the snakefile template, both
        of which need a recipe's default asset name outside a build. The name
        is the version of the tool that builds the asset, and the only way to
        learn a tool's version is to run it.

        Where it runs matters. A recipe that declares a `docker_image` is
        saying its tools live in that image, not on this host, and probing the
        host for `bowtie2-build` on a machine that (correctly) does not have it
        was how the build form came to refuse every version-named recipe. So
        the probe runs in the declared image when there is one -- through
        `_probe_custom_seek_keys`, the same helper a `--docker` build already
        uses, not a second docker invocation. A recipe with no image still
        resolves on the host, exactly as before.

        Cached for the life of the process, keyed by recipe name, version and
        image: the build form preflights on a 400 ms debounce, and one
        container per keystroke is not a price worth paying for an answer that
        cannot change. Failures are not cached, so a retry after `docker pull`
        works. `run` deliberately does NOT come through here: what a build
        records has to be probed at build time.

        Args:
            recipe: The recipe whose `custom_seek_keys` to resolve.

        Returns:
            dict[str, Any]: seek key name -> resolved value. Empty, with
            nothing run and no container started, for a recipe that declares no
            custom seek keys -- which is every recipe naming its asset the
            literal `default`.

        Raises:
            CustomSeekKeyError: A probe failed, could not start, or timed out.
                `preflight` turns this into a field-scoped problem; the
                snakefile template turns it into a fatal parse-time error.
        """
        commands = recipe.custom_seek_keys or {}
        if not commands:
            return {}
        cache_key = (recipe.name, recipe.version, recipe.docker_image)
        if cache_key not in self._custom_seek_key_cache:
            self._custom_seek_key_cache[cache_key] = self._probe_custom_seek_keys(
                commands, recipe.docker_image
            )
        return dict(self._custom_seek_key_cache[cache_key])

    def default_asset_name(
        self, recipe: Any, command_values: BuildCommandValues | None = None
    ) -> str:
        """A recipe's default asset name, outside a build.

        The one naming path shared by ``preflight`` and the generated
        Snakefile: resolve the recipe's custom seek keys (cached), put them in
        the template namespace, and render ``default_asset``.

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
        return self.resolve_default_asset(recipe.default_asset, command_values)

    def _probe_custom_seek_keys(
        self, custom_seek_keys: dict[str, str], docker_image: str | None = None
    ) -> dict[str, Any]:
        """
        Resolve custom seek keys by executing shell commands.

        The single place in the package that runs one of these commands. The
        naming path (`resolve_custom_seek_keys` on this class, which is what
        the build form's preflight and the snakefile template reach) routes
        through here too, so a recipe's version is probed in exactly one way.

        Args:
            custom_seek_keys: Dict mapping seek key names to shell commands.
            docker_image: If provided, run commands in this docker container.

        Returns:
            dict[str, Any]: A dictionary of resolved custom seek key values.

        Raises:
            CustomSeekKeyError: A command failed, could not be started (no
                docker binary, daemon down, image unpullable), or ran past
                `CUSTOM_SEEK_KEY_TIMEOUT`.
        """
        resolved: dict[str, Any] = {}
        for key, command in (custom_seek_keys or {}).items():
            try:
                output = (
                    self._run_in_docker(command, docker_image)
                    if docker_image
                    else check_output(
                        command, shell=True, stdin=DEVNULL, timeout=CUSTOM_SEEK_KEY_TIMEOUT
                    )
                )
            except TimeoutExpired as exc:
                hint = (
                    f" If the image is not on this computer yet, fetch it first with"
                    f" `docker pull {docker_image}`."
                    if docker_image
                    else ""
                )
                raise CustomSeekKeyError(
                    key,
                    command,
                    f"That command was still running after {CUSTOM_SEEK_KEY_TIMEOUT} seconds, "
                    f"so refgenie stopped it.{hint}",
                    docker_image,
                ) from exc
            except CalledProcessError as exc:
                raise CustomSeekKeyError(
                    key,
                    command,
                    f"That command failed, exiting with status {exc.returncode}.",
                    docker_image,
                ) from exc
            except OSError as exc:
                # No `docker` binary, or no shell -- the command never started.
                raise CustomSeekKeyError(
                    key, command, f"That command could not be started at all: {exc}.", docker_image
                ) from exc
            value = output.decode("utf-8").strip()
            if not value:
                # Exit 0, nothing on stdout. These probes are all
                # `tool --version | grep ...` pipelines and a grep that matches
                # nothing still exits 0 through the pipe, so the usual cause is
                # a tool whose version string has drifted away from the pattern
                # the recipe was written against -- not a missing tool. This is
                # the last point that still knows the command and the image it
                # ran in, which is the difference between a fixable report and
                # "check that the commands run in this environment".
                raise CustomSeekKeyError(
                    key,
                    command,
                    "That command ran but printed nothing back, so there is no value to "
                    "name the asset with. The tool itself is probably fine: the usual "
                    "cause is that the tool changed how it prints its version, and the "
                    "recipe no longer recognises it.",
                    docker_image,
                )
            resolved[key] = value
        return resolved

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
        genome_digest: GenomeDigest,
        input_asset_digests: dict[str, str],
        input_file_digests: dict[str, str],
        docker: bool = False,
        docker_image: str | None = None,
    ) -> BuildProvenance:
        """
        Record what this build was, for its ``AssetName`` row.

        Provenance belongs to the build, not to the bytes it produced, so it
        goes on the name row. Stored on the asset row (for example as seek
        keys), a second build producing identical content would silently lose
        it.
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
    def resolve_default_asset(default_asset: str, namespaces: BuildCommandValues) -> str:
        """
        Resolve the default asset name using template rendering.

        Args:
            default_asset: The default asset template string.
            namespaces: The build command values.

        Returns:
            str: The resolved default asset name.

        Raises:
            ValueError: If the template resolves to an empty name.
        """
        # Deliberately NO `or "default"` fallback. Assets are named after the
        # version of the tool that built them, and that version is produced by a
        # shell pipeline in the recipe's custom_seek_keys. Those pipelines end in
        # grep/awk/head, so a tool that fails or stops reporting its version
        # yields an empty string with a zero exit status -- silently.
        #
        # Never fall back to "default": it would bake a mis-named asset into
        # build paths and published keys.
        #
        # Recipes that legitimately have no tool version (fasta, fasta_index)
        # declare the literal `default_asset: "default"`, which renders truthy
        # and is unaffected.
        resolved = jinja_render_template_strictly(default_asset, namespaces)
        if not resolved or not resolved.strip():
            raise ValueError(
                f"Refgenie could not work out a name for this asset. The recipe builds one "
                f"from the pattern {default_asset!r}, and that came out empty. You can type "
                f"your own name instead (the Asset name box, or --asset on the command "
                f"line) and the build will go ahead. Refgenie will not quietly call it "
                f"'default', because that name would say nothing about what is in it and "
                f"would then be baked into every path and published key."
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

        The single formula for the flag path. The advertised template, the
        write in ``_build``, and the skip-build guard all derive from here, so
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
        wildcard substitution. Substituting it yields exactly the path
        ``_build`` writes the completion flag to — the advertised path and the
        real path are the same file, with no symlink reconciliation between them.

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
        # resolve remaining input assets. Ones that were not provided by the user (= None)
        for input_asset_name, default in asset_by_input_asset_name.items():
            if default is not None:
                # skip if the input asset class has been provided by the user
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
        docker: bool = False,
        docker_volumes: list[str] | None = None,
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
            docker: Whether to use docker to build the asset.
            docker_volumes: The docker volumes to mount.
            asset_description: The description of the asset, used only if the asset does not exist.
            pull_parents: Whether to pull the parent genome if not found.
            pipeline_kwargs: Additional kwargs for PipelineManager.
            push_to: Optional remote names or ids to create push intent records for.
        Returns:
            Asset: The asset that was built, or None if skipped/failed
        """
        # Fail before the build, not after it. A misconfigured stage folder used
        # to surface only once the build had finished, which from a browser (or
        # a multi-hour bowtie2 index) means losing the whole run to a check that
        # costs nothing up front.
        if stage and self._stage_folder_getter() is None:
            raise ValueError(
                "Can't stage the asset: genome_stage_folder is not set. "
                "Set it in the refgenie config or build without --stage."
            )

        # Resolve genome digest, with optional pull_parents fallback
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
                docker=docker,
                docker_volumes=docker_volumes,
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
            # The stage folder was validated at the top of this method.
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

        The seam the web build form's preflight endpoint calls: it runs the
        same resolution and validation `run` would (genome, recipe, input
        assets, required files/params, staging config, default asset name) but
        starts no pipeline and writes nothing.

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
            self._validate_build_inputs(recipe, command_values)
        except MissingBuildInputError as exc:
            problem("params", "missing_build_input", str(exc))

        if asset_name:
            resolved["asset_name"] = asset_name
        else:
            # The default-asset template needs the recipe's custom seek keys
            # (tool versions, produced by shell one-liners run wherever the
            # recipe says its tools live). Failing here is not a failure of
            # the build: refgenie tried to fill the name in and could not, and
            # the same request with ``asset_name`` set preflights ok. So it is
            # reported under its own code, which a form reads as "make the
            # field required" rather than as a red error.
            try:
                resolved["asset_name"] = self.default_asset_name(recipe, command_values)
            except Exception as exc:  # noqa: BLE001 - report, don't crash a preflight
                # Passed through unprefixed. Both failures underneath already
                # open with "Refgenie could not work out a name for this
                # asset"; wrapping that in a second sentence saying the same
                # thing in refgenie's own vocabulary is what turned this into a
                # wall of text nobody could act on.
                problem("asset", "asset_name_required", str(exc))

        return {"ok": not errors, "errors": errors, "resolved": resolved}

    @update_scope
    def initialize_and_build(
        self,
        fasta_file_path: Path,
        genome_names: list[GenomeAlias],
        description: str = "",
        species_name: str | None = None,
        use_existing: bool = False,
        build_fasta: bool = True,
    ) -> tuple[str, bool]:
        """Initialize a genome from FASTA and optionally build the fasta asset.

        This is the recommended way to add a new genome. It ingests the FASTA
        into the RefgetStore and builds the fasta asset (fa, fai, chrom.sizes)
        in one step.

        Args:
            fasta_file_path: Path to the FASTA file.
            genome_names: Alias names for the genome.
            description: Genome description.
            species_name: Species name.
            use_existing: Allow re-init of existing genome.
            build_fasta: Whether to build the fasta asset (default True).

        Returns:
            Tuple of (genome_digest, was_created).
        """
        digest, created = self._genome_manager.initialize_genome(
            fasta_file_path=fasta_file_path,
            description=description,
            alias_names=genome_names,
            use_existing=use_existing,
            species_name=species_name,
        )

        if build_fasta:
            self.run(
                recipe_name="fasta",
                genome_alias=genome_names[0],
                asset_group_name="fasta",
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
        The values a recipe's command templates are rendered with.

        The one place a ``BuildCommandValues`` is made, for a build and for
        ``preflight`` alike.
        """
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
        and validates the required inputs. Each step raises its own exception
        unchanged. ``preflight`` calls the same steps one by one
        instead, because it reports every problem rather than the first.

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
        self._validate_build_inputs(recipe, command_values)
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
        docker: bool = False,
        docker_volumes: list[str] | None = None,
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
            docker: Whether to use docker.
            docker_volumes: Docker volumes to mount.
            asset_description: Description of the asset.
            pipeline_kwargs: Additional kwargs for PipelineManager.

        Returns:
            Asset | None: The built asset (or the existing asset when the
                build is skipped), or None if the build failed.
        """
        recipe = self._recipe_manager.get(recipe_name, recipe_version)
        logger.info(f"Building '{genome_name}/{asset_group_name}' using recipe '{recipe}'")

        custom_seek_keys = (
            self._probe_custom_seek_keys(
                recipe.custom_seek_keys or {}, recipe.docker_image if docker else None
            )
            or {}
        )
        logger.debug(f"Build parameters: {params}")
        logger.debug(f"Custom seek keys: {custom_seek_keys}")
        input_assets, command_values = self._prepare(
            recipe, genome_digest, asset_group_name, params, custom_seek_keys
        )
        # resolve default asset if not provided using command build values.
        # Not `default_asset_name`: that path uses the cached probe, and what a
        # build records has to be probed at build time (done just above).
        asset_name = asset_name or self.resolve_default_asset(recipe.default_asset, command_values)

        # check if the asset already exists, and if so, skip the build
        if self._asset.exists(asset_group_name, asset_name, genome_digest=genome_digest):
            logger.warning(
                f"Asset '{genome_name}/{asset_group_name}:{asset_name}' already exists. "
                f"Skipping build"
            )
            # The build is skipped, but snakemake still checks the completion
            # flag declared by target_template(). It normally exists
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
            return self._asset.get(asset_group_name, asset_name, genome_digest=genome_digest)

        # compose build output folder and update command values with the result
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
            docker=docker,
            docker_volumes=docker_volumes,
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
            docker=docker,
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
        docker: bool,
        docker_volumes: list[str] | None,
        pipeline_kwargs: dict[str, Any | None] | None,
    ) -> bool:
        """
        Run the recipe's commands through pypiper, into ``command_values.output_folder``.

        Returns:
            bool: Whether the pipeline succeeded.
        """
        build_target_string = f"{genome_name}/{asset_group_name}:{asset_name}"
        build_output_folder = command_values.output_folder

        # Create colocation symlinks BEFORE recipe commands run.
        # Tools like BWA need the parent file (e.g., .fa) present in the output
        # directory when running commands like `bwa index output/genome.fa`.
        create_colocation_symlinks(
            output_folder=build_output_folder,
            genome_folder=self._asset.genome_folder,
            input_assets=recipe.input_assets,
            resolved_assets=input_assets,
        )

        commands = self._populate_commands(
            command_templates=recipe.command_templates,
            command_values=command_values,
        )
        # Build bookkeeping lives outside the asset directory, in the builds/
        # tree. This path is identical to target_template() with
        # {genome_name} substituted, so snakemake's declared output is the
        # real flag rather than a symlink to one written elsewhere.
        build_stats_output_folder = self._get_build_dir(genome_name, asset_group_name, asset_name)
        build_flag = self._get_build_flag(genome_name, asset_group_name, asset_name)
        target = build_flag.as_posix()
        # The flag is written by _record, after content.add commits -- never as
        # a recipe command: a flag written before the catalog row means a build
        # killed in the gap leaves snakemake and pypiper skipping the rebuild
        # forever. Its meaning is "this asset is in the catalog".
        logger.debug(f"Commands to be executed: {commands}")

        # We only get here when the catalog has no such asset, so any flag left
        # on disk is debris from an interrupted build. Left in place it would
        # make pypiper skip every command and hand an empty output folder to
        # content.add.
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
            genome_folder_mount = self._asset.genome_folder.as_posix()
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
        docker: bool,
    ) -> Asset | None:
        """
        Record a finished build: the asset, its completion flag, parents and
        provenance. The alias tree is rendered inside ``content.add``.

        Returns:
            Asset | None: The added asset.
        """
        # The same parent digests collected for set_parents below; resolved once.
        input_asset_digests = {
            name: asset.digest
            for name, asset in (input_assets or {}).items()
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
        # Post-operation: write the completion flag. It follows the commit that
        # made the asset real (which also rendered the alias tree, inside
        # content.add); anything that reads the flag to decide whether to
        # rebuild is therefore reading the catalog's answer, not the recipe's.
        if added_asset is not None:
            build_flag = self._get_build_flag(genome_name, asset_group_name, asset_name)
            build_flag.parent.mkdir(parents=True, exist_ok=True)
            build_flag.touch()

        # Set parent relationships if there are input assets
        if input_asset_digests:
            self._asset.links.set_parents(
                genome_digest=genome_digest,
                asset_group_name=asset_group_name,
                asset_name=asset_name,
                parent_asset_digests=list(input_asset_digests.values()),
            )

        return added_asset
