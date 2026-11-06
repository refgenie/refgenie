"""
Build steps that need no BuildManager state: input checks, custom seek key
probes, command rendering and build provenance.
"""

from datetime import datetime, timezone
from pathlib import Path
from subprocess import DEVNULL, CalledProcessError, TimeoutExpired, check_output
from typing import Any

from refgenie.exceptions import CustomSeekKeyError, MissingBuildInputError
from refgenie.managers import containers
from refgenie.models import BuildCommandValues, GenomeDigest
from refgenie.utils.digest import (
    BUILD_DIGEST_SCHEME,
    BuildProvenance,
    build_digest,
    build_level2,
    checksum,
)
from refgenie.utils.templating import jinja_render_template_strictly

#: Seconds one custom-seek-key probe may run before it is killed. Generous
#: because a container probe may pull the image first; bounded because the web
#: build form runs probes on a request thread.
CUSTOM_SEEK_KEY_TIMEOUT = 120


def populate_commands(
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


def validate_build_inputs(recipe, command_values: BuildCommandValues) -> None:
    """Validate that all required build inputs are available before template rendering.

    Raises MissingBuildInputError with a clear message if any required input is missing.
    """
    errors = []

    if recipe.input_files:
        for file_name, file_spec in recipe.input_files.items():
            if command_values.files is None or file_name not in command_values.files:
                if file_spec.get("default") is None:
                    errors.append(
                        f"Missing required file '{file_name}': "
                        f"{file_spec.get('description', '')}. "
                        f"Provide it with: --files {file_name}=/path/to/file"
                    )

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


def probe_custom_seek_keys(
    custom_seek_keys: dict[str, str],
    image: str | None = None,
    engine: "containers.ContainerEngine | None" = None,
) -> dict[str, Any]:
    """
    Resolve custom seek keys by executing shell commands.

    Args:
        custom_seek_keys: Dict mapping seek key names to shell commands.
        image: If provided, run commands inside this container image.
        engine: The engine that runs ``image``; None means this host has
            neither Docker nor Apptainer.

    Returns:
        dict[str, Any]: A dictionary of resolved custom seek key values.

    Raises:
        CustomSeekKeyError: A command failed, could not be started (no
            container engine, daemon down, image unpullable), or ran past
            `CUSTOM_SEEK_KEY_TIMEOUT`.
    """
    resolved: dict[str, Any] = {}
    for key, command in (custom_seek_keys or {}).items():
        if image and engine is None:
            raise CustomSeekKeyError(
                key,
                command,
                "This computer has neither Docker nor Apptainer, so refgenie cannot "
                "start that image. Install one of them.",
                image,
            )
        try:
            output = (
                containers.one_shot(
                    engine, image, command, mounts=[], timeout=CUSTOM_SEEK_KEY_TIMEOUT
                )
                if image
                else check_output(
                    command, shell=True, stdin=DEVNULL, timeout=CUSTOM_SEEK_KEY_TIMEOUT
                ).decode("utf-8")
            )
        except TimeoutExpired as exc:
            hint = (
                f" If the image is not on this computer yet, fetch it first with"
                f" `{containers.pull_command(engine, image)}`."
                if image
                else ""
            )
            raise CustomSeekKeyError(
                key,
                command,
                f"That command was still running after {CUSTOM_SEEK_KEY_TIMEOUT} seconds, "
                f"so refgenie stopped it.{hint}",
                image,
                engine,
            ) from exc
        except CalledProcessError as exc:
            raise CustomSeekKeyError(
                key,
                command,
                f"That command failed, exiting with status {exc.returncode}.",
                image,
                engine,
            ) from exc
        except OSError as exc:
            raise CustomSeekKeyError(
                key,
                command,
                f"That command could not be started at all: {exc}.",
                image,
                engine,
            ) from exc
        value = output.strip()
        if not value:
            # A `tool --version | grep` pipeline exits 0 even when grep matches nothing.
            raise CustomSeekKeyError(
                key,
                command,
                "That command ran but printed nothing back, so there is no value to "
                "name the asset with. The tool itself is probably fine: the usual "
                "cause is that the tool changed how it prints its version, and the "
                "recipe no longer recognises it.",
                image,
                engine,
            )
        resolved[key] = value
    return resolved


def input_file_digests(command_values: BuildCommandValues) -> dict[str, str]:
    """sha256 per input file that exists on disk, for the build digest."""
    return {
        name: checksum(str(path))
        for name, path in (command_values.files or {}).items()
        if Path(path).is_file()
    }


def build_provenance(
    command_values: BuildCommandValues,
    recipe: Any,
    genome_digest: GenomeDigest,
    input_asset_digests: dict[str, str],
    input_file_digests: dict[str, str],
    image: str | None = None,
    engine: containers.ContainerEngine | None = None,
) -> BuildProvenance:
    """
    Record what this build was, for its ``AssetName`` row.

    ``image`` and ``engine`` are the container the recipe ran in, or None
    for a build on the host.
    """
    from importlib.metadata import version

    try:
        refgenie_version = version("refgenie")
    except Exception:
        refgenie_version = "unknown"

    resolved_image_digest = containers.image_digest(engine, image) if image else None

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
        docker_image=image,
        docker_image_digest=resolved_image_digest,
        container_engine=engine,
    )
    digest, level1 = build_digest(level2)
    return BuildProvenance(
        build_digest=digest,
        build_level1=level1,
        build_digest_scheme=BUILD_DIGEST_SCHEME,
        build_timestamp=build_timestamp,
        refgenie_version=refgenie_version,
        inputs=level2,
        docker_image=image,
        docker_image_digest=resolved_image_digest,
        recipe_id=recipe.id,
    )


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
