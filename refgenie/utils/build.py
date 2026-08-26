"""
Build utilities for refgenie.

NOTE: Keep imports in utils/ modules lightweight. These modules are imported
eagerly at package load time. Heavy dependencies (pypiper, pandas, etc.) should
be confined to manager classes which support lazy loading.
"""

import hashlib
import sys
from dataclasses import dataclass
from datetime import datetime
from fnmatch import fnmatch
from functools import partial
from pathlib import Path
from typing import Any

from ubiquerg import checksum as _checksum

from refgenie.const import BUILDS_DIR, DEFAULT_PULL_SIZE_CUTOFF_GB
from refgenie.utils.console import CONSOLE
from refgenie.logger import logger
from refgenie.utils.prompt import Confirmer, resolve_confirmer

checksum = partial(_checksum, algorithm="sha256")

#: Version tag mixed into every directory digest preimage. Bump this whenever
#: the set of hashed inputs or the framing changes, so that a reinterpretation
#: of the inputs is explicit rather than silent.
DIGEST_SCHEME = b"refgenie-dir-digest-v2"

#: Chunk size for streaming file reads, so multi-GB indexes are not loaded
#: into memory to be hashed.
_DIGEST_CHUNK_SIZE = 1024 * 1024

#: The ``inherent`` list a recipe gets when it does not declare one: everything
#: the build produced is part of the asset's identity.
#:
#: Include-all is deliberate, because the two failure directions are not
#: symmetric. Forgetting to *include* a load-bearing file makes two genuinely
#: different assets share a digest -- a silent identity collision on the asset
#: primary key. Accidentally *including* an incidental file only produces a
#: spurious digest difference, which is noisy but visible. Over-sensitive is
#: recoverable; under-sensitive is not, so declarations should mostly subtract.
DEFAULT_INHERENT = ["*"]

#: Appended after every recipe-declared ``inherent`` list, so it always has the
#: last word. macOS metadata sidecars are never part of an asset's identity and
#: a recipe must not be able to re-include them.
ALWAYS_EXCLUDED = ["!._*"]


def get_build_dir(
    genome_folder: Path,
    genome_name: str,
    asset_group_name: str | None = None,
    asset_name: str | None = None,
) -> Path:
    """
    Get the build bookkeeping directory for a build invocation.

    This is the single formula for the ``builds/`` layout; every call site
    derives its path from here rather than joining the components itself.
    The tree is keyed by build invocation (alias/genome name, group, asset),
    NOT by any asset addressing scheme — it must be computable before the
    build runs, and it must not follow asset content when content moves.

    Args:
        genome_folder: The refgenie genome folder.
        genome_name: The genome alias name (or a ``{genome_name}`` template
            placeholder for snakemake substitution).
        asset_group_name: The name of the asset group (optional).
        asset_name: The name of the asset (optional).

    Returns:
        Path: The build directory, resolved as deeply as the given arguments allow.
    """
    build_dir = genome_folder / BUILDS_DIR / genome_name
    if asset_group_name is None:
        return build_dir
    build_dir = build_dir / asset_group_name
    if asset_name is None:
        return build_dir
    return build_dir / asset_name


def _file_digest(path: Path) -> bytes:
    """
    SHA-256 of a single file's bytes, read in chunks.

    Args:
        path: The path to the file.

    Returns:
        bytes: The raw 32-byte digest (not hex).
    """
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(_DIGEST_CHUNK_SIZE):
            h.update(chunk)
    return h.digest()


def directory_size(path: Path) -> int:
    """
    Total size in bytes of every regular file under ``path``, recursively.

    Deliberately callable from outside a transaction: it walks the whole asset,
    which is far too much work to do inside an insert flush holding the write
    lock.

    Args:
        path: The directory to measure.

    Returns:
        int: The summed size of the files it contains.
    """
    return sum(f.stat().st_size for f in Path(path).glob("**/*") if f.is_file())


def resolve_inherent(inherent: list[str] | None) -> list[str]:
    """
    Build the effective glob list for a recipe's declared ``inherent`` set.

    Args:
        inherent: The recipe's declared list, or None if it declared none.

    Returns:
        list[str]: The declared list (or the include-all default) with the
        always-excluded patterns appended, so they match last.
    """
    return [*(inherent if inherent is not None else DEFAULT_INHERENT), *ALWAYS_EXCLUDED]


def _is_inherent(rel: str, patterns: list[str]) -> bool:
    """
    Whether a file belongs to the asset's identity, per an ordered glob list.

    Semantics are gitignore-style **last match wins**: each pattern is tested
    in order, and the last one that matches decides. A ``!`` prefix negates,
    marking the match incidental rather than inherent. A path no pattern
    matches is not inherent, so a pure include-list ("only these files") works
    without a leading exclude-everything entry.

    The glob dialect is :func:`fnmatch.fnmatch` against the full normalized
    relative path (POSIX separators, no leading ``./``) rather than the
    basename, so directory-scoped patterns have somewhere to land. Note that
    ``*`` and ``?`` cross ``/`` in this dialect, so ``*.log`` matches
    ``sub/dir/x.log`` and ``**`` is merely a redundant spelling of ``*``.

    Args:
        rel: The normalized path, relative to the digested directory.
        patterns: The effective glob list from :func:`resolve_inherent`.

    Returns:
        bool: True if the file should contribute to the digest.
    """
    inherent = False
    for pattern in patterns:
        negated = pattern.startswith("!")
        if fnmatch(rel, pattern[1:] if negated else pattern):
            inherent = not negated
    return inherent


def get_dir_digest(path, inherent: list[str] | None = None) -> str:
    """
    Generate a SHA-256 digest over the inherent contents of a directory.

    Which files count is **recipe-declarable** rather than hard-coded. The
    recipe's compute commands are the only thing that structurally knows what
    files a build produces, so the recipe is what declares which of them are
    inherent to the asset's identity and which are incidental (logs, scratch,
    caches). See :func:`_is_inherent` for the glob dialect and match rules,
    and :data:`DEFAULT_INHERENT` for why omitting the declaration includes
    everything.

    The digest covers each file's **relative path and its content**, so a
    rename or a relocation into a subdirectory changes the digest. This
    matters because refgenie's payload is largely aligner indexes, whose
    filenames are load-bearing: bwa locates index parts by suffix off a common
    prefix, bowtie2 expects ``<base>.{1,2,3,4}.bt2``, and seek keys resolve to
    specific filenames. An index with correct bytes under wrong names is
    broken, and must not be digest-identical to a working one.

    Entries are sorted byte-wise on the encoded relative path, so the result
    does not depend on the builder's locale. Each entry is length-prefixed in
    the preimage, so ``("ab", "c")`` cannot collide with ``("a", "bc")``.

    Deliberate design decisions, all of them choices rather than consequences:

    - **Symlinks are excluded.** Colocation symlinks are external references
      to parent assets, created before the build runs (``managers/asset/colocation.py``)
      and deliberately stripped from archives. They are not this asset's
      content. Consequence: an asset with a colocation symlink and one without
      hash identically. That is intended.
    - **Empty directories are invisible.** Only files are hashed, not
      directory structure, so every content-free asset shares one digest.
      Accepted: there should be no content-free assets, and if several exist,
      sharing one empty directory is what content addressing does anyway.
    - **Mode and the executable bit are excluded.** Genomic index files are
      not executables, and including mode would make digests depend on the
      builder's umask. Git includes mode; refgenie deliberately does not.
    - **The scheme is versioned** via ``DIGEST_SCHEME`` in the preimage.

    Args:
        path: The path to the directory.
        inherent: The recipe's declared ordered glob list. None means the
            recipe declared none, which includes everything.

    Returns:
        str: The hex-encoded SHA-256 digest.

    Raises:
        Exception: If the directory cannot be read. A digest that cannot be
            computed must not be silently substituted -- it is written
            straight into the asset primary key.
    """
    root = Path(path)
    patterns = resolve_inherent(inherent)

    # rglob() on a missing directory yields nothing rather than raising, which
    # would quietly produce the empty-directory digest for a path that is not
    # there. Fail loudly instead.
    if not root.is_dir():
        raise NotADirectoryError(f"Cannot calculate digest; not a directory: {path}")

    with CONSOLE.status(f"[bold]Calculating digest for {path}..."):
        entries = []
        for p in root.rglob("*"):
            # is_symlink() must be checked first: is_file() follows symlinks.
            if p.is_symlink() or not p.is_file():
                continue
            rel = p.relative_to(root).as_posix()
            if not _is_inherent(rel, patterns):
                continue
            entries.append((rel.encode("utf-8"), _file_digest(p)))

        h = hashlib.sha256()
        h.update(DIGEST_SCHEME)
        for rel, file_digest in sorted(entries):
            h.update(len(rel).to_bytes(4, "big"))
            h.update(rel)
            h.update(file_digest)
    return h.hexdigest()


#: The attributes that identify *a build*, as opposed to its output. Mirrors
#: ``refget.const.DEFAULT_INHERENT_ATTRS`` one level down: seqcol declares
#: which attributes of a sequence collection are inherent to its identity;
#: this declares which attributes of a build are inherent to its identity.
#: Everything else -- ``input_asset_names``, ``build_timestamp``,
#: ``refgenie_version``, ``docker_image``, ``docker_image_digest`` -- is
#: digested into level 1 for re-derivability, but excluded here because it
#: varies between runs of the same build or is redundant with an inherent
#: attribute.
BUILD_INHERENT_ATTRS = [
    "genome",
    "recipe",
    "recipe_version",
    "input_assets",
    "input_files",
    "params",
]

#: Version tag for the build-digest algorithm. Persisted alongside every
#: digest so a later change to ``BUILD_INHERENT_ATTRS`` (or the level-2 shape)
#: is visible on old rows as a scheme mismatch rather than a silent
#: incomparability.
BUILD_DIGEST_SCHEME = "refgenie-build-1"


def _canonicalize_value(value: Any) -> Any:
    """
    Coerce a level-2 attribute to the form its digest is computed from.

    Params are interpolated as text into the recipe's shell command template,
    so ``30`` and ``"30"`` render the identical command and are the same
    build -- coercing every leaf to a string makes that hold at the digest
    layer. It also keeps every leaf inside ``canonical_str``'s documented
    domain, which raises on floats: a float param digests here instead of
    failing the build. ``None`` is dropped rather than stringified to
    ``"None"``, since a value that was never provided and one that was
    explicitly given as the string ``"None"`` must not collide.

    Args:
        value: One level-2 attribute's raw value (as returned by
            :func:`build_level2`), or a value nested inside one.

    Returns:
        Any: The same shape with every scalar leaf coerced to ``str`` (or
        omitted, if it was ``None``).
    """
    if isinstance(value, dict):
        return {k: _canonicalize_value(v) for k, v in value.items() if v is not None}
    if isinstance(value, (list, tuple)):
        return [_canonicalize_value(v) for v in value if v is not None]
    if value is None:
        return None
    return str(value)


def build_level2(
    genome_digest: str,
    recipe_name: str,
    recipe_version: str | None,
    input_asset_digests: dict[str, str] | None,
    input_file_digests: dict[str, str] | None,
    params: dict[str, Any] | None,
    build_timestamp: datetime | None = None,
    refgenie_version: str | None = None,
    docker_image: str | None = None,
    docker_image_digest: str | None = None,
) -> dict[str, Any]:
    """
    Assemble the level-2 record for one build: every attribute it produced,
    inherent and non-inherent alike, in its raw (uncoerced) form.

    ``input_assets`` is the sorted array of parent *content* digests, not the
    name -> digest map: it is what makes ``build_digest`` reconstructible from
    the catalog, since ``AssetLink`` stores parent digests with no name.
    ``input_asset_names`` keeps the map for humans, as a non-inherent
    attribute -- renaming which slot an input filled without changing what
    content fills it is not a different build.

    ``genome_digest`` and ``input_file_digests`` are load-bearing, not
    redundant with the input assets: ``fasta`` declares no params, files or
    input assets and reads from the RefgetStore by genome digest, so without
    the genome every genome's fasta build would collide.

    Which params are output-affecting is a per-recipe question, so callers
    pass the whole resolved dict. An incidental param (a thread count)
    therefore splits one build into two identities -- the safe direction,
    since a spurious identity is visible and a missing one is silent.

    Returns:
        dict[str, Any]: The level-2 record, with values in the form they were
        passed -- see :func:`build_level2_to_level1` for the coercion applied
        before digesting.
    """
    input_asset_digests = input_asset_digests or {}
    return {
        "genome": genome_digest,
        "recipe": recipe_name,
        "recipe_version": recipe_version,
        "input_assets": sorted(input_asset_digests.values()),
        "input_asset_names": dict(input_asset_digests),
        "input_files": dict(input_file_digests or {}),
        "params": dict(params or {}),
        "build_timestamp": build_timestamp.isoformat() if build_timestamp else None,
        "refgenie_version": refgenie_version,
        "docker_image": docker_image,
        "docker_image_digest": docker_image_digest,
    }


def build_level2_to_level1(level2: dict[str, Any]) -> dict[str, str]:
    """
    Digest every attribute of a level-2 record individually.

    Mirrors ``refget.utils.seqcol_dict_to_level1_dict``: level 1 holds a digest
    of *every* attribute, not just the inherent ones, so the inherent set can
    be changed later by re-filtering stored level 1 rather than needing the
    original level 2 values again.

    Args:
        level2: The level-2 record, as built by :func:`build_level2`.

    Returns:
        dict[str, str]: Per-attribute digests -- the level-1 object.
    """
    from refget import canonical_str
    from refget.digests import sha512t24u_digest

    return {
        name: sha512t24u_digest(canonical_str(_canonicalize_value(value)))
        for name, value in level2.items()
    }


def build_level1_to_digest(
    level1: dict[str, str], inherent_attrs: list[str] = BUILD_INHERENT_ATTRS
) -> str:
    """
    Digest a level-1 object down to the single build digest.

    Filters to ``inherent_attrs`` first, so a build digest can be recomputed
    from a *stored* level 1 under a different inherent set without ever
    touching level 2 again -- the property the whole restructure exists for.

    Args:
        level1: Per-attribute digests, as returned by
            :func:`build_level2_to_level1`.
        inherent_attrs: Which attributes identify the build. Defaults to
            :data:`BUILD_INHERENT_ATTRS`.

    Returns:
        str: The build digest.
    """
    from refget import canonical_str
    from refget.digests import sha512t24u_digest

    filtered = {k: v for k, v in level1.items() if k in inherent_attrs}
    return sha512t24u_digest(canonical_str(filtered))


def build_digest(level2: dict[str, Any]) -> tuple[str, dict[str, str]]:
    """
    A deterministic digest identifying *a build*, as opposed to its output.

    ``Asset.digest`` addresses the bytes a build produced; this addresses the
    build that produced them, so two builds yielding byte-identical content
    stay distinguishable. The level-1 object this computes from is returned
    too, and persisted, so the digest is re-derivable from the row alone.

    Args:
        level2: The level-2 record, as built by :func:`build_level2`.

    Returns:
        tuple[str, dict[str, str]]: The build digest, and the level-1 object
            of per-attribute digests it was computed from.
    """
    level1 = build_level2_to_level1(level2)
    return build_level1_to_digest(level1), level1


@dataclass
class BuildProvenance:
    """
    What one build recorded about itself, bound for its ``AssetName`` row.

    Threaded from :meth:`AssetBuilder.build` through ``add_from_path`` to
    whichever name-insert path runs. A caller that did not build passes nothing
    and the columns stay NULL -- correct, not a value awaiting backfill.
    """

    build_digest: str | None = None
    #: The level-1 object :func:`build_digest` (the module function) computed
    #: the digest from -- every attribute's individual digest, so
    #: ``build_digest`` is re-derivable from the row with no other input.
    build_level1: dict[str, str] | None = None
    #: Which digest algorithm ``build_digest``/``build_level1`` were computed
    #: under. See :data:`BUILD_DIGEST_SCHEME`.
    build_digest_scheme: str | None = None
    build_timestamp: datetime | None = None
    refgenie_version: str | None = None
    #: The level-2 record :func:`build_level2` assembled: every attribute of
    #: the build, inherent and non-inherent, in its raw form.
    inputs: dict[str, Any] | None = None
    docker_image: str | None = None
    docker_image_digest: str | None = None
    #: The recipe this build ran -- the asset is shared across builds, the recipe is not.
    recipe_id: int | None = None

    def as_columns(self) -> dict[str, Any]:
        """The provenance as ``AssetName`` column keyword arguments."""
        return {
            "build_digest": self.build_digest,
            "build_level1": self.build_level1,
            "build_digest_scheme": self.build_digest_scheme,
            "build_timestamp": self.build_timestamp,
            "refgenie_version": self.refgenie_version,
            "inputs": self.inputs,
            "docker_image": self.docker_image,
            "docker_image_digest": self.docker_image_digest,
            "recipe_id": self.recipe_id,
        }


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


def handle_sigint_pull(filepath: Path):
    """Report an interrupted download and exit.

    Removing the partial file is not this handler's job: the download writes to
    a ``.part`` file beside ``filepath`` and unlinks it on any exception,
    including the ``SystemExit`` raised here. ``filepath`` itself never exists
    until the download has completed and been renamed into place.
    """

    def handle(sig, frame):
        logger.warning(f"\nThe download was interrupted: {filepath}")
        sys.exit(0)

    return handle


def should_pull_large_archive(
    archive_size: float | int,
    asset_registry_path: str,
    size_cutoff: float | int | None = None,
    force: bool | None = None,
    confirm: Confirmer | None = None,
) -> bool:
    """
    Check whether an archive over the size cutoff should be downloaded.

    Args:
        archive_size: The size of the archive, in bytes.
        asset_registry_path: The asset registry path.
        size_cutoff: The size cutoff, in GB.
        force: True to download regardless of size, False to skip, None to ask.
        confirm: Confirmation callback. Defaults to a refusal unless the CLI
            has enabled interactive prompts; see `refgenie.utils.prompt`.

    Returns:
        bool: Whether to pull the large archive.
    """
    size_cutoff_gb = size_cutoff if size_cutoff is not None else DEFAULT_PULL_SIZE_CUTOFF_GB
    size_cutoff_bytes = size_cutoff_gb * 1000**3
    logger.debug(f"'{asset_registry_path}' archive size: {archive_size}")
    if not force and archive_size > size_cutoff_bytes:
        if force is False:
            logger.info(f"Skipping pull of '{asset_registry_path}'; size: {archive_size}")
            return False
        # force is None, ask
        if not resolve_confirmer(confirm)(
            f"This archive exceeds the size cutoff "
            f"({archive_size / 1000**3:.1f}GB > {size_cutoff_gb:.1f}GB). Do you want to proceed?"
        ):
            logger.info(f"Skipping pull of '{asset_registry_path}'. Size too large: {archive_size}")
            return False
    return True
