"""
Asset identity: the content digest scheme and build provenance.

``get_dir_digest`` computes an asset's content digest (its primary key);
``build_digest`` and ``BuildProvenance`` identify the build that produced it.
``checksum`` is the SHA-256 file checksum used for downloads and staging.
"""

import hashlib
from dataclasses import dataclass
from datetime import datetime
from fnmatch import fnmatch
from functools import partial
from pathlib import Path
from typing import Any

from ubiquerg import checksum as _checksum

from refgenie.utils.console import CONSOLE

checksum = partial(_checksum, algorithm="sha256")

#: Version tag mixed into every directory digest preimage. Bump it whenever the
#: hashed inputs or the framing change.
DIGEST_SCHEME = b"refgenie-dir-digest-v2"

_DIGEST_CHUNK_SIZE = 1024 * 1024

#: The ``inherent`` list a recipe gets when it declares none: every file the
#: build produced is part of the asset's identity. See design-notes.md.
DEFAULT_INHERENT = ["*"]

#: Appended after every ``inherent`` list so it matches last; macOS metadata
#: sidecars are never part of an asset's identity.
ALWAYS_EXCLUDED = ["!._*"]


def _file_digest(path: Path) -> bytes:
    """Return the raw 32-byte SHA-256 of a file, read in chunks."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(_DIGEST_CHUNK_SIZE):
            h.update(chunk)
    return h.digest()


def directory_size(path: Path) -> int:
    """
    Total size in bytes of every regular file under ``path``, recursively.

    Call it outside a transaction: walking a whole asset is too slow to do
    while holding the write lock.
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

    The last matching pattern wins; a ``!`` prefix negates. A path no pattern
    matches is not inherent. Patterns are :func:`fnmatch.fnmatch` globs on the
    full relative POSIX path, so ``*`` crosses ``/`` (``*.log`` matches
    ``sub/x.log``).

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
    SHA-256 digest over the inherent files of a directory.

    Hashes each inherent file's relative path and content, sorted byte-wise
    and length-prefixed, under ``DIGEST_SCHEME``. Symlinks, empty directories
    and file modes are excluded. See "Asset content digest" in design-notes.md.

    Args:
        path: The path to the directory.
        inherent: The recipe's ordered glob list; None includes everything.

    Returns:
        str: The hex-encoded SHA-256 digest.

    Raises:
        NotADirectoryError: If ``path`` is not a directory.
    """
    root = Path(path)
    patterns = resolve_inherent(inherent)

    # rglob() on a missing directory yields nothing; fail instead.
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


#: The attributes that identify a build, like seqcol's inherent attributes.
#: The others are kept in level 1 but vary between runs or are redundant.
BUILD_INHERENT_ATTRS = [
    "genome",
    "recipe",
    "recipe_version",
    "input_assets",
    "input_files",
    "params",
]

#: Version tag for the build-digest algorithm, stored with every digest so a
#: change to ``BUILD_INHERENT_ATTRS`` shows up as a scheme mismatch.
BUILD_DIGEST_SCHEME = "refgenie-build-1"


def _canonicalize_value(value: Any) -> Any:
    """
    Coerce a level-2 attribute to the form its digest is computed from.

    Every leaf becomes a string, since params render as text into the shell
    command (``30`` and ``"30"`` are the same build) and ``canonical_str``
    rejects floats. ``None`` is dropped so it cannot collide with ``"None"``.

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
    container_engine: str | None = None,
) -> dict[str, Any]:
    """
    Assemble the level-2 record for one build: every attribute, raw.

    ``input_assets`` is the sorted list of parent content digests, so the
    build digest can be rebuilt from ``AssetLink`` rows, which store no names;
    ``input_asset_names`` keeps the name map as a non-inherent attribute. The
    genome digest is needed because ``fasta`` has no other inputs.
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
        "container_engine": container_engine,
    }


def build_level2_to_level1(level2: dict[str, Any]) -> dict[str, str]:
    """
    Digest every attribute of a level-2 record, returning the level-1 object.

    Level 1 covers every attribute, not just the inherent ones, so the
    inherent set can change later by re-filtering stored level 1.
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
    """Digest the ``inherent_attrs`` of a level-1 object into the build digest."""
    from refget import canonical_str
    from refget.digests import sha512t24u_digest

    filtered = {k: v for k, v in level1.items() if k in inherent_attrs}
    return sha512t24u_digest(canonical_str(filtered))


def build_digest(level2: dict[str, Any]) -> tuple[str, dict[str, str]]:
    """
    Digest identifying a build, as opposed to its output (``Asset.digest``).

    Returns:
        tuple[str, dict[str, str]]: The build digest and the level-1 object it
            was computed from, which is stored so the digest can be re-derived.
    """
    level1 = build_level2_to_level1(level2)
    return build_level1_to_digest(level1), level1


@dataclass
class BuildProvenance:
    """
    What one build recorded about itself, bound for its ``AssetName`` row.

    Assets that were not built locally leave every field None.
    """

    build_digest: str | None = None
    #: The level-1 object ``build_digest`` was computed from.
    build_level1: dict[str, str] | None = None
    build_digest_scheme: str | None = None
    build_timestamp: datetime | None = None
    refgenie_version: str | None = None
    #: The raw level-2 record from :func:`build_level2`.
    inputs: dict[str, Any] | None = None
    docker_image: str | None = None
    docker_image_digest: str | None = None
    #: The recipe this build ran; an asset can be shared across builds.
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
