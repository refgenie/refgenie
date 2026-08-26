"""Version resolution utilities for asset classes and recipes."""

import re

from sqlalchemy.exc import NoResultFound
from sqlmodel import Session, select


# Semver pattern: major.minor.patch with optional pre-release/build metadata
SEMVER_PATTERN = re.compile(
    r"^(0|[1-9]\d*)\.(0|[1-9]\d*)\.(0|[1-9]\d*)"
    r"(?:-((?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*)"
    r"(?:\.(?:0|[1-9]\d*|\d*[a-zA-Z-][0-9a-zA-Z-]*))*))?"
    r"(?:\+([0-9a-zA-Z-]+(?:\.[0-9a-zA-Z-]+)*))?$"
)


def validate_semver(version: str) -> bool:
    """Check if a version string is valid semver (X.Y.Z with optional pre-release/build)."""
    return bool(SEMVER_PATTERN.match(version))


def parse_name_version(identifier: str) -> tuple[str, str | None]:
    """Parse a 'name:version' identifier into (name, version).

    Args:
        identifier: A string like "fasta" or "fasta:0.1.0"

    Returns:
        Tuple of (name, version) where version is None if not specified.
    """
    if ":" in identifier:
        name, version = identifier.rsplit(":", 1)
        return name, version
    return identifier, None


def _semver_sort_key(version: str) -> tuple:
    """Create a sort key for semver comparison.

    Splits on '.' and converts numeric parts to ints for proper ordering.
    Pre-release versions (with '-') sort before the release.
    """
    # Strip build metadata
    base = version.split("+")[0]
    # Split pre-release
    parts = base.split("-", 1)
    core = parts[0]
    pre_release = parts[1] if len(parts) > 1 else None

    core_tuple = tuple(int(x) for x in core.split("."))
    # Releases sort after pre-releases: (core, 1, "") vs (core, 0, pre)
    if pre_release is None:
        return (core_tuple, 1, "")
    return (core_tuple, 0, pre_release)


def resolve_latest(session: Session, model_class, name: str):
    """Query the DB for the latest version of a given name using semver ordering.

    Works for any SQLModel with `name` and `version` fields (AssetClass, Recipe).

    Args:
        session: The database session.
        model_class: The SQLModel class (e.g., AssetClass or Recipe).
        name: The name to look up.

    Returns:
        The model instance with the latest version.

    Raises:
        NoResultFound: If no entries exist with that name.
    """
    results = session.exec(select(model_class).where(model_class.name == name)).unique().all()

    if not results:
        raise NoResultFound(f"No {model_class.__name__} found with name '{name}'")

    if len(results) == 1:
        return results[0]

    return max(results, key=lambda r: _semver_sort_key(r.version))
