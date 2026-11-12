"""Queries on the SQL ``alias`` table (server mode, and federated stores)."""

from sqlmodel import Session, delete, select

from refgenie.db.queries.common import one_or_raise
from refgenie.db.tables import Alias
from refgenie.exceptions import MissingAliasError


def get(session: Session, name: str) -> Alias:
    """
    The alias row with this exact name.

    Raises:
        MissingAliasError: If there is none.
    """
    return one_or_raise(session, select(Alias).where(Alias.name == name), MissingAliasError(name))


def find(session: Session, name: str) -> Alias | None:
    """The alias row with this exact name, or None."""
    return session.exec(select(Alias).where(Alias.name == name)).one_or_none()


def exists(session: Session, name: str) -> bool:
    """Whether an alias with this exact name exists."""
    return session.exec(select(Alias.name).where(Alias.name == name)).first() is not None


def list_all(session: Session, genome_digest: str | None = None) -> list[Alias]:
    """Every alias, or every alias of one genome, by name."""
    statement = select(Alias)
    if genome_digest:
        statement = statement.where(Alias.genome_digest == genome_digest)
    return list(session.exec(statement.order_by(Alias.name)).all())


def list_by_store(session: Session, store_name: str) -> list[Alias]:
    """The alias rows a federated store wrote (``store_name``), collision losers included."""
    return list(session.exec(select(Alias).where(Alias.store_name == store_name)).all())


def delete_names(session: Session, names: list[str]) -> None:
    """Delete these alias rows. A bulk delete: no per-row events, so no files are touched."""
    if names:
        session.exec(delete(Alias).where(Alias.name.in_(names)))


def names_by_genome(session: Session, digests: list[str] | None = None) -> dict[str, list[str]]:
    """
    Alias names keyed by genome digest, in one query.

    Args:
        session: The session to execute in.
        digests: Only these genomes. None reads every alias.

    Returns:
        Names by digest, each list sorted. A genome with no alias is absent.
    """
    statement = select(Alias.name, Alias.genome_digest).order_by(Alias.name)
    if digests is not None:
        if not digests:
            return {}
        statement = statement.where(Alias.genome_digest.in_(digests))
    by_digest: dict[str, list[str]] = {}
    for name, digest in session.exec(statement).all():
        by_digest.setdefault(digest, []).append(name)
    return by_digest
