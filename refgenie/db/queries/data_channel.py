"""Queries on data channels."""

from sqlmodel import Session, select

from refgenie.db.tables import DataChannel


def find(session: Session, name: str) -> DataChannel | None:
    """The data channel with this name, or None."""
    return session.exec(select(DataChannel).where(DataChannel.name == name)).first()


def list_all(session: Session) -> list[DataChannel]:
    """Every data channel, by name."""
    return list(session.exec(select(DataChannel).order_by(DataChannel.name)).unique().all())
