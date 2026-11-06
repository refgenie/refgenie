"""Counting, ordering and slicing a list statement in one place."""

from typing import Any

from sqlalchemy import Select
from sqlmodel import Session, func, select


def paginate(
    session: Session,
    statement: Select[Any],
    *,
    order_by: tuple[Any, ...],
    offset: int,
    limit: int,
) -> tuple[list[Any], int]:
    """
    Run ``statement`` as one page and count every row it matches.

    Args:
        session: The session to execute in.
        statement: The filtered SELECT, without ordering or slicing.
        order_by: The ORDER BY columns. Required, and must end in a unique
            column, so consecutive pages neither overlap nor skip rows.
        offset: Rows to skip.
        limit: Most rows to return.

    Returns:
        The page's rows and the total number of matching rows.
    """
    total = session.exec(select(func.count()).select_from(statement.subquery())).one()
    rows = session.exec(statement.order_by(*order_by).offset(offset).limit(limit)).unique().all()
    return list(rows), total
