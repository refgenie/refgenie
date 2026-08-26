"""
One way to build the Rich tables the managers print.

Every ``.table()`` method in ``refgenie.managers`` used to open a `Table`, call
``add_column`` once per heading and ``add_row`` once per record, which is the
same eleven lines eight times over with only the headings differing. `build_table`
takes the headings and the rows and does the rest, so a manager's table method
is left holding only the query and the row shaping that are actually its own.

Anything Rich itself accepts on the constructor (``caption``, ``box``, ...)
passes straight through as a keyword argument.
"""

from collections.abc import Iterable, Sequence
from typing import Any

from rich.table import Table

#: Yield this in place of a row to draw a separator before the next one.
SECTION = object()


def build_table(
    title: str,
    columns: Sequence[str | tuple[str, str]],
    rows: Iterable[Sequence[Any] | object],
    *,
    end_section: bool = False,
    **table_kwargs: Any,
) -> Table:
    """
    Build a Rich table from headings and rows.

    Args:
        title: The table title.
        columns: One heading per column, either a name or a ``(name, style)``
            pair.
        rows: The rows, each a sequence of cell values with one entry per
            column. A row that is the `SECTION` sentinel draws a separator
            instead.
        end_section: Draw a separator after every row. Use for tables whose
            rows are tall enough to run together otherwise.
        **table_kwargs: Passed to `rich.table.Table`.

    Returns:
        Table: The built table.
    """
    table = Table(title=title, **table_kwargs)
    for column in columns:
        if isinstance(column, tuple):
            name, style = column
            table.add_column(name, style=style)
        else:
            table.add_column(column)

    for row in rows:
        if row is SECTION:
            table.add_section()
        else:
            table.add_row(*row, end_section=end_section)
    return table
