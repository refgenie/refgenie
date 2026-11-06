"""
Every SQL statement refgenie runs, as named functions, one module per entity.

Each function takes the session first and never opens, commits or rolls one
back; the manager that called it owns the transaction. It returns rows, counts,
or ``(rows, total)`` for a page, never a statement. Every list has a fixed
``ORDER BY`` ending in a unique column, and every function eager-loads what its
callers read (picked by boolean keywords such as ``with_links=True``), so rows
stay usable after the session closes.
"""
