"""Add species tables

``species`` and ``speciesname`` hold the species in the stores' ``species.json``
files, loaded by ``store sync``: one row per taxon, and one row per name of it.
No foreign key to or from ``genome``.

Revision ID: 18dacacdcc92
Revises: 3e6f9e817a1f
Create Date: 2026-10-31 14:57:49.363028

"""

from typing import Sequence, Union

import sqlalchemy as sa
import sqlmodel
from alembic import op

# revision identifiers, used by Alembic.
revision: str = "18dacacdcc92"
down_revision: Union[str, None] = "3e6f9e817a1f"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        "species",
        sa.Column("taxon_id", sa.Integer(), autoincrement=False, nullable=False),
        sa.Column("scientific_name", sqlmodel.sql.sqltypes.AutoString(), nullable=False),
        sa.Column("common_name", sqlmodel.sql.sqltypes.AutoString(), nullable=True),
        sa.PrimaryKeyConstraint("taxon_id"),
    )
    op.create_table(
        "speciesname",
        sa.Column("taxon_id", sa.Integer(), nullable=False),
        sa.Column("name", sqlmodel.sql.sqltypes.AutoString(), nullable=False),
        sa.ForeignKeyConstraint(["taxon_id"], ["species.taxon_id"]),
        sa.PrimaryKeyConstraint("taxon_id", "name"),
    )
    op.create_index(op.f("ix_speciesname_name"), "speciesname", ["name"], unique=False)


def downgrade() -> None:
    op.drop_index(op.f("ix_speciesname_name"), table_name="speciesname")
    op.drop_table("speciesname")
    op.drop_table("species")
