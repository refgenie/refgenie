"""Replace NCBI taxonomy tables with genome synonyms

Drops ``taxonname`` and ``taxonlineage``, which were never filled, and adds
``genome.synonyms``: other names of the genome, read from FHR ``genomeSynonym``.

Revision ID: 3e6f9e817a1f
Revises: d5d9db835bcd
Create Date: 2026-10-30 21:38:42.017875

"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa
import sqlmodel


# revision identifiers, used by Alembic.
revision: str = "3e6f9e817a1f"
down_revision: Union[str, None] = "d5d9db835bcd"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.add_column("genome", sa.Column("synonyms", sa.JSON(), nullable=True))
    op.drop_index(op.f("ix_taxonname_name"), table_name="taxonname")
    op.drop_table("taxonname")
    op.drop_index(op.f("ix_taxonlineage_ancestor_id"), table_name="taxonlineage")
    op.drop_table("taxonlineage")


def downgrade() -> None:
    op.create_table(
        "taxonlineage",
        sa.Column("taxon_id", sa.Integer(), nullable=False),
        sa.Column("ancestor_id", sa.Integer(), nullable=False),
        sa.PrimaryKeyConstraint("taxon_id", "ancestor_id"),
    )
    op.create_index(
        op.f("ix_taxonlineage_ancestor_id"), "taxonlineage", ["ancestor_id"], unique=False
    )
    op.create_table(
        "taxonname",
        sa.Column("taxon_id", sa.Integer(), nullable=False),
        sa.Column("name", sqlmodel.sql.sqltypes.AutoString(), nullable=False),
        sa.PrimaryKeyConstraint("taxon_id", "name"),
    )
    op.create_index(op.f("ix_taxonname_name"), "taxonname", ["name"], unique=False)
    op.drop_column("genome", "synonyms")
