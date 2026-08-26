"""Federated multi-store serving

Adds the ``store`` registry table (the single source of truth for what the
service serves) plus ``store_name`` owner pointers on ``genome`` and ``alias``,
so a single refgenie service can federate over several physically separate
refget stores. See ``refgenie/core/store_router.py``.

Revision ID: b2f1c3d4e5a6
Revises: 1a0e5c7b9d21
Create Date: 2026-08-25 08:10:00.000000

"""

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa
import sqlmodel


# revision identifiers, used by Alembic.
revision: str = "b2f1c3d4e5a6"
down_revision: Union[str, None] = "1a0e5c7b9d21"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        "store",
        sa.Column("name", sqlmodel.sql.sqltypes.AutoString(), nullable=False),
        sa.Column("url", sqlmodel.sql.sqltypes.AutoString(), nullable=False),
        sa.Column("type", sqlmodel.sql.sqltypes.AutoString(), nullable=False),
        sa.Column("priority", sa.Integer(), nullable=False),
        sa.Column("enabled", sa.Boolean(), nullable=False),
        sa.Column("description", sqlmodel.sql.sqltypes.AutoString(), nullable=True),
        sa.Column("updated_at", sa.DateTime(), nullable=True),
        sa.Column("created_at", sa.DateTime(), nullable=True),
        sa.PrimaryKeyConstraint("name"),
    )

    with op.batch_alter_table("genome", schema=None) as batch_op:
        batch_op.add_column(sa.Column("store_name", sqlmodel.sql.sqltypes.AutoString(), nullable=True))
        batch_op.create_index("ix_genome_store_name", ["store_name"], unique=False)
        batch_op.create_foreign_key("fk_genome_store_name", "store", ["store_name"], ["name"])

    with op.batch_alter_table("alias", schema=None) as batch_op:
        batch_op.add_column(sa.Column("store_name", sqlmodel.sql.sqltypes.AutoString(), nullable=True))
        batch_op.create_index("ix_alias_store_name", ["store_name"], unique=False)
        batch_op.create_foreign_key("fk_alias_store_name", "store", ["store_name"], ["name"])


def downgrade() -> None:
    with op.batch_alter_table("alias", schema=None) as batch_op:
        batch_op.drop_constraint("fk_alias_store_name", type_="foreignkey")
        batch_op.drop_index("ix_alias_store_name")
        batch_op.drop_column("store_name")

    with op.batch_alter_table("genome", schema=None) as batch_op:
        batch_op.drop_constraint("fk_genome_store_name", type_="foreignkey")
        batch_op.drop_index("ix_genome_store_name")
        batch_op.drop_column("store_name")

    op.drop_table("store")
