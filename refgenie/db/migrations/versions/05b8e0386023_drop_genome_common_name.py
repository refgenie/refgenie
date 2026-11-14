"""Drop genome.common_name

Common names come from the species list (``species.common_name``), looked up
by ``taxon_id``. FHR ``commonName`` is no longer read, so nothing fills the
column.

Revision ID: 05b8e0386023
Revises: 18dacacdcc92
Create Date: 2026-10-31 15:42:58

"""

from typing import Sequence, Union

import sqlalchemy as sa
import sqlmodel
from alembic import op

# revision identifiers, used by Alembic.
revision: str = "05b8e0386023"
down_revision: Union[str, None] = "18dacacdcc92"
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.drop_column("genome", "common_name")


def downgrade() -> None:
    op.add_column(
        "genome",
        sa.Column("common_name", sqlmodel.sql.sqltypes.AutoString(), nullable=True),
    )
