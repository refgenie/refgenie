"""Alias coverage in the published catalog export.

A build node splits its aliases between its RefgetStore (genomes it built) and
the SQL ``alias`` table (genomes it federates over). The export must carry both
halves, and must fail loudly rather than publish a catalog whose genomes are
served as bare digests.
"""

import pytest
from sqlmodel import select

from refgenie import Refgenie
from refgenie.catalog_transfer import export_publish_catalog
from refgenie.db.tables import Alias, StoreType
from refgenie.managers.alias import AliasManager, StoreAliasManager

from .helpers import make_engine

pytestmark = pytest.mark.component

HTTPS_PREFIX = "https://example.org/assets"

#: A genome this node did not build. Only its SQL rows exist here.
FEDERATED_DIGEST = "F" * 32


@pytest.fixture
def build_node(tmp_path, fixtures_path):
    """A local-mode catalog holding one locally built and one federated genome."""
    rg = Refgenie(database_engine=make_engine(), suppress_migrations=True)
    rg.init(genome_folder=tmp_path / "genomes", genome_stage_folder=tmp_path / "stage")
    local_digest, _ = rg.genome.initialize_genome(
        fasta_file_path=fixtures_path / "rCRSd.fa",
        alias_names=["rCRSd"],
        description="built here",
    )

    rg.store.add("vgp", "url://vgp", StoreType.remote, priority=20)
    rg.genome.add(FEDERATED_DIGEST, "federated genome", [], store_name="vgp")
    AliasManager(rg.database_engine).federated_sync(
        "rAllMis1", FEDERATED_DIGEST, "vgp", {"vgp": 20}
    )
    return rg, local_digest


def _exported_aliases(path):
    from sqlalchemy import create_engine

    engine = create_engine(f"sqlite:///{path}")
    with engine.connect() as conn:
        return {(a.name, a.genome_digest) for a in conn.execute(select(Alias.__table__))}


class TestAliasCoverage:
    def test_export_carries_both_alias_halves(self, build_node, tmp_path):
        rg, local_digest = build_node
        dest = tmp_path / "publish.sqlite"

        summary = export_publish_catalog(rg, dest, HTTPS_PREFIX)

        assert summary["genome"] == 2
        assert summary["alias"] == 2
        assert _exported_aliases(dest) == {
            ("rCRSd", local_digest),
            ("rAllMis1", FEDERATED_DIGEST),
        }

    def test_export_refuses_to_drop_federated_names(self, build_node, tmp_path):
        """The guard fires when the alias manager reads only the store half.

        This is the exact regression it exists for: every table still has rows,
        so nothing else in the export notices.
        """
        rg, _local_digest = build_node
        store_only = StoreAliasManager(refget_store_getter=lambda: rg.refget_store)
        rg._alias_manager = store_only

        with pytest.raises(RuntimeError, match="alias rows the catalog holds"):
            export_publish_catalog(rg, tmp_path / "publish.sqlite", HTTPS_PREFIX)
