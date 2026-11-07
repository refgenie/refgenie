"""`seek -e` accepts prefix seek keys whose files exist (#296)."""

from pathlib import Path

import pytest

from refgenie.db.tables import SeekKeyType
from tests.support.factories import GENOME, GROUP, stage_copy
from tests.support.snapshots import fasta_asset

pytestmark = pytest.mark.component


def _add_with_prefix_key(r, prefix):
    _, rel = stage_copy(r, fasta_asset(r), "staging_prefix", mutate=b">p\nACGT\n")
    r.asset.content.add(
        asset_class_name="fasta",
        path=rel,
        asset_group_name=GROUP,
        genome_digest=r.alias.resolve(GENOME),
        asset_name="with_prefix",
        custom_seek_keys={"idx": (prefix, SeekKeyType.prefix)},
    )


@pytest.mark.parametrize("abs_path", [False, True])
def test_force_exists_accepts_a_prefix_with_files(refgenie_built, abs_path):
    r = refgenie_built
    digest = r.alias.resolve(GENOME)
    _add_with_prefix_key(r, digest)  # files: <digest>.fa, <digest>.fa.fai, ...
    path = r.asset.seek(digest, GROUP, "with_prefix", "idx", force_exists=True, abs_path=abs_path)
    assert not Path(path).exists()  # the prefix itself is not a file
    assert list(Path(path).parent.glob(f"{Path(path).name}*"))


def test_force_exists_rejects_a_prefix_without_files(refgenie_built):
    r = refgenie_built
    _add_with_prefix_key(r, "no_such_prefix")
    with pytest.raises(FileNotFoundError):
        r.asset.seek(r.alias.resolve(GENOME), GROUP, "with_prefix", "idx", force_exists=True)
