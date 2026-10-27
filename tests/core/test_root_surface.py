"""``Refgenie`` stays a thin aggregate: build and channel sync live on managers.

Guards against the removed wrappers creeping back onto the root. Build is
``rgc.build``; data-channel sync is ``rgc.sources.sync_channel``.
"""

import pytest

from refgenie.managers.build import BuildManager

pytestmark = pytest.mark.unit

REMOVED = [
    "build_asset",
    "preflight_build",
    "initialize_and_build",
    "sync_data_channel",
    "resolve_custom_seek_keys",
    "resolve_default_asset",
    "get_asset_build_target_template",
    "get_genome_init_target_template",
    "_builder",
    "_asset_builder_instance",
    "_custom_seek_key_cache",
]


@pytest.mark.parametrize("name", REMOVED)
def test_root_has_no_build_or_sync_wrappers(refgenie_minimal, name):
    assert not hasattr(refgenie_minimal, name)


def test_build_is_the_public_build_manager(refgenie_minimal):
    assert isinstance(refgenie_minimal.build, BuildManager)
    assert refgenie_minimal.build is refgenie_minimal.build
    assert callable(refgenie_minimal.sources.sync_channel)
