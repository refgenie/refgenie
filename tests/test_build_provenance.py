"""Tests for refgenie.utils.build's build-digest algorithm, which identifies
*a build* as opposed to its output.

``Asset.digest`` addresses the bytes a build produced; ``build_digest``
addresses the build. The failure worth guarding is a digest that ignores a
load-bearing input, which makes two different builds look like one.
"""

from datetime import datetime

import pytest

from refgenie.utils.build import build_digest, build_level1_to_digest, build_level2

pytestmark = pytest.mark.unit


#: The arguments every test perturbs one attribute of.
BASE = dict(
    genome_digest="genome-A",
    recipe_name="fasta",
    recipe_version="0.0.1",
    input_asset_digests={"parent": "digest-1"},
    input_file_digests={"source": "sha-1"},
    params={"threads": "4"},
)


def digest_of(**overrides) -> str:
    return build_digest(build_level2(**{**BASE, **overrides}))[0]


@pytest.mark.parametrize(
    "override",
    [
        pytest.param({"genome_digest": "genome-B"}, id="genome"),
        pytest.param({"recipe_name": "bwa_index"}, id="recipe"),
        pytest.param({"recipe_version": "0.0.2"}, id="recipe_version"),
        pytest.param({"input_asset_digests": {"parent": "digest-2"}}, id="input_asset"),
        pytest.param({"input_file_digests": {"source": "sha-2"}}, id="input_file"),
        pytest.param({"params": {"threads": "8"}}, id="param"),
    ],
)
def test_every_identity_attribute_changes_the_digest(override):
    assert digest_of(**override) != digest_of()


def test_is_reproducible_from_its_level1_object():
    """A stored row verifies without re-reading the files the build consumed."""
    digest, level1 = build_digest(build_level2(**BASE))
    assert build_level1_to_digest(level1) == digest


def test_recipes_with_no_declared_inputs_still_differ_by_genome():
    """The regression guard for the `fasta` recipe.

    It declares no input files, params or assets -- it reads the RefgetStore by
    genome digest. Drop the genome from the digest and every genome's fasta
    build computes the same value.
    """
    empty = dict(
        recipe_name="fasta",
        recipe_version="0.0.1",
        input_asset_digests={},
        input_file_digests={},
        params={},
    )
    assert (
        build_digest(build_level2(genome_digest="genome-A", **empty))[0]
        != build_digest(build_level2(genome_digest="genome-B", **empty))[0]
    )


def test_non_inherent_attributes_do_not_move_the_digest():
    """When a build ran, and under which version or image, is not what it built.

    The input-asset *names* are non-inherent for the same reason, and that is
    what makes a row's digest reconstructible from AssetLink, which stores none.
    """
    varied = build_level2(
        **BASE,
        build_timestamp=datetime(2026, 1, 1),
        refgenie_version="9.9.9",
        docker_image="my/image",
    )
    assert build_digest(build_level2(**BASE))[0] == build_digest(varied)[0]
    assert digest_of(input_asset_digests={"first": "d1"}) == digest_of(
        input_asset_digests={"second": "d1"}
    )


def test_param_type_does_not_change_the_digest():
    """Params are interpolated as text into the recipe's shell command, so 10
    and "10" are one build. A float must digest rather than raise: canonical_str
    rejects floats, and the coercion is what keeps params inside its domain.
    """
    assert digest_of(params={"size": 10}) == digest_of(params={"size": "10"})
    assert digest_of(params={"fraction": 0.5})
