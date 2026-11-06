"""The one exception classifier, ``refgenie.server.errors.classify_exception``."""

import pytest

from tests.support.server import requires_server

requires_server()

from refgenie.exceptions import (  # noqa: E402  (must follow the extras guard)
    AssetExistsError,
    MissingAliasError,
    MissingAssetError,
    MissingBuildInputError,
    MissingGenomeError,
    MissingRecipeError,
    NoArchiveError,
    NoSubscriptionsError,
    PullFailedError,
    PullSkipped,
    RefgenieError,
)
from refgenie.server.errors import (  # noqa: E402
    _CLASSIFICATION,
    ErrorCode,
    WebError,
    action_error,
    classify_exception,
)


@pytest.mark.parametrize(
    "exc, status, code",
    [
        # Subclasses of PullFailedError keep their own codes.
        (AssetExistsError("asset already there"), 409, ErrorCode.ASSET_EXISTS),
        (NoArchiveError("no archive on server"), 404, ErrorCode.NO_ARCHIVE),
        (NoSubscriptionsError(), 409, ErrorCode.NO_SUBSCRIPTIONS),
        (PullFailedError("download died"), 502, ErrorCode.PULL_FAILED),
        # So do the RefgenieError lookups.
        (MissingGenomeError(genome="x"), 404, ErrorCode.GENOME_NOT_FOUND),
        (MissingAliasError("hg38"), 404, ErrorCode.ALIAS_NOT_FOUND),
        (MissingAssetError(digest="abc"), 404, ErrorCode.ASSET_NOT_FOUND),
        (MissingRecipeError("fasta"), 404, ErrorCode.RECIPE_NOT_FOUND),
        (MissingBuildInputError("need a file"), 400, ErrorCode.MISSING_BUILD_INPUT),
        (PullSkipped("declined"), 409, ErrorCode.PULL_SKIPPED),
        (RefgenieError("something domain-y"), 500, ErrorCode.REFGENIE_ERROR),
        (ValueError("asset has children"), 409, ErrorCode.CONFLICT),
        (KeyError("nope"), 500, ErrorCode.INTERNAL_ERROR),
        (WebError("pipeline failed", ErrorCode.BUILD_FAILED, 500), 500, ErrorCode.BUILD_FAILED),
    ],
    ids=lambda v: type(v).__name__ if isinstance(v, Exception) else None,
)
def test_classification(exc, status, code):
    assert classify_exception(exc)[:2] == (status, code)


def test_messages_pass_through():
    assert "subscribe" in classify_exception(NoSubscriptionsError())[2].lower()
    exc = WebError("no servers", ErrorCode.NO_SUBSCRIPTIONS, 409)
    assert classify_exception(exc) == (409, ErrorCode.NO_SUBSCRIPTIONS, "no servers")


def test_no_subclass_is_listed_after_its_parent():
    """The table is scanned in order, so a subclass after its parent never matches."""
    types = [entry[0] for entry in _CLASSIFICATION]
    for i, earlier in enumerate(types):
        for later in types[i + 1 :]:
            if later is not earlier:
                assert not issubclass(later, earlier), (
                    f"{later.__name__} is listed after its parent {earlier.__name__}"
                )


def test_action_error_carries_classified_status_and_envelope_detail():
    exc = action_error(MissingGenomeError(genome="abc123"))
    assert exc.status_code == 404
    assert exc.detail["code"] == str(ErrorCode.GENOME_NOT_FOUND)
    assert "abc123" in exc.detail["message"]
