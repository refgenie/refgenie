"""The hooks refgenie fires and the objects a plugin receives.

Standard library only: a plugin imports these types, and importing them must
cost nothing.
"""

from dataclasses import dataclass, field
from typing import Literal

#: Entry-point groups are named ``refgenie.hooks.<hook>``.
ENTRY_POINT_PREFIX = "refgenie.hooks."

PRE_PULL = "pre_pull"
POST_PULL = "post_pull"
PRE_BUILD = "pre_build"
POST_BUILD = "post_build"
POST_UPDATE = "post_update"
HOOKS: tuple[str, ...] = (PRE_PULL, POST_PULL, PRE_BUILD, POST_BUILD, POST_UPDATE)

ChangeAction = Literal[
    "asset_added",
    "asset_removed",
    "asset_renamed",
    "default_changed",
    "genome_added",
    "genome_removed",
    "alias_added",
    "alias_removed",
]


@dataclass(frozen=True)
class Change:
    """One committed change to local state. Carried by ``post_update``."""

    action: ChangeAction
    #: The genome digest.
    genome: str | None = None
    asset_group: str | None = None
    #: The asset name (the new name, for a rename).
    asset: str | None = None
    #: The old asset name for a rename; the old default for ``default_changed``.
    previous: str | None = None
    #: The asset content digest, when known.
    digest: str | None = None
    #: The alias, for ``alias_added`` / ``alias_removed``.
    alias: str | None = None


@dataclass(frozen=True)
class HookEvent:
    """What a plugin receives as its second argument.

    New fields may be added later; a plugin should read only the ones it needs.
    """

    hook: str
    #: The genome as the caller named it (alias or digest).
    genome: str | None = None
    asset_group: str | None = None
    asset: str | None = None
    #: ``post_pull`` / ``post_build`` only.
    succeeded: bool | None = None
    #: ``post_update`` only.
    changes: tuple[Change, ...] = field(default_factory=tuple)
