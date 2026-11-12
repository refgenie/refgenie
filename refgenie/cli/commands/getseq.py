"""The `getseq` command: model and handler."""

import sys

from pydantic import AliasChoices, BaseModel, Field

from refgenie.cli.commands.helpers import (
    _GENOME_DESCRIPTION,
    _GENOME_DIGEST_DESCRIPTION,
    genome_ref,
)
from refgenie.cli.errors import EXIT_INVALID_INPUT, EXIT_NOT_FOUND, fail
from refgenie.exceptions import MissingSequenceError
from refgenie.managers.sequence import parse_locus


class GetseqModel(BaseModel):
    genome: str | None = Field(
        default=None,
        description=_GENOME_DESCRIPTION
        + " A genome not registered here is first registered from the subscribed "
        "servers (metadata only, like `pull --init`).",
        validation_alias=AliasChoices("g", "genome"),
    )
    genome_digest: str | None = Field(
        default=None,
        description=_GENOME_DIGEST_DESCRIPTION,
        validation_alias=AliasChoices("genome-digest"),
    )
    locus: str = Field(
        description="Coordinates of desired sequence, 0-based and half-open unless "
        "--one-based is given; e.g. 'chr1:50000-50200'.",
        validation_alias=AliasChoices("l", "locus"),
    )
    one_based: bool = Field(
        default=False,
        description="Read the locus as 1-based, closed coordinates (samtools style).",
        validation_alias=AliasChoices("one-based"),
    )
    fasta: bool = Field(
        default=False,
        description="Print a FASTA record: a header line, then 60 bases per line.",
        validation_alias=AliasChoices("fasta"),
    )


def handle_getseq(cmd, refgenie) -> None:
    ref = genome_ref(cmd.genome, cmd.genome_digest)
    digest = refgenie.genome.init_from_remote(ref, require_collection=True)
    if digest is None:
        fail(
            f"Genome '{ref}' is not registered here, and no subscribed server could provide it.",
            EXIT_NOT_FOUND,
        )
    try:
        name, start, end = parse_locus(cmd.locus, one_based=cmd.one_based)
    except ValueError as e:
        fail(str(e), EXIT_INVALID_INPUT)
    try:
        seq = refgenie.sequence.get(digest, name, start, end)
    except MissingSequenceError as e:
        fail(str(e), EXIT_NOT_FOUND)
    except ValueError as e:
        fail(str(e), EXIT_INVALID_INPUT)
    if cmd.fasta:
        lines = [seq[i : i + 60] for i in range(0, len(seq), 60)]
        sys.stdout.write(f">{cmd.locus}\n" + "\n".join(lines) + "\n")
    else:
        sys.stdout.write(seq + "\n")
