"""The `getseq` command: model and handler."""

import sys

from pydantic import AliasChoices, BaseModel, Field

from refgenie.cli.commands.helpers import (
    _GENOME_DESCRIPTION,
    _GENOME_DIGEST_DESCRIPTION,
    genome_ref,
    note,
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
        description="The region, 1-based and inclusive like samtools and genome "
        "browsers; e.g. 'chr17:43,044,295-43,044,394' (commas are fine), or a "
        "whole sequence, 'chrM'.",
        validation_alias=AliasChoices("l", "locus"),
    )
    zero_based: bool = Field(
        default=False,
        description="Read the region's numbers as 0-based and half-open (BED style), "
        "as the HTTP API, MCP and Python API take them.",
        validation_alias=AliasChoices("zero-based"),
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
        name, start, end = parse_locus(cmd.locus, zero_based=cmd.zero_based)
    except ValueError as e:
        fail(str(e), EXIT_INVALID_INPUT)
    try:
        seq = refgenie.sequence.get(digest, name, start, end)
    except MissingSequenceError as e:
        fail(str(e), EXIT_NOT_FOUND)
    except ValueError as e:
        fail(str(e), EXIT_INVALID_INPUT)
    # Say how the region was read, 1-based like the text people write, so an
    # off-by-one is visible.
    first = (start or 0) + 1
    region = f"{name}:{first}-{first + len(seq) - 1}" if start is not None or end else name
    note(f"{name}:{first:,}-{first + len(seq) - 1:,} (1-based, inclusive): {len(seq):,} bases")
    if cmd.fasta:
        lines = [seq[i : i + 60] for i in range(0, len(seq), 60)]
        sys.stdout.write(f">{region}\n" + "\n".join(lines) + "\n")
    else:
        sys.stdout.write(seq + "\n")
