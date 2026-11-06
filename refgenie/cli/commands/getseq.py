"""The `getseq` command: model and handler."""

from pydantic import AliasChoices, BaseModel, Field
from rich import print as rprint

from refgenie.cli.commands.helpers import (
    _GENOME_DESCRIPTION,
    _GENOME_DIGEST_DESCRIPTION,
    genome_digest_from,
)


class GetseqModel(BaseModel):
    genome: str | None = Field(
        default=None,
        description=_GENOME_DESCRIPTION,
        validation_alias=AliasChoices("g", "genome"),
    )
    genome_digest: str | None = Field(
        default=None,
        description=_GENOME_DIGEST_DESCRIPTION,
        validation_alias=AliasChoices("genome-digest"),
    )
    locus: str = Field(
        description="Coordinates of desired sequence (0-based, half-open); "
        "e.g. 'chr1:50000-50200'.",
        validation_alias=AliasChoices("l", "locus"),
    )


def handle_getseq(cmd, refgenie) -> None:
    sequence = refgenie.sequence.get(
        genome_digest=genome_digest_from(refgenie, cmd.genome, cmd.genome_digest),
        locus=cmd.locus,
    )
    rprint(sequence)
