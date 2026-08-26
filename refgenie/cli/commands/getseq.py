"""The `getseq` command: model and handler."""

from pydantic import AliasChoices, BaseModel, Field
from rich import print as rprint

from refgenie.cli.commands.helpers import _GENOME_DESCRIPTION


class GetseqModel(BaseModel):
    genome: str = Field(
        description=_GENOME_DESCRIPTION,
        validation_alias=AliasChoices("g", "genome"),
    )
    locus: str = Field(
        description="Coordinates of desired sequence (0-based, half-open); "
        "e.g. 'chr1:50000-50200'.",
        validation_alias=AliasChoices("l", "locus"),
    )


def handle_getseq(cmd, refgenie) -> None:
    sequence = refgenie.getseq(
        genome_name=cmd.genome,
        locus=cmd.locus,
    )
    rprint(sequence)
