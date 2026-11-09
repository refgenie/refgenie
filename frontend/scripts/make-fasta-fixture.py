"""Build the refget store fixture for src/services/fastaExport.test.ts.

Run from the repository root:

    uv run --with gtars python frontend/scripts/make-fasta-fixture.py

Writes src/test/fixtures/fasta-store/{source.fa,digest.txt,store/}. The store
is checked to export back to source.fa byte for byte with native gtars, so the
browser converter test compares against a known-good FASTA.
"""

import random
import shutil
from pathlib import Path

from gtars.refget import RefgetStore

FIXTURE = Path(__file__).resolve().parent.parent / "src" / "test" / "fixtures" / "fasta-store"
LINE_WIDTH = 80


def bases(rng: random.Random, alphabet: str, n: int) -> str:
    return "".join(rng.choice(alphabet) for _ in range(n))


def records() -> list[tuple[str, str]]:
    rng = random.Random(7)
    chr_b = (
        bases(rng, "ACGT", 40)
        + "N" * 30
        + bases(rng, "ACGT", 50)
        + "N" * 7
        + bases(rng, "ACGT", 33)
    )
    return [
        ("chrA", bases(rng, "ACGT", 203)),
        ("chrB test contig B", chr_b),
        ("chrC", "ARYKMGT"),
        ("chrD long one", bases(rng, "ACGT", 1000)),
        ("chrE", bases(rng, "ACGTN", 80) + "N"),
    ]


def write_fasta(path: Path, recs: list[tuple[str, str]]) -> None:
    with path.open("w") as out:
        for header, seq in recs:
            out.write(f">{header}\n")
            for i in range(0, len(seq), LINE_WIDTH):
                out.write(seq[i : i + LINE_WIDTH] + "\n")


def main() -> None:
    if FIXTURE.exists():
        shutil.rmtree(FIXTURE)
    FIXTURE.mkdir(parents=True)
    source = FIXTURE / "source.fa"
    recs = records()
    assert len(recs[1][1]) == 160 and len(recs[4][1]) == 81
    write_fasta(source, recs)

    store = RefgetStore.on_disk(str(FIXTURE / "store"))
    meta, _ = store.add_sequence_collection_from_fasta(str(source))
    digest = meta.digest
    (FIXTURE / "digest.txt").write_text(digest + "\n")

    native = FIXTURE / "native.fa"
    store.export_fasta(digest, str(native))
    assert native.read_bytes() == source.read_bytes(), "native export differs from source.fa"
    native.unlink()
    # gtars caches a digest sidecar next to the input FASTA.
    (FIXTURE / "source.rgsi").unlink(missing_ok=True)

    for leftover in (FIXTURE / "store").rglob("*"):
        temp = leftover.suffix in {".lock", ".tmp"} or leftover.name.startswith(".")
        if leftover.is_file() and temp:
            leftover.unlink()

    rgsi = FIXTURE / "store" / "collections" / f"{digest}.rgsi"
    lines = [ln for ln in rgsi.read_text().splitlines() if not ln.startswith("#")]
    print(digest, [ln.split("\t")[2] for ln in lines])


if __name__ == "__main__":
    main()
