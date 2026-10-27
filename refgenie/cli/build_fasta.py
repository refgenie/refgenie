"""CLI command to build fasta asset files from RefgetStore.

Produces .fa, .fa.fai, and .chrom.sizes from a RefgetStore collection.
Designed to be called from recipe shell templates. No samtools required.
"""

import argparse
from pathlib import Path


def build_fasta_main(args=None):
    """Build fasta asset files from a RefgetStore."""
    parser = argparse.ArgumentParser(
        description="Build fasta asset files (.fa, .fai, .chrom.sizes) from RefgetStore"
    )
    parser.add_argument("store_path", type=Path, help="Path to the RefgetStore directory")
    parser.add_argument("digest", help="Sequence collection digest to export")
    parser.add_argument("output_folder", type=Path, help="Output folder for fasta asset files")
    parser.add_argument(
        "--line-width", type=int, default=80, help="Line width for FASTA output (default: 80)"
    )
    parsed = parser.parse_args(args)

    from gtars.refget import RefgetStore, compute_fai

    store = RefgetStore.open_local(parsed.store_path)
    digest = parsed.digest
    # Load only the collection being exported: the store opens as lazy stubs
    # but export_fasta() needs Full records, and load_all_*() makes memory
    # scale with the whole store rather than the genome.
    store.load_collection(digest)
    for seq_digest in store.get_collection_level2(digest)["sequences"]:
        store.load_sequence(seq_digest.removeprefix("SQ."))
    output = parsed.output_folder
    output.mkdir(parents=True, exist_ok=True)

    fa_path = output / f"{digest}.fa"
    store.export_fasta(digest, fa_path, None, parsed.line_width)

    # Compute FAI index (no samtools needed)
    fai_records = compute_fai(fa_path)
    fai_path = output / f"{digest}.fa.fai"
    with open(fai_path, "w") as f:
        for record in fai_records:
            if record.fai is None:
                raise RuntimeError(
                    f"No FAI metadata for sequence '{record.name}' -- is the FASTA compressed?"
                )
            f.write(
                f"{record.name}\t{record.length}\t{record.fai.offset}\t"
                f"{record.fai.line_bases}\t{record.fai.line_bytes}\n"
            )

    level2 = store.get_collection_level2(digest)
    chrom_sizes_path = output / f"{digest}.chrom.sizes"
    with open(chrom_sizes_path, "w") as f:
        for name, length in zip(level2["names"], level2["lengths"]):
            f.write(f"{name}\t{length}\n")


if __name__ == "__main__":
    build_fasta_main()
