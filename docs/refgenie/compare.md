# Genome compatibility

Many analyses need two genomes to be compatible, not identical. For example, to annotate reads aligned with a `bowtie2_index` using a feature annotation file, the two only need the same coordinate system (sequence names and lengths); the sequences themselves may differ.

`refgenie compare` answers this. When you run `refgenie genome init`, refgenie loads the genome into the RefgetStore as a GA4GH [sequence collection](https://ga4gh.github.io/seqcol-spec/) (seqcol): its sequence names, lengths, and sequence digests. `refgenie compare` runs the standard seqcol comparison on two of these collections, so no asset has to be built first.

## Usage

Give any two genomes, local or on your subscribed servers. A genome you have not pulled is compared on the server that has it:

```console
refgenie compare hg38 hg38_primary
```

You can also name genomes by digest, or mix an alias with a digest. Either way, you must give exactly two genomes:

```console
refgenie compare --genome-digest <digest_a> --genome-digest <digest_b>
refgenie compare hg38 --genome-digest <digest_b>
```

`refgenie compare` prints one sentence on how the two genomes differ, then every sequence that has a different name in the second genome:

```console
$ refgenie compare hg38 hg38_primary
Every name and length of hg38_primary is also in hg38, but 18 of those 25 sequences differ in content.

$ refgenie compare sacCer3 R64-1-1
Same 17 sequences in the same order; all 17 names differ: chrI->I, chrII->II, chrIII->III, ..., chrXVI->XVI.
  chrI -> I
  chrII -> II
  ...
  chrM -> Mito
```

The first pair shares a coordinate system (the same names and lengths for hg38_primary's 25 sequences), but the bases differ, so they are not interchangeable for alignment. The second pair holds the same sequences under different names, so files made for one can be used with the other after renaming.

## The full comparison

`refgenie compare --json` prints the summary, the rename map (`name_map`) and the GA4GH [sequence collection](https://ga4gh.github.io/seqcol-spec/) comparison (`comparison`):

```console
{
  "summary": "Same 17 sequences in the same order; all 17 names differ: ...",
  "name_map": {"chrI": "I", "chrII": "II", "...": "...", "chrM": "Mito"},
  "comparison": {
    "digests": {"a": "8zSJ7M2oerPcCKc8x310hjyfB5u8SJNJ", "b": "ia9g4Myony-xxLsSc7eQOC_NdClGgPGu"},
    "attributes": {"a_only": [], "b_only": [], "a_and_b": ["lengths", "names", "sequences", "..."]},
    "array_elements": {
      "a_count": {"names": 17, "lengths": 17, "sequences": 17, "...": 17},
      "b_count": {"names": 17, "lengths": 17, "sequences": 17, "...": 17},
      "a_and_b_count": {"names": 0, "lengths": 17, "sequences": 17, "...": 0},
      "a_and_b_same_order": {"names": null, "lengths": true, "sequences": true, "...": null}
    }
  }
}
```

How to read `comparison`:

- `attributes` lists which attributes (such as `names`, `lengths`, `sequences`) exist in both collections, or in only one.
- `array_elements` counts the elements of each attribute: how many are in `a` (`a_count`), in `b` (`b_count`), and in both (`a_and_b_count`).
- `a_and_b_same_order` tells you whether the shared elements appear in the same order in both genomes. It is `null` when order cannot be judged, for example when nothing is shared.

For the coordinate-system question above, look at `name_length_pairs`: if every name-length pair of one genome is shared with the other, the two share a coordinate system for those sequences, even if the sequences themselves differ.

## Other ways to compare

- **Website**: the Compare tab on a genome's page at [refgenie.org](https://refgenie.org) shows the same summary, against another genome of the species or your own FASTA file.
- **HTTP API**: `GET /v4/genomes/{digest}/compare/{other}` returns the GA4GH comparison plus `summary` and `name_map`, for any two genomes on the server. Each store is also a GA4GH seqcol service (`/seqcol`, or `/seqcol/<store>`; `/service-info` lists them) with the standard `GET .../comparison/{digest_a}/{digest_b}` for two genomes in that store.
- **AI assistants**: the [MCP server](mcp.md) provides the same comparison as the `compare_genomes` tool.
