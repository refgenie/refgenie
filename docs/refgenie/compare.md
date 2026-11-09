# Genome compatibility

Many analyses need two genomes to be compatible, not identical. For example, to annotate reads aligned with a `bowtie2_index` using a feature annotation file, the two only need the same coordinate system (sequence names and lengths); the sequences themselves may differ.

`refgenie compare` answers this. When you run `refgenie genome init`, refgenie loads the genome into the RefgetStore as a GA4GH [sequence collection](https://ga4gh.github.io/seqcol-spec/) (seqcol): its sequence names, lengths, and sequence digests. `refgenie compare` runs the standard seqcol comparison on two of these collections, so no asset has to be built first.

## Usage

To compare two initialized genomes, give their aliases:

```console
refgenie compare hg38 hg38_plus
```

You can also name genomes by digest, or mix an alias with a digest. Either way, you must give exactly two genomes:

```console
refgenie compare --genome-digest <digest_a> --genome-digest <digest_b>
refgenie compare hg38 --genome-digest <digest_b>
```

The result is a seqcol comparison object, for example:

```console
{
    'digests': {'a': 'pGauuPa6DYFFjNGTjVBbPsYgAmS5f3y2', 'b': 'hP1GNDAYRYvMmuEBp94uTC10v9W406oS'},
    'attributes': {
        'a_only': [],
        'b_only': [],
        'a_and_b': ['lengths', 'name_length_pairs', 'names', 'sequences', 'sorted_name_length_pairs', 'sorted_sequences']
    },
    'array_elements': {
        'a_count': {'names': 2, 'lengths': 2, 'sequences': 2, 'name_length_pairs': 2, ...},
        'b_count': {'names': 3, 'lengths': 3, 'sequences': 3, 'name_length_pairs': 3, ...},
        'a_and_b_count': {'names': 2, 'lengths': 2, 'sequences': 2, 'name_length_pairs': 2, ...},
        'a_and_b_same_order': {'names': True, 'lengths': None, 'sequences': True, 'name_length_pairs': True, ...}
    }
}
```

How to read it:

- `attributes` lists which attributes (such as `names`, `lengths`, `sequences`) exist in both collections, or in only one.
- `array_elements` counts the elements of each attribute: how many are in `a` (`a_count`), in `b` (`b_count`), and in both (`a_and_b_count`).
- `a_and_b_same_order` tells you whether the shared elements appear in the same order in both genomes. It is `None` when order cannot be judged, for example when an attribute has repeated values.

Based on the output above we can conclude that genome `hg38_plus` is a superset of `hg38`: every name, length, and sequence of `hg38` is also in `hg38_plus` (`a_and_b_count` equals `a_count`), in the same order, and `hg38_plus` has one extra sequence.

For the coordinate-system question above, look at `name_length_pairs`: if every name-length pair of one genome is shared with the other, the two share a coordinate system for those sequences, even if the sequences themselves differ.

The result has this same layout in every case: when both genomes are in one RefgetStore, when they live in different RefgetStores, and when one or both come from a remote server.

## Other ways to compare

- **Server**: a server-mode instance (`refgenie serve`) mounts the GA4GH seqcol service at `/seqcol`, which offers the standard comparison endpoint `GET /seqcol/comparison/{digest_a}/{digest_b}`.
- **AI assistants**: the [MCP server](mcp.md) provides the same comparison as the `compare_genomes` tool.
