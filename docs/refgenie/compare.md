# Genomes compatibility

## Motivation

Many genomic analyses require less stringent comparison: simple compatibility between assets that are not necessarily identical.
For example, if we wanted to annotate genomic regions based on aligned BAM file using feature annotations we would need the reads to share just the coordinate system (the sequence lengths and names). Importantly, it does not require sequence identity at all. In this case, we would like to confirm that the given `bowtie2_index` asset is at least compatible with the coordinate system of feature annotation file.

To sum up, we need a more detailed comparison that may ignore sequences but allow to compare other aspects of the reference genome.

## Solution

Refgenie facilitates such fine-grained comparison of genomes via the `refgenie compare` command. When you initialize a genome with `refgenie genome init`, refgenie loads its sequences into the RefgetStore as a GA4GH [sequence collection](https://ga4gh.github.io/seqcol-spec/) (seqcol): the sequence names, lengths, and sequence digests. `refgenie compare` runs the standard seqcol comparison on two of these collections, so it needs no asset to be built first.

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

For the coordinate-system question in the motivation above, look at `name_length_pairs`: if every name-length pair of one genome is shared with the other, the two share a coordinate system for those sequences, even if the sequences themselves differ.

If one or both genomes come from a remote server, or the two genomes live in different RefgetStores, refgenie compares them from their sequence collections and the result may use a simpler layout (counts and the shared or unique elements for `names`, `lengths`, and `sequences`).

## Other ways to compare

- **Server**: a server-mode instance (`refgenie serve`) mounts the GA4GH seqcol service at `/seqcol`, which offers the standard comparison endpoint `GET /seqcol/comparison/{digest_a}/{digest_b}`.
- **AI assistants**: the [MCP server](mcp.md) provides the same comparison as the `compare_genomes` tool.
