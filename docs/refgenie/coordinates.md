# 1-based or 0-based coordinates?

Refgenie follows one rule everywhere it reads a region of a genome:

- **A region written as text is 1-based and inclusive**, like samtools, UCSC and IGV. `chr17:43,044,295-43,044,394` means base 43,044,295 through base 43,044,394, which is 100 bases. Commas are allowed. This is how `refgenie getseq -l` and the Extract tab on [refgenie.org](https://refgenie.org) read a region.
- **Separate start and end numbers are 0-based and half-open**, like the GA4GH refget API. `start=43044294, end=43044394` is the same 100 bases: the first base is number 0, and the end base is not included. This is how the MCP tool `get_sequence` and the Python API read a region. BED files are 0-based and half-open too.

To turn a text region into start and end numbers, subtract 1 from the start and keep the end.

| You have | Example | Read as |
|---|---|---|
| A region as text | `chr17:43,044,295-43,044,394` | 1-based, inclusive |
| Start and end numbers | `start=43044294`, `end=43044394` | 0-based, half-open |
| A BED line | `chr17 43044294 43044394` | 0-based, half-open |

## The same 100 bases, four ways

These all return the first 100 bases of BRCA1 in `hg38`, starting `TGGAAGTGTTTGCTACCAAG`.

On the command line:

```bash
refgenie getseq -g hg38 -l chr17:43,044,295-43,044,394
```

If your text region is already 0-based, for example copied from a BED file, add `--zero-based`:

```bash
refgenie getseq -g hg38 -l chr17:43044294-43044394 --zero-based
```

Either way, `getseq` prints how it read the region on stderr, so you can check it: `chr17:43,044,295-43,044,394 (1-based, inclusive): 100 bases`.

From an AI assistant, over [MCP](mcp.md):

```
get_sequence(alias="hg38", name="chr17", start=43044294, end=43044394)
```

From Python, once the genome is registered (any `refgenie getseq` or `refgenie pull` does that):

```python
from refgenie import Refgenie

rgc = Refgenie()
rgc.sequence.get(rgc.alias.resolve("hg38"), "chr17", 43044294, 43044394)
```
