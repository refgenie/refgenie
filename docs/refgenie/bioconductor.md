# Use refgenie genomes in R and Bioconductor

Many [Bioconductor](https://bioconductor.org) tools take a reference genome as a [BSgenome](https://bioconductor.org/packages/BSgenome/) object. A BSgenome object comes from an R data package, such as `BSgenome.Hsapiens.UCSC.hg38`, that holds a genome's sequences and reads any region with `getSeq()`. The R package [BiocRefgetStore](https://github.com/databio/BiocRefgetStore) builds such a package from any refgenie genome, so your R code uses exactly the genome your other tools use.

BiocRefgetStore reads a genome from the [RefgetStore](https://docs.refgenie.org/refget/refgetstore-explained/) where refgenie keeps its sequences. It can also read sequences straight from the store with `getSeq()`; see the [BiocRefgetStore documentation](https://docs.refgenie.org/biocrefgetstore/getting-started/).

## Install

BiocRefgetStore depends on [gtars](https://github.com/databio/gtars), which comes as a prebuilt R package for Linux (x86_64) and macOS (Apple Silicon) on the [gtars releases page](https://github.com/databio/gtars/releases). Install gtars, then BiocRefgetStore and [BSgenomeForge](https://bioconductor.org/packages/BSgenomeForge/):

```r
install.packages(c("data.table", "BiocManager", "remotes"))
BiocManager::install(c("BiocGenerics", "IRanges", "BSgenomeForge"))

# Linux x86_64:
install.packages(
  "https://github.com/databio/gtars/releases/download/gtars-r-v0.10.1/gtars_0.10.1_R_x86_64-pc-linux-gnu.tar.gz",
  repos = NULL
)
# macOS Apple Silicon (arm64):
install.packages(
  "https://github.com/databio/gtars/releases/download/gtars-r-v0.10.1/gtars_0.10.1_R_arm64-apple-darwin.tgz",
  repos = NULL
)

remotes::install_github("databio/BiocRefgetStore")
```

The prebuilt gtars packages are built for R 4.4. For another R version or platform, clone the [gtars repository](https://github.com/databio/gtars) and build its R package from source, which needs a [Rust toolchain](https://www.rust-lang.org/tools/install):

```r
# git clone https://github.com/databio/gtars.git
install.packages("path/to/gtars/gtars-r", repos = NULL, type = "source")
```

## Make a BSgenome package from a refgenie genome

Three calls turn a genome on the public refgenie server into an installed BSgenome package. This example uses the yeast genome `R64-1-1`:

```r
library(BiocRefgetStore)

genome <- RefgetGenome.from_refgenie("R64-1-1")
pkg <- forgeBSgenome(genome, "BSgenome.Scerevisiae.refgenie.R6411",
                     organism = "Saccharomyces cerevisiae", circ_seqs = "Mito",
                     dest = tempdir())
install.packages(pkg, repos = NULL, type = "source")

library(BSgenome.Scerevisiae.refgenie.R6411)
getSeq(Scerevisiae, "Mito", 1, 60)
```

`RefgetGenome.from_refgenie()` asks the refgenie server (by default `https://api.refgenie.org`) for the genome and for the public RefgetStore that holds it, and downloads the genome's metadata. `forgeBSgenome()` then downloads each sequence once and writes the package. For yeast this downloads about 3 MB and takes about 20 seconds, most of it package installation. For human `hg38` the download is about 1.4 GB.

- **Genome names.** Use any name listed at [refgenie.org](https://refgenie.org/genomes), or by `refgenie genome browse --server-url https://api.refgenie.org`.
- **Package names.** A BSgenome package is named `BSgenome.<organism>.<provider>.<build>`, and each part holds only letters and numbers, so `R64-1-1` becomes `R6411`. The package's BSgenome object is named after the second part, here `Scerevisiae`.
- **Circular sequences.** List them in `circ_seqs`, using the names from `seqnames(genome)`, such as `Mito` for this yeast genome or `chrM` for `hg38`.

## Use your local genome folder

A genome you initialized from a FASTA file with `refgenie genome init --fasta` has its sequences in your genome folder, so BiocRefgetStore reads it from disk. Refgenie keeps the store in `<genome folder>/.refget_store`, with genome names under the alias namespace `"refgenie"`:

```r
genome <- RefgetGenome.from_directory("~/.refgenie/genomes/.refget_store",
                                      namespace = "refgenie", alias = "mygenome")
```

The genome folder is `~/.refgenie/genomes` unless you set `REFGENIE_GENOME_FOLDER` (see [Configure](configuration.md)). Then call `forgeBSgenome()` as above.

A genome registered from a server with `refgenie genome init --server ... --no-build` has only its metadata in your genome folder. Open it with `RefgetGenome.from_refgenie()` and the same name.

## Trace a package back to its genome

Refgenie identifies each genome by its [sequence collection](https://ga4gh.github.io/refget/seqcols/) digest, computed from its sequences. `forgeBSgenome()` records that digest as the package's `genome` field:

```r
GenomeInfoDb::genome(Scerevisiae)[[1]]
#> [1] "ia9g4Myony-xxLsSc7eQOC_NdClGgPGu"
```

This is the digest that `refgenie genome list` shows, and the server describes the genome at `https://api.refgenie.org/v4/genomes/<digest>` (for this yeast genome, [refgenie.org/genomes/ia9g4Myony-xxLsSc7eQOC_NdClGgPGu](https://refgenie.org/genomes/ia9g4Myony-xxLsSc7eQOC_NdClGgPGu)). To check a package against its source, recompute the digest from its sequences:

```r
collection_digest(RefgetGenome.from_bsgenome(Scerevisiae))
#> [1] "ia9g4Myony-xxLsSc7eQOC_NdClGgPGu"
```

For more detail, including sizes and times, see the [Forge a BSgenome package from refgenie](https://docs.refgenie.org/biocrefgetstore/refgenie-bsgenome/) guide.
