# Genome aliases

Refgenie identifies each genome by a *digest*: a 32-character ID computed from its sequences, such as `pGauuPa6DYFFjNGTjVBbPsYgAmS5f3y2`. Digests are exact but hard to type, so refgenie also gives genomes *aliases*: readable names such as `hg38` or `GRCh38`. One genome can have several aliases, and you can use an alias anywhere refgenie asks for a genome.

## List aliases

```console
refgenie alias get
```

```console
┏━━━━━━━━━━━┳━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┓
┃ Alias     ┃ Genome Digest                    ┃
┡━━━━━━━━━━━╇━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━┩
│ hg38      │ pGauuPa6DYFFjNGTjVBbPsYgAmS5f3y2 │
│ GRCh38    │ pGauuPa6DYFFjNGTjVBbPsYgAmS5f3y2 │
│ mm10      │ hP1GNDAYRYvMmuEBp94uTC10v9W406oS │
└───────────┴──────────────────────────────────┘
```

Look up the digest of some aliases with `-a`, or the aliases of a digest with `-g` (use one or the other):

```console
refgenie alias get -a hg38,GRCh38
refgenie alias get -g pGauuPa6DYFFjNGTjVBbPsYgAmS5f3y2
```

## Set aliases

Give one or more aliases (comma-separated) and the genome's digest:

```console
refgenie alias set -a hg38 -d pGauuPa6DYFFjNGTjVBbPsYgAmS5f3y2
refgenie alias set -a hg38,GRCh38,human -d pGauuPa6DYFFjNGTjVBbPsYgAmS5f3y2
```

New aliases are added to the genome's existing ones. Add `-r`/`--reset` to remove the existing ones first, for example to rename a genome:

```console
refgenie alias set -a hg38_v2 -d pGauuPa6DYFFjNGTjVBbPsYgAmS5f3y2 --reset
```

The genome must already be in your database (from `refgenie genome init` or a pull). `-f`/`--force` sets the alias anyway.

## Remove aliases

```console
refgenie alias remove -a hg38
refgenie alias remove -a hg38,GRCh38
```

The genome and its assets stay, reachable by digest or by any alias it still has.

## Use aliases

```console
refgenie pull hg38/fasta
refgenie list -g hg38
refgenie getseq -g hg38 -l chrM:1-60
refgenie compare hg38 mm10
```
