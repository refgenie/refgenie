# Migrating from legacy refgenie (0.x) to refgenie1 (1.0)

Refgenie 1.0 is a complete rewrite. The CLI, configuration, asset
storage, and seek-key contract all changed. This page documents the
breaking changes and how to migrate downstream code.

## Environment variables

See the [README "Environment variables"](../README.md#environment-variables) section for the
full mapping. The main rename:

- `$REFGENIE` → `$REFGENIE_DB_CONFIG_PATH`

If `$REFGENIE` is set in your shell but `$REFGENIE_DB_CONFIG_PATH` is not,
refgenie1 prints a one-shot warning to stderr and ignores `$REFGENIE`.

## Seek keys: the `.dir` convention is gone

### What changed

Legacy refgenconf gave every asset an implicit `dir` seek key that
returned the asset's containing folder. Refgenie1 removed this.
Asset classes now declare what they expose via typed seek keys
(`file`, `directory`, `prefix`, `string`, `json`).

### Why

`.dir` was a workaround for legacy refgenconf having no type system.
For aligner indexes — the most common use of `.dir` — consumers had
to ask for the directory and then hand-construct the prefix
(`{dir}/{genome}`). Refgenie1 has a `prefix`-typed seek key that
returns the prefix directly: no concatenation, no string hacks.

### How to migrate

For every consumer that used `<asset>.dir`, switch to the asset's
declared seek key. In most cases each asset class declares a single
canonical seek key with the same name as the asset class, and the
default-seek-key resolver picks it up automatically — so you can drop
the `.<seek_key>` suffix entirely.

Examples (registry-path syntax):

| Legacy | Refgenie1 |
|---|---|
| `hg38/bowtie2_index.dir` (returns directory; caller appends `/{genome}`) | `hg38/bowtie2_index` (returns prefix directly) |
| `hg38/bwa_index.dir` (caller appends `/{genome}.fa`) | `hg38/bwa_index` (returns prefix directly) |
| `hg38/fasta.dir` | `hg38/fasta` (returns the FASTA file path; default seek key) |

For pipeline-interface populator templates (e.g. PEPATAC):

```diff
- --genome-index { refgenie[sample.genome].bowtie2_index.dir }
+ --genome-index { refgenie[sample.genome].bowtie2_index }
```

### Tools that genuinely need a directory

If a tool requires a directory path (e.g., STAR's `--genomeDir`), the
asset class should declare a `directory`-typed seek key explicitly.
This is intentional — directories are first-class asset interfaces,
not a universal escape hatch.

## `seek` returns human-readable alias paths

Asset content is now addressed by digest and stored under
`data/<genome_digest>/<asset_group>/<content_digest>/`, with human-readable
names living in an `assetname` table and rendered as a tree of per-file
symlinks under `alias/<alias>/<asset_group>/<name>/`. One content digest may
carry more than one name (for example, a toolchain bump that produces
byte-identical output).

Because the content path now contains a digest rather than a name, `seek`
returns the **alias path** by default — the human-readable view named for the
alias you asked about, with filenames rewritten to that alias:

```console
$ refgenie seek hg38/bwa_index:0.7.19
<genome_folder>/alias/hg38/bwa_index/0.7.19/hg38.fa

$ refgenie seek --abs hg38/bwa_index:0.7.19
<genome_folder>/data/<genome_digest>/bwa_index/<content_digest>/<genome_digest>.fa
```

Pass `--abs` (library: `seek(..., abs_path=True)`) when you need the stable,
digest-addressed content path instead of the alias view. Both point at the
same bytes.

## Aligner-index recipes: prefix typing

The `bowtie2_index`, `bwa_index`, and similar aligner-index asset
classes now declare their canonical seek key as `type: prefix`. The
seek-key value resolves to `<asset_folder>/<genome>` — exactly what
`bowtie2 -x`, `bwa mem`, etc. expect on their command line.

If you maintain a custom recipe for an aligner-index-style asset:

- Declare `type: prefix` (not `file` pointing at one of the index
  files).
- Set `value: "{genome}"` so the resolved path is
  `<asset_folder>/<genome>`.
- In the recipe's command template, use the aligner's prefix flag
  (e.g., `bwa index -p <prefix>`) so the output files are written
  with the prefix you declared.

## Removed APIs and CLI subcommands

- `refgenie pull <genome>/<asset>` — replaced by `refgenie genome init`
  + `refgenie add` against a subscribed source.
- `refgenie asset list` — use `refgenie list -g <genome>`.
- `Refgenie.list_seek_keys_values()` — restored on `AssetManager` in
  refgenie1; signature returns
  `{genome: {asset_group: {asset: {seek_key: str}}}}`.
- `Refgenie.asset.seek(...)` — now returns `str` (was `pathlib.Path`).
  JSON-serialization and Jinja-templating callers no longer need to
  coerce.

## Build bookkeeping moved out of the asset directory

### What changed

Build bookkeeping — the pypiper log, commands file, profile, cleanup
script, `stats.yaml`, status flags, and the build-completion flag —
used to be written to a `_refgenie_build/` directory *inside* the asset
output folder. It now lives in a separate top-level tree:

```
<genome_folder>/builds/<genome_alias>/<asset_group>/<asset_name>/
```

The tree is keyed by build invocation, not by asset digest, so it does
not follow asset content if that content moves or is deduplicated.

### Why

Co-locating bookkeeping with content meant every operation over an asset
directory had to exclude it by name — the content digest, tarballs, rsync,
the served-file listing, and prefix seek-key matching each carried their
own exclusion, and file-mode staging was missing one, publishing build
logs to served buckets. Separating the trees removes all five exclusions
and the leak, and it removes the ~100 lines that reconciled the
completion flag between the alias path snakemake was told about and the
digest path the flag was actually written to. Those are now one path.

### Impact on asset digests

`asset.digest` is a hash over the asset directory's files. Previously
`_refgenie_build/` was excluded from that hash by an explicit rule;
now it is excluded structurally, because it is not there. For a genome
folder with no stale `_refgenie_build/` directories the resulting digest
is unchanged — but **digests computed before this change are not
guaranteed to compare against digests computed after it**, and no
migration reads the old location. Rebuild rather than reconcile.

### Impact on Snakemake consumers

`get_asset_build_target_template()` now returns a path under `builds/`
rather than under `alias/`. Consumers that call it need no change.
Consumers that hardcoded the old `_refgenie_build/` flag path must
switch to calling `get_asset_build_target_template()`. The flag basename
(`{genome}_{group}__{asset}.flag`) is unchanged.

## Config object

Legacy refgenconf's `RGC` is gone. Use the `Refgenie` class directly.
Constructor accepts `database_config_path: str | Path | None`.
