---
name: refgenie-api
description: >-
  Query and download reference genome assets from a live refgenie instance
  over its HTTP API: list the genomes a server holds, resolve a name like
  hg38 to its GA4GH sequence-collection digest, find the asset carrying a
  FASTA, chromosome sizes, or a bwa/bowtie2/salmon index, and download the
  individual files or the whole archive. Use whenever you have a refgenie
  server or dashboard URL, or need reference genome files over HTTP without
  installing anything. To drive the refgenie command-line tool on a machine
  you control, read https://docs.refgenie.org/SKILL.md instead.
homepage: https://refgenie.org
---

# Refgenie HTTP API

You are reading this from a running **refgenie** instance. Refgenie manages
reference genome assets (FASTA files, chromosome sizes, aligner indexes,
annotations) and identifies every genome by a **GA4GH digest computed from its
sequences**, so you can check that two genomes are the same.

Everything here is a plain HTTP `GET`, so you need nothing installed to browse
or download from this server. This file covers one instance's HTTP API and MCP
endpoint. To install refgenie and run its command-line tool (`pull`, `build`,
`seek`, `getseq`), read <https://docs.refgenie.org/SKILL.md>.

---

## Step 0: discover this instance

Always start here. It gives the API base, the instance's mode, and which of the
endpoints below exist on it.

```bash
BASE=https://api.refgenie.org        # or http://localhost:8080 for a local dash
curl -s "$BASE/service-info" | jq '.refgenie'
```

The `refgenie` block:

| Field | Meaning |
|---|---|
| `mode` | `"server"` (public, read-only, has archives and downloads) or `"local"` (one user's `refgenie dash`, has the command surface) |
| `api_base` | `"/v4"`. Join it to the origin. |
| `root_path` | Sub-path prefix when mounted behind a proxy; usually `""` |
| `refgenie_version` | The running package version |
| `service_name` | `"refgenie server"` or `"refgenie local dashboard"` |
| `capabilities` | Booleans. **Check these before using a gated endpoint.** |
| `links` | `docs`, `github`, `openapi` |
| `web_ui` | Build stamp of the bundled web interface |

Capability keys, and what each gates:

| Key | True in | Gates |
|---|---|---|
| `downloads` | server | individual file downloads under an asset |
| `archives` | server | the archive listing and `.tgz` downloads |
| `seqcol` | server | the `/seqcol` GA4GH sequence-collections service |
| `drs` | server | the GA4GH DRS endpoints |
| `pull` `build` `delete` `genome_init` `aliases_write` `subscriptions` `data_channels` `remote_browse` `jobs` `jobs_cancel` `local_files` | local | the `/v1` command surface |
| `recipes_write` `asset_classes_write` | never | reserved; these stay CLI-only |

A missing key reads as `false`.

In server mode, `/service-info` also has a top-level `seqcol` block: `url`
(the seqcol mount, `/seqcol`), `service_info`, and `refget_store.url`, the
public RefgetStore that holds the sequences.

**If `/service-info` 404s or returns HTML**, this origin is a statically hosted
UI whose API lives elsewhere. Use `https://api.refgenie.org` as the base.

The machine-readable schema is at `GET $BASE/openapi.json`, with browsable
docs at `GET $BASE/docs`.

---

## Core concepts

- **Genome digest**: the GA4GH sequence-collection digest, computed from the
  sequences. It is the genome's identity: the same sequences give the same
  digest anywhere.
- **Alias**: a human name like `hg38` or `GRCh38` that points to a digest.
  Each instance keeps its own aliases. **An alias is not a digest.**
- **Asset group**: a named slot for one kind of asset on one genome, like
  `hg38`'s `bwa_index`. A group holds one or more assets and has a default.
- **Asset**: one built artifact in a group, with a name (`default` for
  `fasta`; the tool version for an index, like `bowtie2_index:2.5.5`)
  and its own asset digest. This is the thing you download.
- **Asset class**: the type of an asset, which defines its seek keys.
- **Seek key**: a named file or value inside an asset. The `fasta` class has
  seek keys `fasta` (the `.fa`), `fai` (the `.fa.fai`), and `chrom_sizes` (the
  `.chrom.sizes`). Read the seek key's `value` off the asset record instead of
  guessing a filename.
- **Registry path**: the CLI's way to name an asset,
  `genome/asset_group.seek_key:asset`, e.g. `hg38/fasta.fai:default`. Only
  `genome/asset_group` is required; the seek key defaults to the asset's
  default seek key and the asset to the group's default asset. Registry paths
  appear in CLI output, local-mode job results, and DRS object names. They are
  not API paths and not a field on the `/v4` asset record.
- **Serving mode**: a served asset is staged as `file` (one request per file)
  and/or `archive` (one `.tgz`). The asset record lists these in
  `serving_modes`.

---

## Endpoint reference

All paths are relative to the origin. `{braces}` are path parameters.

### Every instance, both modes

| Method + path | Returns |
|---|---|
| `GET /service-info` | discovery document (above) |
| `GET /ping` | localhost-bridge handshake (see below) |
| `GET /openapi.json` | full OpenAPI schema |
| `GET /docs` | Swagger UI |
| `GET /v4/genomes` | page of genomes, most assets first |
| `GET /v4/genomes/{digest}` | one genome, by **digest only** |
| `GET /v4/aliases` | page of `{name, genome_digest}` |
| `GET /v4/aliases/{name}` | **resolve an alias to a digest**, plus its FHR metadata and level-2 sequence collection |
| `GET /v4/asset_groups` | page of asset groups |
| `GET /v4/asset_groups/{id}` | one asset group, by integer id |
| `GET /v4/assets` | page of assets |
| `GET /v4/assets/{digest}` | one asset, with its seek keys |
| `GET /v4/assets/{asset_digest}/files` | `{asset_digest, files: [...]}`, the downloadable file list |
| `GET /v4/asset_classes` | page of asset classes |
| `GET /v4/asset_classes/{id}` | one asset class, by integer id |
| `GET /v4/recipes` | page of recipes |
| `GET /v4/recipes/{id}` | one recipe, by integer id |
| `GET /v4/relationships/{asset_digest}` | `{asset_digest, parents, children}`; add `?expand=true` for full objects |
| `GET /v4/staged_assets` | which assets are staged, and how |
| `GET /v4/staged_assets/{id}` | one staged-asset record |
| `GET /v4/configurations` | instance configuration |
| `GET /v4/configurations/{id}` | one configuration |

### Server mode only

| Method + path | Capability | Returns |
|---|---|---|
| `GET /v4/archives` | `archives` | page of archive records |
| `GET /v4/archives/{asset_digest}/download` | `archives` | the `.tgz`, **or a 307 to cloud storage** |
| `GET /v4/assets/{asset_digest}/files/{file_path}` | `downloads` | one file, **or a 307 to cloud storage** |
| `GET /v4/summary` | none | `{genomes, asset_groups, assets}` |
| `GET /v4/species/summary` | none | per-species `{genomes, asset_classes, assets}` |
| `GET /ga4gh/drs/v1/service-info` | `drs` | GA4GH DRS service info |
| `GET /ga4gh/drs/v1/objects/{object_id}` | `drs` | a DRS object |
| `GET /ga4gh/drs/v1/objects/{object_id}/access/{access_id}` | `drs` | a DRS access URL |
| `GET /ga4gh/drs/v1/objects/{object_id}/access/{access_id}/bytes` | `drs` | the bytes, for the `file` and `archive` access methods |
| `GET /data_channel/` | none | the data channels (recipe and asset-class sources) this server relays |

The DRS routes follow the GA4GH spec path `/ga4gh/drs/v1/...`, outside `/v4`.
The GA4GH sequence-collections service is at `/seqcol`, with its own
`/seqcol/service-info`. The seqcol service and the MCP endpoint are mounted
sub-applications, so they are missing from `/openapi.json`.

### Local mode only (`refgenie dash`)

A local dash has a command surface under `/v1`: pull, build, genome init,
alias and subscription edits, file viewing, and a job queue with a live SSE
stream. Read `GET /openapi.json` for the request and response shapes.

Two rules apply to all of it:

- **Every mutation needs the `X-Refgenie-Action` header** (any value).
  Without it you get `403 missing_action_header`. The header guards against
  cross-site requests.
- **A non-loopback `Host` header gets `421 forbidden_host`.** The dash
  answers only on loopback.

For changes, prefer the CLI: read <https://docs.refgenie.org/SKILL.md>.

---

## Query parameters

Every list endpoint takes the same paging and search parameters, in the
GA4GH page-based style.

| Param | Default | Notes |
|---|---|---|
| `page` | `0` | the page to return, counted from `0`; a page past the last one is a **400** |
| `page_size` | `100` | results per page; over `1000` is a **400**, under `1` a **422** |
| `q` | none | search term |
| `search_fields` | all searchable | comma-separated; an unknown field is a **422** |
| `operator` | `contains` | `eq`, `contains`, `starts_with`, `ends_with` |

Searchable fields, per endpoint:

| Endpoint | Fields |
|---|---|
| `/v4/genomes` | `digest`, `description`, `species` (every species name: scientific, common, and other common names like brown rat), `synonyms` (other names of the genome, like GRCh38), `assembly_source`, `assembly_accession`, `aliases` |
| `/v4/assets` | `name`, `digest`, `path` |
| `/v4/asset_groups` | `name` |
| `/v4/asset_classes` | `name`, `version`, `description` |
| `/v4/recipes` | `name`, `version`, `description` |
| `/v4/aliases` | `name` |
| `/v4/staged_assets` | `asset_digest`, `mode` |

Genome rows carry `species_name`, `common_name`, `taxon_id`, and a `species`
object (`taxon_id`, `scientific_name`, `common_name`, `synonyms`), which is
`null` when the species list has no row for that taxon.

Exact-match filters:

| Endpoint | Filters |
|---|---|
| `/v4/genomes` | `digest`, `alias`, `has_assets` (`true`/`false`), `taxon_id` (NCBI taxon id) |
| `/v4/assets` | `genome_digest`, `asset_group_name`, `asset_group_id`, `name`, `recipe_name` |
| `/v4/asset_groups` | `genome_digest`, `asset_class`, `asset_group_name`, `asset_group_id` |
| `/v4/asset_classes` | `name`, `version` |
| `/v4/recipes` | `name`, `version`, `output_asset_class` |
| `/v4/aliases` | `name`, `genome_digest` |
| `/v4/staged_assets` | `asset_digest`, `mode` |

**Every list returns this envelope:**

```json
{
  "results": [ ... ],
  "pagination": { "page": 0, "page_size": 100, "total": 1234 }
}
```

The count is `pagination.total`; there is no top-level `total` or `items` key.
The last page is `(total - 1) // page_size`; stop when
`(page + 1) * page_size >= total`. Page `0` of an empty list is an empty
`results`, not an error. Every list has a fixed order, so walking the pages
returns each row once. To count a collection cheaply, request `page_size=1`
and read `pagination.total`.

---

## Workflow 1: from a genome name to a file on disk

Do not skip step 2.

```bash
BASE=https://api.refgenie.org

# 1. Discover. Confirm capabilities.downloads / .archives before step 5.
curl -s "$BASE/service-info" | jq '.refgenie | {mode, api_base, capabilities}'

# 2. Resolve the human name to a digest. This is the ONLY correct way.
curl -s "$BASE/v4/aliases/hg38" | jq '{alias, digest}'

# 3. List that genome's assets.
curl -s "$BASE/v4/assets?genome_digest=<GENOME_DIGEST>&page_size=100" \
  | jq '.results[] | {digest, asset_class_name, asset_group_name, serving_modes, size}'

# 4. Read the asset's seek keys -- these name the files inside it.
curl -s "$BASE/v4/assets/<ASSET_DIGEST>" | jq '{asset_class_name, seek_keys}'

# 5a. Download individual files (needs capabilities.downloads).
curl -s "$BASE/v4/assets/<ASSET_DIGEST>/files" | jq '.files'
curl -L -O "$BASE/v4/assets/<ASSET_DIGEST>/files/<FILE_PATH_FROM_THAT_LIST>"

# 5b. Or take the whole asset as one tarball (needs capabilities.archives).
curl -L -O "$BASE/v4/archives/<ASSET_DIGEST>/download"
```

**`-L` is required.** Both download routes answer with a `307` to cloud
storage when the instance has a remote configured. Without `-L` you get an
empty body.

Take the file path in step 5a from the `files` list. A path that is not in the
list is a `404`, because the list is the path-traversal guard.

Step 3 reports each asset's `serving_modes`: `["file"]`, `["archive"]`, or
both. Use 5a only for an asset staged as `file`, and 5b only for one staged as
`archive`.

## Workflow 2: find out what a server holds

```bash
# Count everything, cheaply.
curl -s "$BASE/v4/genomes?page_size=1" | jq '.pagination.total'

# Only the genomes that have assets (most have none), or one species by taxon.
curl -s "$BASE/v4/genomes?has_assets=true&page_size=20" \
  | jq '.results[] | {digest, aliases, species_name, common_name, asset_count}'
curl -s "$BASE/v4/genomes?taxon_id=10090&page_size=20" | jq '.results[] | {digest, aliases}'

# Search across species, aliases, descriptions and accessions.
curl -s "$BASE/v4/genomes?q=mouse&page_size=20" \
  | jq '.results[] | {digest, aliases, species_name, common_name, asset_count}'

# What kinds of asset exist on this instance?
curl -s "$BASE/v4/asset_classes?page_size=100" | jq '.results[] | {name, description}'

# Which genomes have a bowtie2 index?
curl -s "$BASE/v4/assets?asset_group_name=bowtie2_index&page_size=100" \
  | jq '.results[] | {genome_digest, digest, size}'
```

On a server-mode instance, `GET /v4/summary` gives
`{genomes, asset_groups, assets}` in one request, and `GET /v4/species/summary`
breaks the counts down per species.

## Workflow 3: verify genome identity

Two genomes are the same genome if their sequence-collection digests match,
whatever they are called. `GET /v4/aliases/{name}` returns the level-2
sequence collection (`names`, `lengths`, `sequences`) in its `collection`
field, which is often enough to compare. For full comparisons on a
server-mode instance, use the GA4GH sequence-collections service at `/seqcol`
(start from `GET /seqcol/service-info`).

---

## The MCP server

Refgenie has a **read-only** Model Context Protocol server with thirteen
tools. It reads through the same code as the HTTP API, so both give the same
answers.

**Over HTTP, against a server-mode instance** (Streamable HTTP transport,
nothing to install):

```
https://api.refgenie.org/mcp
```

To add the public server to Claude Code:

```bash
claude mcp add --transport http refgenie-public https://api.refgenie.org/mcp
```

Use `/mcp` exactly, with no trailing slash. Only server mode has `/mcp`;
check `"mode": "server"` in `/service-info` first.

**Over stdio, against the user's own database**, with the `refgenie-mcp`
command from the refgenie package:

```bash
pip install "refgenie[mcp]>=1.0.0a1"
claude mcp add refgenie refgenie-mcp
```

`refgenie-mcp` takes no arguments. It reads the user's config from
`REFGENIE_HOME_PATH` (default `~/.refgenie`).

The tools:

| Tool | Purpose |
|---|---|
| `catalog_summary` | counts of genomes, asset groups and assets, and the species that have assets |
| `list_genomes` | genomes, most assets first; filters `has_assets`, `taxon_id` |
| `search_genomes` | substring search over digest, alias, any species name (scientific, common, other common names like brown rat), other genome names (GRCh38), description, accession; same filters |
| `get_genome` | one genome by `alias` or `digest`, with its FHR metadata and its species' names |
| `list_assets` | assets, optionally of one genome (`alias` or `digest`) and/or `asset_class` |
| `get_asset` | one asset by digest: names, seek keys, parents, children, and how to get it |
| `list_asset_classes` | asset classes and their seek key names |
| `list_recipes` | recipes, optionally by `output_asset_class` |
| `get_recipe` | one recipe in full: inputs, command templates, docker image, the build command |
| `lookup_digest` | what a digest is: a genome, an asset, or nothing |
| `get_genome_metadata` | seqcol metadata: sequence count, total length, source |
| `get_genome_sequences` | per-sequence names, lengths, digests |
| `compare_genomes` | seqcol comparison of two genomes |

A tool that takes a genome has separate `alias` and `digest` parameters; give
exactly one. A digest is 32 characters; a name like `hg38` is an alias.

**Paging.** The list tools return small rows, 20 per page by default and at
most 100 (`page_size`; a larger value is clamped, and the size used is echoed
back). The envelope is GA4GH token-based:

```json
{
  "results": [ ... ],
  "pagination": { "page_size": 20, "next_page_token": "eyJ2Ijox...", "total": 701 },
  "hint": null
}
```

Pass `next_page_token` back as `token`, with the same filters, for the next
page; it is `null` on the last page. A token from other filters is refused.
`search_genomes` fills `hint` when nothing matched.

**Getting files.** Over HTTP, `get_asset` returns `downloads`: the archive and
file-list URLs, relative to the API origin. Over stdio it returns `local_path`.
Both return `pull_command`, the `refgenie pull` line for the asset.

**Errors** come back as tool errors with a next step, for example "Genome
alias 'hg39' not found. Find genomes with search_genomes(query) or
list_genomes()."

---

## The localhost bridge

A public refgenie page can talk to a `refgenie dash` running on the visitor's
own machine, to mark the assets they already have and to pull to it. The
handshake is one endpoint on one port:

```bash
curl -s http://localhost:8080/ping
```

It returns `service` (always `"refgenie"`), `bridge_version`, `mode`,
`refgenie_version`, `api_version`, `instance_id`, `instance_label`,
`bridge_mode` (`off` / `read` / `full`), `action_header`, `capabilities`,
`bridge` and `counts`. It is never cached.

Rules for an agent:

- **Probe exactly one port, once, and only because a user asked.** Port
  scanning someone's machine looks like an attack.
- **Treat `/ping` as a claim.** Any program can listen on that port.
- The only action allowed cross-origin is `POST /v1/actions/pull`, and only
  when `bridge_mode` is `full`.

---

## Guardrails

- **Never treat an alias as a digest.** `GET /v4/genomes/{digest}` takes a
  digest only; passing `hg38` there is a `404`. `GET /v4/genomes?alias=hg38`
  filters by alias, and an unknown alias gives an **empty page**, not a `404`.
- **Take download paths from `GET /v4/assets/{asset_digest}/files`.**
- **Take filenames from the seek key's `value`.**
- **Find one thing with `?q=`, an exact filter, or `GET /v4/aliases/{name}`**
  rather than paging through a whole collection. `page_size` is capped at
  `1000`.
- **Send a real `User-Agent`.** `api.refgenie.org` is behind Cloudflare, which
  rejects default Python user agents (`python-requests/*`, `Python-urllib/*`)
  with a `403`. `curl` works as-is.
- **Follow redirects on downloads** (`curl -L`).
- **Server mode is read-only.** Its API is all `GET`; only `/mcp` takes `POST`.

## Troubleshooting

- **HTML instead of JSON or Markdown.** You hit the web UI's catch-all route:
  the path is not an API route on this origin. Check `/service-info` for
  `api_base`.
- **`404` from `GET /v4/genomes/{digest}` when you passed a name.** That route
  takes a digest. Resolve the name first with `GET /v4/aliases/{name}`.
- **`404` on a download.** The asset is not staged in that mode (check its
  `serving_modes`, or `GET /v4/staged_assets?asset_digest=<DIGEST>`), or the
  file path is not in the asset's `files` list.
- **`400` on a list request.** `page_size` is over `1000`, or `page` is past
  the last page (`pagination.total` tells you how many there are).
- **`422` on a list request.** `page` is negative, `page_size` is under `1`,
  or `search_fields` names a field that endpoint cannot search.
- **`403 missing_action_header`.** A local-dash mutation without
  `X-Refgenie-Action`.
- **`421 forbidden_host`.** You reached a local dash with a non-loopback
  `Host` header.
- **Empty body on a download.** You did not follow the `307`.

## Where to read more

- The refgenie CLI and Python API: <https://docs.refgenie.org/SKILL.md>
- Human documentation: <https://docs.refgenie.org>
- This instance's OpenAPI schema: `GET /openapi.json`
- Source: <https://github.com/refgenie/refgenie>
