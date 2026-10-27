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

# Refgenie instance — agent capability doc

You are reading this from a running **refgenie** instance. Refgenie manages
reference genome assets — FASTA files, chromosome sizes, aligner indexes,
annotations — and identifies every genome by a **sequence-derived GA4GH
digest** rather than by a name like "hg38", so genome identity is verifiable
instead of conventional.

Everything in this document is a plain HTTP `GET`. **You do not need to
install anything to browse or download from this server.**

**Scope.** This file describes the HTTP API of *one instance*. For installing
the package and driving the `refgenie` command-line tool on a machine you
control — `pull`, `build`, `seek`, `getseq` — read
<https://docs.refgenie.org/SKILL.md> instead. The two are complementary.

---

## Step 0: discover this instance

Always start here. It tells you the API base, which mode the instance is in,
and which of the endpoints below actually exist on it.

```bash
BASE=https://api.refgenie.org        # or http://localhost:8080 for a local dash
curl -s "$BASE/service-info" | jq '.refgenie'
```

The `refgenie` block:

| Field | Meaning |
|---|---|
| `mode` | `"server"` (public, read-only, has archives and downloads) or `"local"` (one user's `refgenie dash`, has the command surface) |
| `api_base` | Always `"/v4"` today. Join it to the origin. |
| `root_path` | Sub-path prefix when mounted behind a proxy; usually `""` |
| `refgenie_version` | The running package version |
| `service_name` | `"refgenie server"` or `"refgenie local dashboard"` |
| `capabilities` | 15 booleans. **Check these before using a gated endpoint.** |
| `links` | `docs`, `github`, `openapi` |
| `web_ui` | Build stamp of the bundled web interface |

Capability keys, and what each gates:

| Key | True in | Gates |
|---|---|---|
| `downloads` | server | individual file downloads under an asset |
| `archives` | server | the archive listing and `.tgz` downloads |
| `seqcol` | server | the `/seqcol` GA4GH sequence-collections service |
| `drs` | server | the GA4GH DRS endpoints |
| `pull` `build` `delete` `genome_init` `aliases_write` `subscriptions` `remote_browse` `jobs` `jobs_cancel` | local | the `/v1` command surface |
| `recipes_write` `asset_classes_write` | never | reserved; these stay CLI-only |

A key that is missing reads as `false`.

**If `/service-info` 404s or returns HTML**, this origin is a statically hosted
UI whose API lives elsewhere. Use the public API base
`https://api.refgenie.org` and continue.

The machine-readable schema for everything below is at
`GET $BASE/openapi.json`, with human-browsable docs at `GET $BASE/docs`.

---

## Core concepts

- **Genome digest** — the GA4GH sequence-collection digest, derived from the
  sequences themselves. This is the genome's real identity. Two genomes built
  from the same sequences get the same digest anywhere in the world.
- **Alias** — a human name like `hg38` or `GRCh38` that maps to a digest.
  Aliases are per-instance and are *not* identity. **Never treat an alias as a
  digest.**
- **Asset group** — the named bundle of assets of one class for one genome,
  e.g. `hg38`'s `bwa_index`.
- **Asset** — one concrete built artifact inside a group, identified by its own
  asset digest. This is the thing you download.
- **Asset class** — the typed shape an asset takes: which files it contains and
  under which seek keys.
- **Seek key** — a named file within an asset. The `fasta` asset class has seek
  keys `fasta` (the `.fa`), `fai` (the `.fa.fai`), and `chrom_sizes` (the
  `.chrom.sizes`). Read the seek key's `value` off the asset record; never
  guess a filename.
- **Registry path** — the CLI's `genome/asset_class:tag` shorthand, e.g.
  `hg38/bwa_index`. It shows up in CLI output, in local-mode job results and in
  DRS object names, but it is **not** an API path and it is **not** a field on
  the v4 asset record.
- **Staging mode** — a served asset is staged as `file` (individual files, one
  request each) and/or `archive` (one `.tgz`). An asset may have one, both, or
  neither. The asset record reports this as `serving_modes`.

---

## Endpoint reference

All paths are relative to the origin. `{braces}` are path parameters.

### Present on every instance, both modes

| Method + path | Returns |
|---|---|
| `GET /service-info` | discovery document (above) |
| `GET /ping` | localhost-bridge handshake (see below) |
| `GET /openapi.json` | full OpenAPI schema |
| `GET /docs` | Swagger UI |
| `GET /v4/genomes` | page of genomes |
| `GET /v4/genomes/{digest}` | one genome, by **digest only** |
| `GET /v4/aliases` | page of `{name, genome_digest}` |
| `GET /v4/aliases/{name}` | **resolve an alias to a digest**, plus its level-2 sequence collection |
| `GET /v4/asset_groups` | page of asset groups |
| `GET /v4/asset_groups/{id}` | one asset group, by integer id |
| `GET /v4/assets` | page of assets |
| `GET /v4/assets/{digest}` | one asset, with its seek keys |
| `GET /v4/assets/{asset_digest}/files` | `{asset_digest, files: [...]}` — the downloadable file list |
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
| `GET /v4/summary` | — | `{genomes, asset_groups, assets}` |
| `GET /v4/species/summary` | — | per-species `{genomes, asset_classes, assets}` |
| `GET /v4/ga4gh/drs/service-info` | `drs` | GA4GH DRS service info |
| `GET /v4/ga4gh/drs/objects/{object_id}` | `drs` | a DRS object |
| `GET /v4/ga4gh/drs/objects/{object_id}/access/{access_id}` | `drs` | a DRS access URL |

The DRS routes are also mounted unprefixed at `/ga4gh/drs/...`. The GA4GH
sequence-collections service is mounted at `/seqcol` and carries its own
`/seqcol/service-info`. Neither the seqcol service nor the MCP endpoint appears
in `/openapi.json`, because both are mounted sub-applications.

### Local mode only (`refgenie dash`)

A local dash exposes a command surface under `/v1` — pull, build, genome init,
alias and subscription edits, and a job queue with a live SSE stream. Read
`GET /openapi.json` for the shapes; they are not enumerated here because they
change with the dash.

Two rules apply to all of it:

- **Every mutation requires the `X-Refgenie-Action` header** (any value).
  Without it you get `403 missing_action_header`. This is the CSRF defence and
  it is not optional.
- **A non-loopback `Host` header gets `421 forbidden_host`.** The dash is
  loopback-only by design.

Prefer the CLI for mutations. Read <https://docs.refgenie.org/SKILL.md>.

---

## Query parameters

Every list endpoint takes the same pagination and search parameters.

| Param | Default | Notes |
|---|---|---|
| `offset` | `0` | must be `>= 0` |
| `limit` | `100` | `1`–`1000`; outside that range is a **422** |
| `q` | — | search term |
| `search_fields` | all searchable | comma-separated; an unknown field is a **422** |
| `operator` | `contains` | `eq`, `contains`, `starts_with`, `ends_with` |

Searchable fields, per endpoint:

| Endpoint | Fields |
|---|---|
| `/v4/genomes` | `digest`, `description`, `species_name`, `common_name`, `assembly_source`, `assembly_accession`, `aliases` |
| `/v4/assets` | `name`, `digest`, `path` |
| `/v4/asset_groups` | `name` |
| `/v4/asset_classes` | `name`, `version`, `description` |
| `/v4/recipes` | `name`, `version`, `description` |
| `/v4/aliases` | `name` |
| `/v4/staged_assets` | `asset_digest`, `mode` |

Additional exact-match filters:

| Endpoint | Filters |
|---|---|
| `/v4/genomes` | `digest`, `alias` |
| `/v4/assets` | `genome_digest`, `asset_group_name`, `asset_group_id`, `name`, `recipe_name` |
| `/v4/asset_groups` | `genome_digest`, `asset_class`, `asset_group_name`, `asset_group_id` |
| `/v4/asset_classes` | `name`, `version` |
| `/v4/recipes` | `name`, `version`, `output_asset_class` |
| `/v4/aliases` | `name`, `genome_digest` |
| `/v4/staged_assets` | `asset_digest`, `mode` |

**The paginated envelope is always:**

```json
{
  "items": [ ... ],
  "pagination": { "offset": 0, "limit": 100, "total": 1234 }
}
```

There is **no top-level `total`**, **no `results` key**, and **no `page`
parameter**. To count a collection without fetching it, request `limit=1` and
read `pagination.total`.

---

## Workflow 1: from a genome name to a file on disk

The canonical flow. Do not skip step 2.

```bash
BASE=https://api.refgenie.org

# 1. Discover. Confirm capabilities.downloads / .archives before step 5.
curl -s "$BASE/service-info" | jq '.refgenie | {mode, api_base, capabilities}'

# 2. Resolve the human name to a digest. This is the ONLY correct way.
curl -s "$BASE/v4/aliases/hg38" | jq '{alias, digest}'

# 3. List that genome's assets.
curl -s "$BASE/v4/assets?genome_digest=<GENOME_DIGEST>&limit=100" \
  | jq '.items[] | {digest, asset_class_name, asset_group_name, serving_modes, size}'

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
empty body and a redirect you did not follow.

The file path in step 5a must come from the `files` list in the response.
**Never construct it by convention** — a path not in that list is a `404` by
design, because the list is the path-traversal guard.

Step 3 reports each asset's `serving_modes`: `["file"]`, `["archive"]`, or
both. Take 5a only for an asset staged as `file`, and 5b only for one staged as
`archive`.

## Workflow 2: find out what a server holds

```bash
# Count everything, cheaply.
curl -s "$BASE/v4/genomes?limit=1" | jq '.pagination.total'

# Search across species, aliases, descriptions and accessions.
curl -s "$BASE/v4/genomes?q=mouse&limit=20" \
  | jq '.items[] | {digest, aliases, species_name, asset_count}'

# What kinds of asset exist at all on this instance?
curl -s "$BASE/v4/asset_classes?limit=100" | jq '.items[] | {name, description}'

# Which genomes have a bowtie2 index?
curl -s "$BASE/v4/assets?asset_group_name=bowtie2_index&limit=100" \
  | jq '.items[] | {genome_digest, digest, size}'
```

On a server-mode instance, `GET /v4/summary` gives
`{genomes, asset_groups, assets}` in one request, and `GET /v4/species/summary`
breaks the same counts down per species.

## Workflow 3: verify genome identity

Two genomes are the same genome if their sequence-collection digests match,
whatever they are called. On a server-mode instance the GA4GH
sequence-collections service is mounted at `/seqcol`; start from
`GET /seqcol/service-info`. `GET /v4/aliases/{name}` also returns the full
level-2 sequence collection in its `collection` field, which is often enough
to compare names and lengths without a second service.

---

## The MCP server

Refgenie ships a **read-only** Model Context Protocol server with eleven
tools. It never modifies anything.

**Over HTTP, against a server-mode instance** — Streamable HTTP transport, and
nothing to install:

```
https://api.refgenie.org/mcp/mcp
```

Note the doubled segment: the MCP sub-application is mounted at `/mcp` and
registers its own `/mcp` route inside itself, so `/mcp/mcp` is the working URL.
`POST /mcp` answers `405`. **A local `refgenie dash` does not mount MCP at
all** — check `"mode": "server"` in `/service-info` first.

**Over stdio, against the user's own database** — this needs the refgenie
Python package and the `refgenie-mcp` command it installs:

```bash
claude mcp add refgenie refgenie-mcp
```

`refgenie-mcp` takes no arguments. It reads the user's config from
`REFGENIE_HOME_PATH` (default `~/.refgenie`).

**About installing refgenie 1.x, honestly:** this generation of refgenie is not
published on any package index yet. PyPI's `refgenie` project stops at the
legacy `0.13.0` line, which is a different codebase with no MCP server, and
there is no `refgenie1` project on PyPI or Test PyPI. So there is currently no
`pip install` command that gets you the code this instance is running. Until a
1.x release is published, use the HTTP API in this document, or the HTTP MCP
endpoint above — neither needs anything installed.

The tools:

| Tool | Purpose |
|---|---|
| `list_genomes` | all genomes with aliases, species, description |
| `search_genomes` | substring search over species, alias, description, accession |
| `get_genome` | one genome by alias or digest, with its asset groups |
| `list_asset_classes` | registered asset classes, seek keys, serving modes |
| `list_recipes` | registered recipes, inputs, command templates |
| `list_assets` | assets, optionally filtered by genome and/or asset class |
| `get_asset` | asset detail by digest: seek keys, parents, children |
| `lookup_digest` | universal digest lookup — genome first, then asset |
| `get_genome_metadata` | seqcol metadata: sequence count, total length, source |
| `get_genome_sequences` | per-sequence names, lengths, digests |
| `compare_genomes` | seqcol comparison of two genomes |

Every tool that takes a genome accepts an alias or a digest and resolves it
for you. Every tool returns a JSON-encoded string.

---

## The localhost bridge

A public refgenie page can talk to a `refgenie dash` running on the visitor's
own machine, to badge the assets they already have and to pull to it. The
handshake is one endpoint on one port:

```bash
curl -s http://localhost:8080/ping
```

It returns `bridge_version`, `mode`, `refgenie_version`, `api_version`,
`instance_id`, `instance_label`, `bridge_mode` (`off` / `read` / `full`),
`action_header`, `capabilities`, `bridge` and `counts`. It is deliberately
separate from `/service-info` and is never cached.

Rules that matter to an agent:

- **Probe exactly one port, once, and only because a user asked.** Scanning
  ports on someone's machine is the attack this must stay distinguishable
  from.
- **`/ping` is validated, never trusted.** Any program can listen on that
  port. Treat what it says as a claim.
- Only `POST /v1/actions/pull` is ever actionable cross-origin, and only when
  `bridge_mode` is `full`.

---

## Guardrails

- **Never treat an alias as a digest.** `GET /v4/genomes/{digest}` takes a
  digest only; passing `hg38` there is a `404`. `GET /v4/genomes?alias=hg38`
  returns an **empty page**, not a `404`, so an empty `items` array there means
  "no such alias", not "no such genome".
- **Never construct a download path by convention.** Read it from
  `GET /v4/assets/{asset_digest}/files`.
- **Never guess a filename from a seek key.** Read the seek key's `value`.
- **Do not page through a whole collection to find one thing.** Use `?q=`, an
  exact filter, or `GET /v4/aliases/{name}`. `limit` is capped at `1000`.
- **Send a real `User-Agent`.** `api.refgenie.org` is behind Cloudflare, which
  rejects default Python user agents (`python-requests/*`, `Python-urllib/*`)
  with a `403`. `curl` is fine as-is.
- **Follow redirects on downloads** (`curl -L`).
- **This API is read-only in server mode.** There is no write surface on a
  public instance — no POST, no PUT, no DELETE. Do not look for one.

## Troubleshooting

- **A page of HTML instead of JSON or Markdown.** You hit the SPA catch-all.
  The path is not an API route on this origin; check `/service-info` for
  `api_base`.
- **`404` from `GET /v4/genomes/{digest}` when you passed a name.** That route
  takes a digest. Resolve the name first with `GET /v4/aliases/{name}`.
- **`404` on a download.** Either the asset is not staged in that mode — check
  its `serving_modes`, or `GET /v4/staged_assets?asset_digest=<DIGEST>` — or
  the file path is not in the asset's `files` list.
- **`422` on a list request.** `limit` is outside `1`–`1000`, `offset` is
  negative, or `search_fields` names a field that endpoint cannot search.
- **`403 missing_action_header`.** A local-dash mutation without
  `X-Refgenie-Action`.
- **`421 forbidden_host`.** You reached a local dash with a non-loopback
  `Host` header.
- **Empty body on a download.** You did not follow the `307`.

## Where to read more

- The refgenie tool, CLI and Python API: <https://docs.refgenie.org/SKILL.md>
- Human documentation: <https://docs.refgenie.org>
- This instance's OpenAPI schema: `GET /openapi.json`
- Source: <https://github.com/refgenie/refgenie>
