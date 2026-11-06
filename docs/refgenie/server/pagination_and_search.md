# Pagination and search

Every server endpoint that lists things (such as `/v4/genomes`, `/v4/asset_groups` and `/v4/recipes`) is paged and can be filtered in two ways:

1. **Exact-match filters**: parameters specific to each endpoint, such as `digest=abc123` or `name=example`
2. **Search**: the `q`, `search_fields` and `operator` parameters

You can use both at once: exact-match filters apply first, then search.

## Pagination

Pagination follows the [GA4GH page-based style](https://github.com/ga4gh/TASC/blob/main/recommendations/API%20pagination%20guide.md). Clients retrieve data in pages, with these query parameters:

- `page`: The page to return, counted from 0 (default: 0). A page past the last one is a `400`.
- `page_size`: Results per page (default: 100, at most 1000). Over 1000 is a `400`.

Every list has a fixed order that ends in a unique column, so walking the pages returns each row exactly once. Genomes come with the most asset groups first, then by species name, then by digest. Page 0 of an empty list is an empty `results`, not an error. The last page is `(total - 1) // page_size`: stop when `(page + 1) * page_size >= total`.

## Search

Search takes these query parameters:

- `q`: The search query string.
- `search_fields`: Comma-separated list of fields to search within. If not specified, searches across all searchable fields for the endpoint.
- `operator`: How to match (default `contains`):
  - `eq`: Exact match (field equals the query string exactly).
  - `contains`: Checks if the field contains the query string (case-insensitive).
  - `starts_with`: Checks if the field starts with the query string (case-insensitive).
  - `ends_with`: Checks if the field ends with the query string (case-insensitive).

### Multiple Search Fields

With several fields in `search_fields`, a row matches if any one of them matches. For example:

- `?q=test&search_fields=name,description` will return items where either the name OR description contains "test"
- `?q=v1&search_fields=name,version&operator=starts_with` will return items where either the name OR version starts with "v1"

## Searchable Fields by Endpoint

- **Genomes** (`/genomes`): `digest`, `description`, `species`, `synonyms`, `assembly_source`, `assembly_accession`, `aliases`
- **Aliases** (`/aliases`): `name`
- **Asset Groups** (`/asset_groups`): `name`
- **Assets** (`/assets`): `name`, `digest`, `path`
- **Asset Classes** (`/asset_classes`): `name`, `version`, `description`
- **Recipes** (`/recipes`): `name`, `version`, `description`
- **Staged assets** (`/staged_assets`): `asset_digest`, `mode`

The genome `species` field matches the genome's scientific name and, when its store publishes a species list, every name of its taxon (scientific, common and other common names). So `brown rat` finds *Rattus norvegicus* even where the listed common name is `Norway rat`. Each name is matched on its own under the usual operators.

The genome `synonyms` field holds other names of the genome itself, from the FHR record's `genomeSynonym` (`GRCh38`). Each name is matched on its own: `eq` must equal one whole name, `starts_with` and `ends_with` test the start and end of each name, and `contains` looks inside them.

In genome responses, `common_name` is the species list's common name for the genome's taxon (genomes have none of their own), and list and detail rows carry a `species` object (`taxon_id`, `scientific_name`, `common_name`, `synonyms`), or `null` when the species list has no row for the genome's taxon.

## Exact-Match Filters by Endpoint

Each endpoint also supports specific exact-match filters:

- **Genomes** (`/genomes`): `digest`, `alias`, `has_assets` (`true` or `false`), `taxon_id` (NCBI taxon id)
- **Aliases** (`/aliases`): `name`, `genome_digest`
- **Asset Groups** (`/asset_groups`): `genome_digest`, `asset_class`, `asset_group_name`, `asset_group_id`
- **Assets** (`/assets`): `name`, `asset_group_name`, `genome_digest`, `recipe_name`, `asset_group_id`
- **Asset Classes** (`/asset_classes`): `name`, `version`
- **Recipes** (`/recipes`): `output_asset_class`, `name`, `version`
- **Archives** (`/archives`): `digest`, `asset_digest`, `genome_digest`

## Example Requests

### Exact-Match Filters

```http
GET /v4/genomes?digest=abc123&alias=hg38
GET /v4/assets?name=fasta&genome_digest=abc123
GET /v4/recipes?output_asset_class=fasta&version=0.0.1
GET /v4/genomes?has_assets=true&taxon_id=9606
```

### Generic Search

```http
GET /v4/recipes?operator=contains&page=0&page_size=100
GET /v4/recipes?q=fasta&search_fields=name&operator=contains&page=0&page_size=100
GET /v4/recipes?q=0.0.1&search_fields=version&operator=contains&page=0&page_size=100
GET /v4/asset_groups?q=fas&search_fields=name&operator=starts_with&page=0&page_size=100
GET /v4/genomes?q=f&search_fields=digest&operator=starts_with&page=1&page_size=100
```

### Combined Approach

```http
GET /v4/genomes?digest=abc123&q=human&search_fields=description&operator=contains
GET /v4/assets?genome_digest=abc123&q=fasta&search_fields=name,path&operator=contains
```

## Example Response Structure

```json
{
  "results": [
    {
      "id": 1,
      "name": "example_recipe",
      "version": "0.0.1"
    },
    {
      "id": 2,
      "name": "another_recipe",
      "version": "0.0.2"
    }
  ],
  "pagination": {
    "page": 0,
    "page_size": 100,
    "total": 150
  }
}
```

