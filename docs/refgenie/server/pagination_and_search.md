# REST API Pagination and Search

## Overview

The REST API supports robust pagination and search functionalities to efficiently handle large datasets. These features are implemented using FastAPI, leveraging its dependency injection system for code reuse and modularity. All endpoints that _list_ resources (e.g., `/recipes`, `/asset_groups`, `/genomes`) support both specific query parameters and generic search capabilities, operating in a consistent manner. More details below.

## Filtering Options

The API provides two complementary approaches for filtering data:

1. **Legacy Query Parameters**: Endpoint-specific parameters for exact matches (e.g., `digest=abc123`, `name=example`)
2. **Generic Search**: Flexible search functionality using `q`, `search_fields`, and `operator` parameters

Both approaches can be used together - legacy filters are applied first, then generic search is applied to the filtered results.

### Pagination

Pagination allows clients to retrieve data in chunks, reducing the load on the server and improving response times. The following query parameters are used for pagination:

- `offset`: The starting point of the data to retrieve (default: 0).
- `limit`: The maximum number of items to return (default: 100, range: 1–1000).

### Search

Search functionality enables filtering and sorting of data based on specific criteria. The following query parameters are supported:

- `q`: The search query string.
- `search_fields`: Comma-separated list of fields to search within. If not specified, searches across all searchable fields for the endpoint.
- `operator`: The search operator. Supported operators include:
  - `eq`: Exact match (field equals the query string exactly).
  - `contains`: Checks if the field contains the query string (case-insensitive).
  - `starts_with`: Checks if the field starts with the query string (case-insensitive).
  - `ends_with`: Checks if the field ends with the query string (case-insensitive).

#### Multiple Search Fields

When multiple fields are specified in `search_fields` (comma-separated), the search operates with an OR logic - meaning the query will match if ANY of the specified fields satisfies the search condition. For example:

- `?q=test&search_fields=name,description` will return items where either the name OR description contains "test"
- `?q=v1&search_fields=name,version&operator=starts_with` will return items where either the name OR version starts with "v1"

### Searchable Fields by Endpoint

- **Genomes** (`/genomes`): `digest`, `description`, `species_name`, `aliases`
- **Aliases** (`/aliases`): `name`
- **Asset Groups** (`/asset_groups`): `name`
- **Assets** (`/assets`): `name`, `digest`, `path`
- **Asset Classes** (`/asset_classes`): `name`, `version`, `description`
- **Recipes** (`/recipes`): `name`, `version`, `description`
- **Archives** (`/archives`): `digest`, `asset_digest`

### Legacy Query Parameters by Endpoint

Each endpoint also supports specific legacy query parameters for exact matching:

- **Genomes** (`/genomes`): `digest`, `alias`
- **Aliases** (`/aliases`): `name`, `genome_digest`
- **Asset Groups** (`/asset_groups`): `genome_digest`, `asset_class`, `asset_group_name`, `asset_group_id`
- **Assets** (`/assets`): `name`, `asset_group_name`, `genome_digest`, `recipe_name`, `asset_group_id`
- **Asset Classes** (`/asset_classes`): `name`, `version`
- **Recipes** (`/recipes`): `output_asset_class`, `name`, `version`
- **Archives** (`/archives`): `digest`, `asset_digest`

### Example Requests

#### Legacy Query Parameters

```http
GET /v4/genomes?digest=abc123&alias=hg38
GET /v4/assets?name=fasta&genome_digest=abc123
GET /v4/recipes?output_asset_class=fasta&version=0.0.1
```

#### Generic Search

```http
GET /v4/recipes?sort_order=asc&operator=contains&offset=0&limit=100
GET /v4/recipes?q=fasta&search_fields=name&sort_order=asc&operator=contains&offset=0&limit=100
GET /v4/recipes?q=0.0.1&search_fields=version&sort_order=asc&operator=contains&offset=0&limit=100
GET /v4/asset_groups?q=fas&search_fields=name&sort_order=asc&operator=starts_with&offset=0&limit=100
GET /v4/genomes?q=f&search_fields=digest&sort_order=asc&operator=starts_with&offset=0&limit=100
```

#### Combined Approach

```http
GET /v4/genomes?digest=abc123&q=human&search_fields=description&operator=contains
GET /v4/assets?genome_digest=abc123&q=fasta&search_fields=name,path&operator=contains
```

### Example Response Structure

```json
{
  "items": [
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
    "offset": 0,
    "limit": 100,
    "total": 150
  }
}
```

### Dependency Injection

FastAPI's dependency injection system is utilized to promote code reuse and maintainability. This approach ensures that shared logic, such as pagination and search handling, can be modularized and reused across multiple endpoints.