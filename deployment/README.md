# Deployment

Files used to run the public refgenie server on AWS ECS.

| File | What it is | Used by |
|---|---|---|
| `config.yaml` | PostgreSQL settings for the server. Each value is an environment variable (`$POSTGRES_HOST` and so on), filled in when the container starts. | The `Dockerfile` copies it to `/config.yaml` in the image. |
| `task_defs/primary.json` | The ECS task definition: the container, its environment (`REFGENIE_COLLECTIONS`, `REFGENIE_CATALOG_URL`) and the secrets that supply the `POSTGRES_*` variables. | `.github/workflows/deploy-api.yml` fills in the new image and deploys it to the `refgenie-service-primary` service. |
| `collections.example.json` | An example value for `REFGENIE_COLLECTIONS`, the collection pages the web UI shows (each a branded slice of the catalog drawn from one or more stores). | Documentation only; a server test checks that it is valid. |

Change `config.yaml` and the task definition together. Every variable that
`config.yaml` names must be set in `task_defs/primary.json`, or the server
starts without a database.

The server learns which RefgetStores to read from the catalog at
`REFGENIE_CATALOG_URL`: it imports the catalog at startup, before it reads its
store list, and again every `CATALOG_IMPORT_INTERVAL_SECONDS`. When an import
changes the list, the server starts serving the new stores and stops serving
removed or disabled ones without a restart. To publish a store, add it in the
registry (refgenie-registry `stores/`), not here.

`REFGENIE_PUBLIC_URL` (default `https://api.refgenie.org`) is the base URL the MCP tools put in download links and in the server sentence.

## The web UI at refgenie.org

refgenie.org serves the frontend as a static site on a Cloudflare Worker, and
the app reads its data from api.refgenie.org. The Worker config is
`frontend/site/` and the build settings are `frontend/.env.site`. Deploy from
`frontend/` with `npm run deploy:site` (needs `wrangler login`).

Cloudflare's Browser Integrity Check must stay off for the refgenie.org zone.
With it on, Python's built-in urllib (and libwww-perl) gets a 403 with
`error code: 1010` before the request reaches the origin.
