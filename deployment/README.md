# Deployment

Files used to run the public refgenie server on AWS ECS.

| File | What it is | Used by |
|---|---|---|
| `config.yaml` | PostgreSQL settings for the server. Each value is an environment variable (`$POSTGRES_HOST` and so on), filled in when the container starts. | The `Dockerfile` copies it to `/config.yaml` in the image. |
| `task_defs/primary.json` | The ECS task definition: the container, its environment (`REFGENIE_STORES`, `REFGENIE_CATALOG_URL`) and the secrets that supply the `POSTGRES_*` variables. | `.github/workflows/deploy-api.yml` fills in the new image and deploys it to the `refgenie-service-primary` service. |
| `stores.example.json` | An example value for `REFGENIE_STORES`, the RefgetStores the server registers at startup. | Documentation only; nothing reads it. |

Change `config.yaml` and the task definition together. Every variable that
`config.yaml` names must be set in `task_defs/primary.json`, or the server
starts without a database.

## The web UI at refgenie.org

refgenie.org serves the frontend as a static site on a Cloudflare Worker, and
the app reads its data from api.refgenie.org. The Worker config is
`frontend/site/` and the build settings are `frontend/.env.site`. Deploy from
`frontend/` with `npm run deploy:site` (needs `wrangler login`).
