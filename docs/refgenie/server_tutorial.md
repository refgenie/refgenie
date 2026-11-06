# Run a server

A refgenie server shares your assets with other people. Clients subscribe to it and `refgenie pull` from it, the same way they pull from the public server at <https://api.refgenie.org>. Run your own to share private or custom assets inside your group, or to publish assets for the community.

This page walks through the steps. For flags, endpoints and environment variables, see the [server reference](server/README.md).

## 1. Install

```bash
pip install "refgenie[server]>=1.0.0a1"
```

## 2. Build and stage assets

The server serves only staged assets. Staging puts a servable copy of each asset (a symlink, a `.tgz` archive, or both, depending on the asset class) in the stage folder, `~/.refgenie/archives` by default.

```bash
refgenie build hg38/fasta --stage
refgenie stage add hg38/bowtie2_index mm10/fasta
refgenie stage list
```

See [Build assets](build.md) and [Stage assets for serving](staging.md).

## 3. Start the server

```bash
refgenie serve
```

The server listens on all network interfaces, port 8000. Open <http://localhost:8000/docs> for the interactive API docs. The server has no web UI of its own; to browse a catalog in a browser, run [`refgenie dash`](dash.md). Use `--port` to pick another port.

## 4. Pull from it

On a client machine, subscribe to your server and pull:

```bash
refgenie subscribe -s http://your-server:8000
refgenie pull hg38/fasta
```

## Run it in production

For a server other people rely on:

- **Use PostgreSQL** instead of the default SQLite database. See [Database backends](database.md).
- **Put it behind a reverse proxy** such as nginx with HTTPS. Do not use `--reload` in production.
- **Host files on cloud storage.** Push staged assets to an S3 bucket or web server with `refgenie push`. When an asset has been pushed to an `http` or `https` remote, the server sends download requests there instead of serving the file itself. See [Configure remote storage](remotes.md) and [Remotes on a server](server/remotes.md).
- **Serve data channels.** The server re-serves the data channels listed in its config file, so clients can get recipes and asset classes from it. See [Serve data channels](server/data_channels.md).

### Docker

The refgenie repository includes a `Dockerfile` that builds a server image:

```bash
git clone https://github.com/refgenie/refgenie.git
cd refgenie
docker build -t refgenie-server .
```

The image reads its database settings from [`deployment/config.yaml`](https://github.com/refgenie/refgenie/blob/master/deployment/config.yaml), a PostgreSQL config filled in from environment variables. It listens on port 80:

```bash
docker run --rm -d -p 8000:80 \
    -e POSTGRES_HOST=db.example.org \
    -e POSTGRES_DB=refgenie \
    -e POSTGRES_USER=refgenie_user \
    -e POSTGRES_PASSWORD=your_secure_password \
    --name refgenie-server \
    refgenie-server
```

Point it at the database you built and staged into, or have it import a published catalog with `REFGENIE_CATALOG_URL` (see the [server reference](server/README.md#environment-variables)).

## Troubleshooting

- **The server will not start.** Install the extras with `pip install "refgenie[server]>=1.0.0a1"`, check that the port is free, and check your database config.
- **Assets do not show up.** Check that they are staged (`refgenie stage list`) and that the server uses the same database and stage folder you staged into.
- **Downloads fail.** Check that the staged files exist in the stage folder, or that the asset was pushed to the remote. The server log names the missing path.

## See also

- [Server reference](server/README.md)
- [Use the dashboard](dash.md), the local web UI for your own refgenie
- [Public servers](servers.md)
