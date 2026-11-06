# Use the dashboard

The refgenie dashboard is a local web interface for your own refgenie instance. You can browse your genomes and assets, pull assets from subscribed servers, build new assets, and manage aliases and subscriptions, all from the browser. To share assets with other machines, run a server with [`refgenie serve`](server_tutorial.md) instead.

## Install

The dashboard requires the `dash` extras:

```bash
pip install --pre "refgenie[dash]"
```

## Launch the dashboard

```bash
refgenie dash
```

This starts the dashboard at <http://127.0.0.1:8080> and opens it in your browser. If the page is blank, refresh after a moment.

Options:

| Option | Meaning |
|---|---|
| `-p`, `--port PORT` | Port to run on (default `8080`) |
| `-b`, `--bridge off\|read\|full` | How much the public web UI may do with this dashboard (see [Connect the web UI to your local refgenie](bridge.md)). Overrides `$REFGENIE_BRIDGE_MODE` for this run. |

```bash
refgenie dash -p 9090
```

## What the dashboard shows

The top bar links to each page.

### Browse

- **Home**: an overview of your instance.
- **Genomes**: every genome you have, with its aliases and assets. Click a genome to see its assets, or an asset to see its seek keys, file paths, and digest.
- **Species**: your genomes grouped by species.
- **Asset classes** and **Recipes**: the asset types refgenie knows about and the recipes that build them.
- **Aliases**: the human-readable names that point at each genome.
- **Servers**: the servers your refgenie pulls from, and whether each one is reachable right now. Subscribe to a new server or unsubscribe from one here. A new install is already subscribed to `https://api.refgenie.org`. If a pull fails because you are not subscribed to the right server, the error links here with the server filled in, ready for you to click **Subscribe**.

### Act

- **Remote**: browse assets on the servers you subscribe to, see which ones you already have, and pull the ones you don't.
- **Build**: pick a genome and a recipe, fill in the inputs, and build an asset.
- **Manage**: add or remove data channels, edit aliases, view registered recipes and asset classes, and initialize a new genome. Register new recipes and asset classes from the command line.
- **Pull**: a confirmation page for a pull handed over from another site, such as the public web UI. The pull starts only when you click to confirm.

You can also delete a genome or an asset from its detail page.

### Jobs

Pulls and builds can take a while, so they run as background **jobs**. When you start one, a dock opens at the bottom of the page and shows live progress and logs for each running job. You can keep browsing while it works.

The **Jobs** page lists recent jobs with their status and duration; click one for details. The list lives in the dashboard process, so it starts empty each time you restart `refgenie dash`.

## Connect from the public site

The public web UI at <https://refgenie.org> can connect to your running dashboard to show which assets you already have and to send pulls to it. See [Connect the web UI to your local refgenie](bridge.md).

## Dashboard vs. server

| Feature | `refgenie dash` | `refgenie serve` |
|---|---|---|
| Purpose | Browse and manage your local refgenie instance | Serve assets to remote clients via API |
| Audience | You (local use) | Other machines and users |
| Default port | 8080 | 8000 |
| Install extra | `pip install --pre "refgenie[dash]"` | `pip install --pre "refgenie[server]"` |

!!! warning
    The dashboard has no login. It listens only on `127.0.0.1`, but anyone with an account on your machine can use it to change your refgenie, so run it only on a computer you do not share.
