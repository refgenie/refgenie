# Use the dashboard

The refgenie dashboard is a local web interface for browsing your genomes, aliases, assets, and remote assets from subscribed servers. It is distinct from `refgenie serve`, which exposes an API for remote clients.

## Install

The dashboard requires the `dash` extras:

```bash
pip install refgenie[dash]
```

## Launch the dashboard

```bash
refgenie dash
```

This starts the dashboard on port 8080 and automatically opens your browser. If the page doesn't load immediately, refresh after a moment.

### Custom port

```bash
refgenie dash -p 9090
```

## What the dashboard shows

The dashboard provides a web UI with the following pages:

### Index page

The landing page lists all genomes in your refgenie instance. Each genome shows its digest, aliases, and a summary of available assets. Click a genome to view its details.

### Genome splash page

Shows details for a single genome, including:

- Genome digest and all associated aliases
- List of all available assets for this genome
- Asset metadata (class, tag, description)

Click an asset to see its full details.

### Asset detail page

Displays metadata and file paths for a specific asset:

- Asset class and tag
- Build date and description
- Seek keys and their resolved file paths
- Asset digest for identity verification

### Manage page (`/page/manage`)

The administration interface provides:

- **Configuration overview**: Current database backend, genome folder, and stage folder paths
- **Remote assets**: Browse assets available on subscribed servers that have not yet been pulled locally
- **Pull actions**: Pull remote assets directly from the manage page interface
- **Server subscriptions**: View and manage your subscribed asset servers

## Dashboard vs. server

| Feature | `refgenie dash` | `refgenie serve` |
|---|---|---|
| Purpose | Browse your local refgenie instance | Serve assets to remote clients via API |
| Audience | You (local use) | Other machines and users |
| Default port | 8080 | 8000 |
| Install extra | `pip install refgenie[dash]` | `pip install refgenie[server]` |
