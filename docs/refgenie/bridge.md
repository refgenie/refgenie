# Connect the public web UI to your local refgenie

The public web UI at <https://refgenie.org> can connect to a `refgenie dash` running on your own computer. Once connected, the page shows which catalog assets you already have (a "Local" column in the genome and asset tables), lets you browse your local genomes, and can send a pull to your local refgenie without a trip to the terminal.

This connection is called the **localhost bridge**. Your browser talks to `http://localhost:8080` directly, so information about your local refgenie stays on your machine.

## Use it

1. Start your dashboard:

        refgenie dash

2. Open <https://refgenie.org> and click **Connect local refgenie** in the top bar, or **Connect** on the "Using refgenie locally?" card on the home page.
3. Confirm the port (default `8080`) and click **Connect**.

Chrome (version 142 and later) asks once whether the site may access devices on your local network. Allow it; that is your own `refgenie dash`. If you denied it by accident, turn it back on from the site settings icon in the address bar.

The page remembers the connection and quietly re-checks it on later visits. While connected, the top bar shows a **connected** badge, a **Browse local** link, and a **Disconnect** button that forgets the connection.

**Browse local** (or the card on the home page) shows your genomes, species, asset classes, recipes, aliases, subscribed servers and job history. These views are read-only and carry a "local refgenie" banner. To subscribe to a server, the **Servers** view links to the same page on your local dashboard, where you confirm it.

A pull you start with **Pull to my refgenie** appears in a **Pulls to your local refgenie** panel at the bottom of the page. It keeps showing progress as you move around the site and after a reload, and links to the asset when it finishes. **Disconnect** clears the list.

## Bridge modes

The **bridge mode** on your local refgenie sets how much the public page may do:

| Mode | What an allowlisted page may do |
|---|---|
| `off` | Nothing. |
| `read` (default) | See that your refgenie is running, and read what it has (genomes, assets, job status). |
| `full` | Also start a **pull** ("Pull to my refgenie" on asset pages). |

Set it for one run with `refgenie dash --bridge full`, or for every run with the environment variable `REFGENIE_BRIDGE_MODE`.

In `full` mode the page can pull only from servers your refgenie subscribes to. For any other server, the page offers **Open in local refgenie**, where you confirm the pull, and a link to your dashboard's **Servers** page with that server filled in. Subscribing, deleting, managing aliases, building and every other change happen on the local dashboard itself (`http://localhost:8080`), and the public page links you there.

Other settings:

- `REFGENIE_BRIDGE_ORIGINS`: comma-separated origins allowed to connect (default `https://refgenie.org,https://api.refgenie.org`).
- `REFGENIE_BRIDGE_ORIGIN_REGEX`: an optional regex for extra origins, for UI development and preview sites. A loose regex lets any matching site connect, so leave it unset unless you need it.
- `REFGENIE_BRIDGE_EXPOSE_PATHS`: default `false`. When `true`, the connection handshake shows allowlisted sites your genome folder path, which can include your OS username.

## Safari

Safari blocks HTTPS pages from talking to servers on your own computer ([WebKit bug 171934](https://bugs.webkit.org/show_bug.cgi?id=171934)), so the page tells you the bridge is unavailable. Open your dashboard at `http://localhost:8080` instead; it runs the same interface.

## Troubleshooting

If the page cannot connect, check that:

- `refgenie dash` is running.
- You allowed the browser's local-network prompt.
- The port in the connect dialog matches your dashboard's port.
- The bridge is on: `REFGENIE_BRIDGE_MODE` (or `--bridge`) is `read` or `full`.

## Security

The bridge keeps a website you visit from taking control of your local refgenie:

- Only allowlisted origins can read anything (CORS). A `Host`-header check blocks DNS-rebinding tricks.
- Requests that change anything need a custom header, which forces a browser preflight. Cross-origin, they are refused unless you chose `--bridge full`, and even then only pull is allowed.
- Remote browsing and bridge pulls accept only servers you subscribe to. A pull from any other server goes through the local `/pull` page, where you confirm it.
- File contents stay off the bridge. The local dashboard can show and download an asset's files (`/v1/files`) for itself only; the public site is refused in every mode.
- No cookies, tokens or credentials are involved.

Two limits to keep in mind:

- **The local dashboard has no login.** Anyone with shell access to your computer controls your refgenie, whatever the bridge settings. Run `refgenie dash` only on a computer you do not share.
- **The page cannot check what answers on the port.** Another local program could pretend to be refgenie, so local data is always labeled as local and never mixed into catalog results.
