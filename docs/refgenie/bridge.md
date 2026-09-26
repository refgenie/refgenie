# Connect the public web UI to your local refgenie

The public refgenie web UI at <https://refgenie.org> can connect to a
`refgenie dash` running on your own computer. When connected, the page shows
which catalog assets you already have on disk ("Local" badges), lets you
browse your local genomes, and can hand a pull directly to your local
refgenie — without you switching to a terminal.

This connection is called the **localhost bridge**. Everything stays on your
machine: the public page talks to `http://localhost:8080` directly from your
browser, and no data about your local instance is sent anywhere else.

## Using it

1. Start your dashboard locally:

        refgenie dash

2. Open <https://refgenie.org>, and click **Connect local refgenie** (in
   the sidebar, or the card on the Genomes page).
3. Confirm the port (default `8080`, matching `refgenie dash`) and click
   **Connect**.

Chrome (version 142 and later) will show a one-time permission prompt asking
whether the site may access devices on your local network — that prompt is
about your own `refgenie dash`, and the connection only works if you allow
it. If you denied it by accident, re-enable it from the site settings icon in
the address bar.

Once connected, the page remembers the connection and re-checks it silently
on later visits. Click **Disconnect** in the sidebar to forget it.

## Bridge modes

How much the public page may do is controlled by the **bridge mode** on your
local refgenie:

| Mode | What an allowlisted page may do |
|---|---|
| `off` | Nothing. No cross-origin access at all. |
| `read` (default) | See that your refgenie is running, and read what it has (genomes, assets, job status). |
| `full` | Additionally submit a **pull** ("Pull to my refgenie" on asset pages). |

Set it for one run:

    refgenie dash --bridge full

or persistently with the environment variable `REFGENIE_BRIDGE_MODE`.

Even under `full`, only pull is available cross-origin, and only from a
server your refgenie subscribes to. For any other server, the page offers
"Open in local refgenie" so you can confirm the pull there. Deleting assets,
managing aliases, building, and every other state change can only be done
from the local dashboard itself (`http://localhost:8080`) — the remote page
links you there instead.

Related settings:

- `REFGENIE_BRIDGE_ORIGINS` — comma-separated origins allowed to connect
  (default `https://refgenie.org`).
- `REFGENIE_BRIDGE_ORIGIN_REGEX` — optional regex for extra origins (for UI
  development and preview deployments). **A careless regex here disables the
  origin allowlist entirely** — leave it unset unless you know why you need
  it.
- `REFGENIE_BRIDGE_EXPOSE_PATHS` — default `false`; when `true`, the
  connection handshake reveals your genome folder path to allowlisted
  origins (paths can leak your OS username, so this is off by default).

## Safari

Safari blocks HTTPS pages from talking to servers on your own computer
([WebKit bug 171934](https://bugs.webkit.org/show_bug.cgi?id=171934)), so
the bridge does not work there and the page will say so instead of trying.
Open your local dashboard directly at `http://localhost:8080` — it runs the
same interface.

## Troubleshooting

If the page cannot connect, check:

- Is `refgenie dash` actually running?
- Did you allow the browser's local-network permission prompt?
- Is the dash on a different port? Enter it in the connect dialog.
- Is your refgenie too old? The bridge needs a version that serves `/ping`.
- Is the bridge disabled? Check `REFGENIE_BRIDGE_MODE` (or `--bridge`).

## Security notes (what the bridge does and does not protect)

The bridge is designed so that visiting a random website never gives that
website control over your local refgenie:

- Only allowlisted origins can read anything (CORS); everything else is
  refused, and a `Host`-header guard blocks DNS-rebinding tricks.
- State-changing requests additionally require a custom header (which forces
  a browser preflight) and are refused cross-origin unless you opted into
  `--bridge full` — and even then, only pull.
- Remote browsing and bridge pulls accept only servers you subscribed to.
  Pulling from any other server goes through the local `/pull` page, where
  you confirm it yourself.
- No cookies, tokens, or credentials exist anywhere in the exchange.

Two honest non-goals to be aware of:

- **The local dashboard is unauthenticated.** The trust boundary is your
  machine: anyone with shell access to it controls your refgenie regardless
  of bridge settings. Do not run `refgenie dash` on a shared multi-user
  machine.
- **The page cannot verify what answers on the port.** Any local process
  could impersonate a refgenie; local data is therefore always displayed
  clearly labeled as local, never silently merged into catalog results.
