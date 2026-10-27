# ADR: Localhost bridge between the public web UI and a local refgenie dash

**Date:** 2026-08-12
**Status:** Accepted

> **Note (2026-08-20):** the hostnames below reflect the layout at the time of
> the decision. The public UI has since moved to the apex `refgenie.org`
> (`ui.refgenie.org` now redirects there) and the docs to
> `docs.refgenie.org`. The decision itself is unchanged: the UI origin is
> allowlisted for the bridge, the docs origin is not.

## Context

The public catalog UI (https://ui.refgenie.org) and a user's local
`refgenie dash` were two disconnected worlds: the public page could not show
which assets already sit on the user's disk, nor hand a pull to the local
instance. Closing that gap means a public HTTPS page talking to a loopback
HTTP server — exactly the shape of a localhost-probing attack, and a pattern
browsers have spent a decade restricting (mixed-content exemptions for
`http://localhost`, Chrome's Local Network Access permission prompt, WebKit
blocking it outright).

Three decisions from the bridge design have cross-repo, long-lived
consequences and are recorded here. The full design (probe behavior, outcome
classification, threat model) lives in the bridge implementation plan.

## Decision 1 (D1): `/ping` is its own endpoint, not an extension of `/service-info`

The local server exposes a dedicated, unversioned, unprefixed `/ping`
endpoint (present in both server and local modes, `Cache-Control: no-store`)
carrying `service`, an integer `bridge_version`, instance identity, the
bridge mode, and the **same capability key set `/service-info` emits**.

**Rationale.** `/service-info` is a GA4GH-shaped public discovery document
with a stable published meaning (`seqcol.refget_store.url` is a client
bootstrap field). Stuffing browser-handshake fields into it would couple an
internal contract to a standards-adjacent document and force every public
consumer to parse fields that only matter to a web page. The capability
vocabulary, however, is shared — not forked — so the UI gates features
identically whichever document it read.

## Decision 2 (D4): cross-origin scope is browse + presence + pull; everything else deep-links to the local UI

The only action reachable cross-origin is `POST /v1/actions/pull` (plus
job-status reads). Delete, alias, subscribe, build, and recipe/asset-class
changes are same-origin only, enforced server-side by an explicit
allowed-path set (and `DELETE` is absent from the CORS method list
entirely). For everything else the remote page deep-links into the local SPA
(`/genomes/{digest}`, `/pull?...` prefilled and never auto-executed).

**Rationale.**

- Pull is the only verb whose natural trigger lives on the remote page (you
  are looking at a remote asset because you do not have it); every other
  verb operates on data the local SPA already browses better (same-origin,
  no CORS, no permission prompts).
- Blast radius: the worst case of an over-trusted bridge is disclosure plus
  an unwanted download — not data loss.
- Auditability: `full` mode's entire cross-origin contract is one POST.
- Chrome's local-network permission is per-origin and persistent, so an XSS
  on the public UI would inherit bridge access forever; that is survivable
  when the bridge can browse and pull, not when it can delete.

## Decision 3 (D5): reads default-on for allowlisted origins; writes default-off

`REFGENIE_BRIDGE_MODE` has three values — `off` / `read` / `full` — with
default `read`. Cross-origin pull requires an explicit opt-in
(`refgenie dash --bridge full`); the server's 403 names that exact remedy so
the UI can render it verbatim.

**Rationale.** `read` gives the feature its no-friction first-run experience
while confining exposure to origins the user already trusts enough to visit
(the origin allowlist defaults to the deployed UI origins only — the docs
site is deliberately not allowlisted).
`full` being opt-in means no web page can talk a user into a
state-changing configuration.

## Consequences

- The `/ping` contract is versioned by a single integer bumped only on
  breaking changes; feature gating is exclusively via capability flags, so
  UI and refgenie versions can skew indefinitely in either direction.
- One security module (`refgenie/server/local/security.py`) owns CORS, the
  host guard, the action header, and the bridge policy; the bridge added
  settings and one middleware (the Local Network Access preflight header) to
  it rather than a parallel path.
- The local dash remains unauthenticated by design; the documented trust
  boundary is the machine, not the user account (do not run `refgenie dash`
  on shared multi-user machines).
- Safari cannot use the bridge at all (WebKit bug 171934); it gets a
  documented fallback to `http://localhost:8080`, not a workaround.

## Related documents

- Bridge how-to: `docs/refgenie/bridge.md`
- Implementation plan: jot `refgenie1_localhost_bridge_plan_v1.md`
- Cross-plan contracts note: jot `refgenie1/web_ui_contracts.md`
