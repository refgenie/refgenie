# Contributing to refgenie

Thanks for your interest in refgenie. This document covers how to propose a
change, what a reviewable pull request looks like, and how to run the tests.

Refgenie 1.0 is currently in **alpha**. Interfaces, the database schema, and the
REST API may change without a deprecation period. Please open an issue before
starting significant work, so you don't build on something that is about to
move.

## Getting set up

See [docs/development.md](./docs/development.md) for the developer quickstart:
prerequisites, `uv sync`, running from source, and the task runner targets.

## Branches and pull requests

- The default branch is `master`. Open pull requests against it unless an issue
  says otherwise.
- Branch from an up-to-date `master`. Name branches for their content —
  `fix-alias-resolution`, `add-drs-checksums`. There is no enforced prefix
  scheme.
- Keep a pull request to one logical change. If you find yourself writing "and
  also" in the description, it is probably two pull requests.
- Write commit subjects in the imperative mood, describing the effect of the
  change rather than the mechanics — see `git log` for the house style
  ("Guard the migration chain against silent drift", not "changes to alembic").
- Rebase rather than merge to bring your branch up to date.

A pull request is ready for review when:

- The unit suite passes (see below).
- `task reformat` has been run, and `ruff check .` is clean.
- New behavior has a test. Bug fixes have a test that fails without the fix.
- User-visible changes are reflected in the README, `docs/`, or both.
- [CHANGELOG.md](./CHANGELOG.md) has an entry under "Unreleased" if the change
  is user-visible.
- Frontend changes (`frontend/`) have `npm run typecheck` and `npm run lint`
  clean (`task web-check` runs both). If `package.json` changed,
  `package-lock.json` is regenerated and committed in the **same** commit --
  a drifted lockfile fails `npm ci` in CI, not `npm install`.

## Running the tests

**Unit tests** — in-memory SQLite, no external services:

```bash
pytest tests/ --ignore=tests/integration
```

Note that `-s` is set in `addopts` and the suite depends on it; do not remove
it. The suite also prints no terminal summary, because something in the piper
path calls `os._exit`. Use `--junitxml=results.xml` if you need reliable counts.

**Integration tests** — require Docker (PostgreSQL) and bulker:

```bash
./tests/scripts/test-integration.sh
```

Do **not** run `pytest tests/integration/` directly; the script manages the
container lifecycle, the HTTP data channel, and the bulker crate. Run it before
submitting anything that touches the server, the data channel, or asset
building.

**Package validation** — builds a wheel, installs it in an isolated venv, and
smoke-tests it:

```bash
./tests/scripts/test-package.sh
```

## Things to know before you change data-layer code

`refgenie/db/events.py` registers SQLAlchemy listeners that delete files from
disk. A plain `session.delete(...); session.commit()` therefore destroys data —
after the `COMMIT`, never before it. The rule everywhere in the data layer is
**commit the catalog first, then clean up the filesystem, and make the cleanup
safe to re-run**. Read
["Filesystem side effects of the ORM"](./docs/development.md#filesystem-side-effects-of-the-orm)
first.

Schema changes need an alembic revision:

```bash
task alembic-revision MESSAGE="<msg>" DB_CONN_STR=<url>
```

`tests/test_alembic_chain.py` guards the migration chain against drift.

## Publishing a release (maintainers)

Releases are published by `.github/workflows/publish.yml`, which runs when a
GitHub Release is *published*. A pre-release goes to Test PyPI only; a full
release goes to Test PyPI and then to PyPI. Both uploads use PyPI Trusted
Publishing (OIDC) — there are no API tokens in this repository.

Trusted publishing is already configured, scoped to the `refgenie/refgenie`
repository. This code currently lives in `refgenie/refgenie1` and is intended to
move into `refgenie/refgenie`; **until that move happens, the publish workflow
will fail OIDC when run from this repository.**

The web UI is built by CI, not by hand: `publish.yml`'s `build` job runs
`.github/workflows/build-frontend.yml` and then
`tests/scripts/check-wheel-assets.py --require-fresh`, which fails the release
unless the bundle's `build-info.json` reports the release commit and a clean
frontend checkout (`dirty: false`). There is no path to a release wheel built
from a maintainer's local frontend build.

To cut a release:

1. Bump `version` in `pyproject.toml` and run `uv lock`.
2. Move the "Unreleased" entries in `CHANGELOG.md` under the new version.
3. Tag and push: `git tag v<version> && git push origin v<version>`.
4. Create a GitHub Release for that tag. Tick "Set as a pre-release" for alpha
   and beta versions.
