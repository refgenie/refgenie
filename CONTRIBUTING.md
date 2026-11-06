# Contributing to refgenie

Thanks for your interest in refgenie. This document covers how to propose a
change, what a reviewable pull request looks like, and how to run the tests.

Refgenie 1.0 is currently in **alpha**. Interfaces, the database schema, and the
REST API may change without a deprecation period. Please open an issue before
starting significant work, so you don't build on something that is about to
move.

## Getting set up

See [docs/refgenie/development.md](./docs/refgenie/development.md) for the developer quickstart:
prerequisites, `uv sync`, running from source, and the task runner targets.

## Branches and pull requests

- Development happens in [refgenie/refgenie](https://github.com/refgenie/refgenie).
  The default branch is `master`. Open pull requests against it unless an issue
  says otherwise.
- Branch from an up-to-date `master`. Name branches for their content, like
  `fix-alias-resolution` or `add-drs-checksums`.
- Keep a pull request to one logical change. If you find yourself writing "and
  also" in the description, it is probably two pull requests.
- Write commit subjects in the imperative mood, describing the effect of the
  change ("Guard the migration chain against silent drift", not "changes to
  alembic").
- Rebase rather than merge to bring your branch up to date.

A pull request is ready for review when:

- The unit and component tiers pass (see below).
- `task reformat` has been run, and `ruff check .` is clean.
- New behavior has a test. Bug fixes have a test that fails without the fix.
- User-visible changes are reflected in the README, `docs/`, or both.
- [CHANGELOG.md](./CHANGELOG.md) has an entry under "Unreleased" if the change
  is user-visible.
- Frontend changes (`frontend/`) have `npm run typecheck` and `npm run lint`
  clean (`task web-check` runs both). If `package.json` changed,
  regenerate `package-lock.json` and commit it in the same commit, since CI
  runs `npm ci`.

## Running the tests

The suite has four tiers:

| Tier | Command | Needs |
|------|---------|-------|
| `unit` | `pytest` | nothing |
| `component` | `pytest -m component` | nothing |
| `e2e` | `pytest -m e2e` | `samtools` on PATH |
| `integration` | `./tests/scripts/test-integration.sh` | Docker, bulker |

Run integration tests only through the script, which starts and stops the
PostgreSQL container, the HTTP data channel, and the bulker crate. Run it
before submitting anything that touches the server, data channels, or asset
building. See [tests/README.md](./tests/README.md) for which tier a new test
belongs in.

To build a wheel, install it in a clean venv, and smoke-test it:

```bash
./tests/scripts/test-package.sh
```

## Things to know before you change data-layer code

`refgenie/db/events.py` registers SQLAlchemy listeners that delete files from
disk, so a plain `session.delete(...); session.commit()` removes files once the
`COMMIT` succeeds. The rule everywhere in the data layer is
**commit the catalog first, then clean up the filesystem, and make the cleanup
safe to re-run**. Read
["Filesystem side effects of the ORM"](./docs/refgenie/development.md#filesystem-side-effects-of-the-orm)
first.

Schema changes need an alembic revision:

```bash
task alembic-revision MESSAGE="<msg>" DB_CONN_STR=<url>
```

`TestAlembicChain` in `tests/db/test_db.py` checks that the migrations match
the models.

## Publishing a release (maintainers)

Releases are published by `.github/workflows/publish.yml`, which runs when a
GitHub Release is *published*. Every release, pre-release or full, goes to
PyPI. A pre-release version such as `1.0.0a1` is skipped by
`pip install refgenie` unless the user passes `--pre` or pins it. The upload
uses PyPI Trusted Publishing (OIDC); there are no API tokens in this repository.

Trusted publishing is scoped to the `refgenie/refgenie` repository, so releases
are tagged on its `master` branch. The `refgenie` project on pypi.org needs a
trusted publisher (owner `refgenie`, repository `refgenie`, workflow
`publish.yml`, environment `pypi`), added once under the project's Publishing
settings.

The web UI is built by CI, not by hand: `publish.yml`'s `build` job runs
`.github/workflows/build-frontend.yml` and then
`tests/scripts/check-wheel-assets.py --require-fresh`, which fails the release
unless the bundle's `build-info.json` reports the release commit and a clean
frontend checkout (`dirty: false`). There is no path to a release wheel built
from a maintainer's local frontend build.

To cut a release:

1. Bump `version` in `pyproject.toml` and run `uv lock`.
2. Move the "Unreleased" entries in `CHANGELOG.md` under the new version.
3. Tag the release commit on `master` with `v<version>` and push the tag.
4. Create a GitHub Release for that tag. Tick "Set as a pre-release" for alpha
   and beta versions.
