<!--
Thanks for contributing. See CONTRIBUTING.md for branch and review conventions.
-->

## What this changes

<!-- One or two sentences. What behavior is different after this merges? -->

## Why

<!-- Link the issue if there is one: "Fixes #123". Otherwise explain the motivation. -->

## How it was verified

<!--
Which suites did you run, and what did you add? Be specific:
- unit: `pytest tests/ --ignore=tests/integration`
- integration: `./tests/scripts/test-integration.sh` (needed for server,
  data channel, or asset-building changes)
- manual steps, if any
-->

## Checklist

- [ ] Unit tests pass
- [ ] Integration tests pass, or this change cannot affect them
- [ ] `ruff check .` is clean and `task reformat` has been run
- [ ] New behavior has a test; a bug fix has a test that fails without it
- [ ] README / `docs/` updated for user-visible changes
- [ ] `CHANGELOG.md` updated under "Unreleased" for user-visible changes
- [ ] Schema changes include an alembic revision
