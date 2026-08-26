"""Tests for alias management."""

from tests.e2e.cli_compat.conftest import assert_exit_ok, assert_exit_error, assert_in_output


class TestAlias:
    def test_alias_set_additional(self, runner_with_genome):
        """Set an additional alias for the same genome digest."""
        runner = runner_with_genome
        # Get the digest for "demo"
        result = runner.id("demo")
        assert_exit_ok(result)
        digest = result.stdout.strip()
        assert len(digest) > 0

        # Set a new alias pointing to the same digest
        result = runner.alias_set("hg38", digest)
        assert_exit_ok(result)

        # Get the new alias
        result = runner.alias_get("hg38")
        assert_exit_ok(result)
        assert_in_output(result, digest)

    def test_alias_remove(self, runner_with_genome):
        """Set an alias, remove it, then verify it is gone."""
        runner = runner_with_genome
        # Get digest
        result = runner.id("demo")
        assert_exit_ok(result)
        digest = result.stdout.strip()

        # Set a new alias
        result = runner.alias_set("toremove", digest)
        assert_exit_ok(result)

        # Remove it
        result = runner.alias_remove("toremove")
        assert_exit_ok(result)

        # Both CLIs exit 0 on a missing alias with an empty result:
        # neither the digest nor the removed alias name may appear.
        result = runner.alias_get("toremove")
        assert_exit_ok(result)
        combined = result.stdout + result.stderr
        assert digest not in combined
        assert "toremove" not in combined

    def test_alias_get_nonexistent(self, initialized_runner):
        """Get an alias that was never set.

        Both CLIs deliberately return exit 0 with an empty result on a miss
        (observed 2026-08: Python prints an empty table, Rust prints nothing).
        """
        result = initialized_runner.alias_get("nonexistent_alias_xyz")
        assert_exit_ok(result)
        assert "nonexistent_alias_xyz" not in (result.stdout + result.stderr)

    def test_alias_remove_nonexistent(self, initialized_runner):
        """Remove an alias that does not exist. Should fail gracefully."""
        result = initialized_runner.alias_remove("never_existed_alias")
        assert_exit_error(result)
