"""Tests for server subscription management."""

from tests.e2e.cli_compat.conftest import assert_exit_ok, assert_in_output


class TestSubscribe:
    def test_subscribe_url(self, initialized_runner):
        """Subscribe to a URL and verify it appears in config."""
        runner = initialized_runner
        result = runner.subscribe("http://srv1.io")
        assert_exit_ok(result)
        result = runner.config_get()
        assert_exit_ok(result)
        assert_in_output(result, "srv1.io")

    def test_subscribe_multiple(self, initialized_runner):
        """Subscribe to two different URLs. Both should appear in config."""
        runner = initialized_runner
        result = runner.subscribe("http://srvA.io")
        assert_exit_ok(result)
        result = runner.subscribe("http://srvB.io")
        assert_exit_ok(result)
        result = runner.config_get()
        assert_exit_ok(result)
        assert_in_output(result, "srvA.io")
        assert_in_output(result, "srvB.io")

    def test_unsubscribe(self, initialized_runner):
        """Subscribe then unsubscribe. URL should no longer appear."""
        runner = initialized_runner
        result = runner.subscribe("http://rmme.io")
        assert_exit_ok(result)
        result = runner.unsubscribe("http://rmme.io")
        assert_exit_ok(result)
        result = runner.config_get()
        assert_exit_ok(result)
        combined = result.stdout + result.stderr
        assert "rmme.io" not in combined

    def test_unsubscribe_nonexistent(self, initialized_runner):
        """Unsubscribe from a URL never subscribed. Should not crash."""
        result = initialized_runner.unsubscribe("http://nosub.io")
        # Should succeed gracefully (Python returns 0)
        assert_exit_ok(result)
