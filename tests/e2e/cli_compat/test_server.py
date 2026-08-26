"""Smoke tests for the serve command."""

import os
import time
import urllib.request

import pytest

from tests.e2e.cli_compat.conftest import find_free_port


class TestServer:
    def test_serve_starts_and_responds(self, initialized_runner):
        """Start serve, wait for it to accept connections, then verify HTTP response."""
        port = find_free_port()
        proc = initialized_runner.serve_start(port)
        try:
            # Wait up to 20 seconds for the server to start (uvicorn can be slow)
            for _ in range(40):
                time.sleep(0.5)
                try:
                    resp = urllib.request.urlopen(f"http://127.0.0.1:{port}/", timeout=2)
                    assert resp.status == 200
                    return  # success
                except Exception:
                    if proc.poll() is not None:
                        # Process died
                        stdout = proc.stdout.read() if proc.stdout else ""
                        stderr = proc.stderr.read() if proc.stderr else ""
                        pytest.fail(
                            f"Server process died before accepting connections.\n"
                            f"stdout: {stdout}\nstderr: {stderr}"
                        )
            pytest.fail("Server did not respond within 20 seconds")
        finally:
            proc.terminate()
            try:
                proc.wait(timeout=5)
            except Exception:
                proc.kill()

    @pytest.mark.skipif(os.geteuid() == 0, reason="port 1 binds as root")
    def test_serve_invalid_port(self, initialized_runner):
        """Serve on a privileged port must exit nonzero with a permission error.

        Both CLIs agree: Python's uvicorn reports "[errno 13] permission
        denied"; the Rust CLI reports "Permission denied (os error 13)".
        """
        proc = initialized_runner.serve_start(1)
        try:
            # uvicorn app startup precedes the bind failure (~2s locally,
            # slower on cold CI). A hang past 20s is a test failure.
            stdout, stderr = proc.communicate(timeout=20)
        finally:
            if proc.poll() is None:
                proc.kill()
                proc.wait(timeout=5)
        assert proc.returncode != 0, (
            f"Expected nonzero exit binding port 1:\nstdout: {stdout}\nstderr: {stderr}"
        )
        assert "permission denied" in (stdout + stderr).lower(), (
            f"Expected a permission-denied message:\nstdout: {stdout}\nstderr: {stderr}"
        )
