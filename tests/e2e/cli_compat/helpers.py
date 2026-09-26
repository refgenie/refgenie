"""Plain helpers for the cli_compat e2e suite: exit/output assertions and the CLI runner.

Fixtures stay in ``conftest.py``; anything a test module imports lives here.
"""

import os
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

from tests.helpers import cli_argv, popen_cli


def assert_exit_ok(result: subprocess.CompletedProcess):
    """Assert exit code 0, print stderr on failure."""
    assert result.returncode == 0, (
        f"Command failed (rc={result.returncode}):\n"
        f"stderr: {result.stderr}\nstdout: {result.stdout}"
    )


def assert_exit_error(result: subprocess.CompletedProcess):
    """Assert nonzero exit code."""
    assert result.returncode != 0, f"Expected error but got rc=0:\nstdout: {result.stdout}"


def assert_in_output(result: subprocess.CompletedProcess, text: str):
    """Assert text appears in stdout or stderr."""
    combined = result.stdout + result.stderr
    assert text in combined, (
        f"Expected '{text}' in output:\nstdout: {result.stdout}\nstderr: {result.stderr}"
    )


# ---------------------------------------------------------------------------
# RefgenieRunner
# ---------------------------------------------------------------------------


@dataclass
class RefgenieRunner:
    """Abstraction over Python or Rust refgenie CLI."""

    mode: str  # "python" or "rust"
    binary: str  # path to binary or "python"
    db_path: Path
    genome_folder: Path
    stage_folder: Path
    home_path: Path
    env: dict = field(default_factory=dict)

    def _build_env(self) -> dict:
        """Build environment dict for subprocess calls."""
        env = {**os.environ, **self.env}
        if self.mode == "python":
            env["REFGENIE_HOME_PATH"] = str(self.home_path)
            env["REFGENIE_GENOME_FOLDER"] = str(self.genome_folder)
            env["REFGENIE_GENOME_STAGE_FOLDER"] = str(self.stage_folder)
        elif self.mode == "rust":
            # Add the directory containing the Rust binary to PATH so that
            # recipe shell templates can find refgenie-build-fasta etc.
            bin_dir = str(Path(self.binary).resolve().parent)
            current_path = env.get("PATH", "")
            env["PATH"] = f"{bin_dir}:{current_path}" if current_path else bin_dir
        return env

    def _run(
        self, args: list, check=False, input_text=None, timeout=60
    ) -> subprocess.CompletedProcess:
        """Run the CLI with args. Returns CompletedProcess."""
        if self.mode == "python":
            cmd = cli_argv(*args)
        else:
            cmd = [self.binary] + self._base_args() + args
        return subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=check,
            env=self._build_env(),
            input=input_text,
            timeout=timeout,
        )

    def _base_args(self) -> list:
        """Return mode-specific database/config args."""
        if self.mode == "rust":
            return ["--database", str(self.db_path)]
        return []

    # -- Commands --

    def init(self) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["init", "-f", str(self.genome_folder)])
        else:
            return self._run(["init", "--genome-folder", str(self.genome_folder)])

    def genome_init(self, name: str, fasta: Path) -> subprocess.CompletedProcess:
        """Initialize a genome from a FASTA file."""
        return self._run(["genome", "init", "-n", name, "--fasta", str(fasta)])

    def alias_set(self, name: str, digest: str) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["alias", "set", "-a", name, "-d", digest])
        else:
            return self._run(["alias", "add", name, digest])

    def alias_get(self, name: str) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["alias", "get", "-a", name])
        else:
            return self._run(["alias", "resolve", name])

    def alias_remove(self, name: str) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["alias", "remove", "-a", name])
        else:
            return self._run(["alias", "remove", name])

    def asset_class_add(self, yaml_path: Path) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["asset-class", "add", "--source", str(yaml_path), "-f"])
        else:
            return self._run(["asset-class", "add", str(yaml_path)])

    def asset_class_list(self) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["asset-class", "list"])
        else:
            return self._run(["list", "asset-classes"])

    def recipe_add(self, yaml_path: Path) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["recipe", "add", "--source", str(yaml_path), "-f"])
        else:
            return self._run(["recipe", "add", str(yaml_path)])

    def recipe_list(self) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["recipe", "list"])
        else:
            return self._run(["list", "recipes"])

    def list_genomes(self) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["genome", "list"])
        else:
            return self._run(["list", "genomes"])

    def list_assets(self, genome: str = None) -> subprocess.CompletedProcess:
        if self.mode == "python":
            args = ["list"]
            if genome:
                args += ["-g", genome]
            return self._run(args)
        else:
            args = ["list", "assets"]
            if genome:
                args += ["--genome", genome]
            return self._run(args)

    def seek(self, path: str, check: bool = False) -> subprocess.CompletedProcess:
        if self.mode == "python":
            args = ["seek", path]
            if check:
                args.append("-e")
            return self._run(args)
        else:
            args = ["seek", path]
            if check:
                args.append("--check")
            return self._run(args)

    def subscribe(self, url: str, reset: bool = False) -> subprocess.CompletedProcess:
        if self.mode == "python":
            args = ["subscribe", "-s", url]
            if reset:
                args.append("-r")
            return self._run(args)
        else:
            args = ["subscribe", url]
            if reset:
                args.append("--reset")
            return self._run(args)

    def unsubscribe(self, url: str) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["unsubscribe", "-s", url])
        else:
            return self._run(["unsubscribe", url])

    def config_get(self) -> subprocess.CompletedProcess:
        """Get configuration (used to verify server subscriptions etc.)."""
        if self.mode == "python":
            return self._run(["config", "get"])
        else:
            return self._run(["config", "get"])

    def build(
        self, registry_path: str, recipe_name: str = None, requirements_only: bool = False
    ) -> subprocess.CompletedProcess:
        args = ["build", registry_path]
        if recipe_name:
            args += ["--recipe-name", recipe_name]
        if requirements_only:
            args += ["-q"]
        return self._run(args)

    def remove(self, registry_path: str) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["remove", registry_path, "-f"])
        else:
            return self._run(["remove", registry_path])

    def tag(self, registry_path: str, new_name: str) -> subprocess.CompletedProcess:
        if self.mode == "python":
            return self._run(["rename", registry_path, "-n", new_name])
        else:
            return self._run(["tag", registry_path, new_name])

    def id(self, registry_path: str) -> subprocess.CompletedProcess:
        return self._run(["id", registry_path])

    def serve_start(self, port: int) -> subprocess.Popen:
        """Start serve in background, returning the Popen handle."""
        if self.mode == "python":
            return popen_cli("serve", "-p", str(port), env=self._build_env(), text=True)
        cmd = [self.binary] + self._base_args() + ["serve", "--port", str(port)]
        return subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=self._build_env(),
        )


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
