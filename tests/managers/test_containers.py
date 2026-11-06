"""How a recipe's container image is run: which engine, and the argv it gets (unit tier).

The same recipe image has to run under Docker on a laptop and under Apptainer
on a cluster with no Docker daemon. Nothing here starts a container: the
subprocess calls are mocked and the argv is what is checked.
"""

import os
import shlex
from pathlib import Path
from subprocess import CalledProcessError
from unittest.mock import patch

import pytest

from refgenie.config import config
from refgenie.managers import containers
from refgenie.managers.containers import ContainerEngine, detect_engine, one_shot, wrap

pytestmark = pytest.mark.unit

IMAGE = "quay.io/biocontainers/bowtie2:2.5.5--ha27dd3b_0"


@pytest.fixture(autouse=True)
def fresh_detection(monkeypatch):
    """Each test sees an empty detection cache and the default settings."""
    monkeypatch.setattr(containers, "_found", {})
    monkeypatch.setattr(config, "container_engine", "auto")
    monkeypatch.setattr(config, "apptainer_cachedir", None)


def host(*, docker_answers: bool, on_path: set[str]):
    """Patch what the host has: a Docker daemon that answers, and binaries on PATH."""

    def which(name):
        return f"/usr/bin/{name}" if name in on_path else None

    def run(argv, **kwargs):
        if not docker_answers:
            raise CalledProcessError(1, argv)

    return (
        patch("refgenie.managers.containers.run", side_effect=run),
        patch("refgenie.managers.containers.shutil.which", side_effect=which),
    )


class TestDetectEngine:
    @pytest.mark.parametrize(
        "docker_answers, on_path, expected",
        [
            pytest.param(True, {"docker", "apptainer"}, ContainerEngine.docker, id="docker-first"),
            pytest.param(
                False, {"docker", "apptainer"}, ContainerEngine.apptainer, id="daemon-down"
            ),
            pytest.param(False, {"apptainer"}, ContainerEngine.apptainer, id="cluster"),
            pytest.param(False, set(), None, id="neither"),
        ],
    )
    def test_auto_prefers_a_docker_that_answers_then_apptainer(
        self, docker_answers, on_path, expected
    ):
        run, which = host(docker_answers=docker_answers, on_path=on_path)
        with run, which:
            assert detect_engine() is expected

    def test_the_setting_picks_apptainer_even_where_docker_answers(self, monkeypatch):
        monkeypatch.setattr(config, "container_engine", "apptainer")
        run, which = host(docker_answers=True, on_path={"docker", "apptainer"})
        with run, which:
            assert detect_engine() is ContainerEngine.apptainer

    def test_a_picked_engine_that_is_not_installed_is_none(self, monkeypatch):
        """So a build preflights `container_engine_missing` instead of failing late."""
        monkeypatch.setattr(config, "container_engine", "apptainer")
        run, which = host(docker_answers=True, on_path={"docker"})
        with run, which:
            assert detect_engine() is None

    def test_a_daemon_started_later_is_found_without_a_restart(self):
        run, which = host(docker_answers=False, on_path={"docker"})
        with run, which:
            assert detect_engine() is None
        run, which = host(docker_answers=True, on_path={"docker"})
        with run, which:
            assert detect_engine() is ContainerEngine.docker


class TestOneShot:
    def test_docker_runs_one_removed_container(self):
        with patch("refgenie.managers.containers.check_output", return_value=b"2.5.5\n") as co:
            out = one_shot(
                ContainerEngine.docker, IMAGE, "bowtie2-build --version", mounts=["/g"], timeout=5
            )
        assert out == "2.5.5\n"
        assert co.call_args.args[0] == [
            "docker", "run", "--rm", "--entrypoint", "sh", "-v", "/g:/g",
            IMAGE, "-c", "bowtie2-build --version",
        ]  # fmt: skip
        assert co.call_args.kwargs["timeout"] == 5

    def test_apptainer_runs_the_same_image_from_its_registry(self):
        with patch("refgenie.managers.containers.check_output", return_value=b"2.5.5\n") as co:
            one_shot(
                ContainerEngine.apptainer,
                IMAGE,
                "bowtie2-build --version",
                mounts=["/g"],
                timeout=5,
            )
        assert co.call_args.args[0] == [
            "apptainer", "exec", "--cleanenv", "--bind", "/g",
            f"docker://{IMAGE}", "sh", "-c", "bowtie2-build --version",
        ]  # fmt: skip
        assert co.call_args.kwargs["env"] is None

    def test_apptainer_gets_the_configured_cache_dir(self, monkeypatch, tmp_path):
        monkeypatch.setattr(config, "apptainer_cachedir", tmp_path)
        with patch("refgenie.managers.containers.check_output", return_value=b"") as co:
            one_shot(ContainerEngine.apptainer, IMAGE, "true", mounts=[], timeout=5)
        assert co.call_args.kwargs["env"]["APPTAINER_CACHEDIR"] == str(tmp_path)


class TestWrap:
    """The shell line a build command runs as. pypiper runs it with shell=True."""

    CMD = "bowtie2-build /g/x.fa /g/out/x | tee log > /g/out/done"

    def test_docker_runs_as_the_caller_in_the_current_folder(self, monkeypatch):
        monkeypatch.chdir("/")
        line = wrap(ContainerEngine.docker, IMAGE, self.CMD, mounts=["/g", "/g/out", "/g"])
        argv = shlex.split(line)
        assert argv[:5] == ["docker", "run", "--rm", "--user", f"{os.getuid()}:{os.getgid()}"]
        assert argv[argv.index("--workdir") + 1] == "/"
        assert [argv[i + 1] for i, a in enumerate(argv) if a == "-v"] == [
            "/g:/g",
            "/g/out:/g/out",
            "/:/",
        ]
        assert argv[-5:] == ["--entrypoint", "sh", IMAGE, "-c", self.CMD]

    def test_apptainer_binds_each_mount_and_cleans_the_environment(self):
        line = wrap(ContainerEngine.apptainer, IMAGE, self.CMD, mounts=["/g", "/g/out"])
        assert shlex.split(line) == [
            "apptainer", "exec", "--cleanenv", "--bind", "/g", "--bind", "/g/out",
            f"docker://{IMAGE}", "sh", "-c", self.CMD,
        ]  # fmt: skip

    def test_apptainer_line_carries_the_configured_cache_dir(self, monkeypatch):
        monkeypatch.setattr(config, "apptainer_cachedir", Path("/scratch/my cache"))
        line = wrap(ContainerEngine.apptainer, IMAGE, "true", mounts=[])
        assert line.startswith("APPTAINER_CACHEDIR='/scratch/my cache' apptainer exec ")
