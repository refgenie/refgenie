"""
Run a recipe's container image with Docker or Apptainer.

A recipe's ``docker_image`` names the image its tools live in. Docker runs it
on a laptop; Apptainer runs the same image on a cluster with no Docker daemon.
"""

import os
import shlex
import shutil
from enum import StrEnum
from subprocess import DEVNULL, CalledProcessError, TimeoutExpired, check_output, run

from refgenie.config import config


class ContainerEngine(StrEnum):
    docker = "docker"
    apptainer = "apptainer"


#: Seconds `docker info` may take to answer before Docker counts as unusable.
DOCKER_INFO_TIMEOUT = 20

#: The engine `detect_engine` found, keyed by `REFGENIE_CONTAINER_ENGINE`. Only
#: a found engine is cached, so a daemon that was down is checked again.
_found: dict[str, ContainerEngine] = {}


def _docker_answers() -> bool:
    """Whether a Docker daemon answers here. A `docker` binary alone is not enough."""
    try:
        run(
            ["docker", "info"],
            stdin=DEVNULL,
            stdout=DEVNULL,
            stderr=DEVNULL,
            timeout=DOCKER_INFO_TIMEOUT,
            check=True,
        )
    except (OSError, CalledProcessError, TimeoutExpired):
        return False
    return True


def detect_engine() -> ContainerEngine | None:
    """The engine to run recipe images with, or None when this host has none.

    ``REFGENIE_CONTAINER_ENGINE`` set to ``docker`` or ``apptainer`` picks that
    engine, as long as its binary is on PATH. ``auto`` (the default) takes
    Docker when ``docker info`` answers, else Apptainer when it is on PATH.

    Returns:
        ContainerEngine | None: The engine, or None if none is usable.
    """
    choice = config.container_engine
    if choice in _found:
        return _found[choice]
    if choice == "auto":
        if _docker_answers():
            engine = ContainerEngine.docker
        elif shutil.which("apptainer"):
            engine = ContainerEngine.apptainer
        else:
            return None
    else:
        engine = ContainerEngine(choice)
        if shutil.which(engine) is None:
            return None
    _found[choice] = engine
    return engine


def _image_uri(image: str) -> str:
    """Apptainer needs a transport; a bare image reference is a Docker registry one."""
    return image if "://" in image else f"docker://{image}"


def _env(engine: ContainerEngine) -> dict[str, str] | None:
    """The environment for an engine process: `REFGENIE_APPTAINER_CACHEDIR` passed on.

    Apptainer's default cache is under ``$HOME``, which is small on a cluster.
    """
    if engine is ContainerEngine.apptainer and config.apptainer_cachedir:
        return {**os.environ, "APPTAINER_CACHEDIR": str(config.apptainer_cachedir)}
    return None


def pull_command(engine: ContainerEngine, image: str) -> str:
    """The command that fetches ``image`` ahead of time, for error messages."""
    if engine is ContainerEngine.docker:
        return f"docker pull {image}"
    return f"apptainer pull {_image_uri(image)}"


def one_shot(
    engine: ContainerEngine, image: str, cmd: str, *, mounts: list[str], timeout: float
) -> str:
    """Run ``sh -c cmd`` once inside ``image`` and return its stdout. Used by probes.

    A single ``docker run --rm`` (or ``apptainer exec``): a detached container
    plus ``docker exec -it`` leaks containers and fails without a terminal.

    Args:
        engine: The engine to run the image with.
        image: The image, as a recipe's ``docker_image`` names it.
        cmd: The command, interpreted by ``sh`` inside the image.
        mounts: Host paths to make visible at the same path inside.
        timeout: Seconds before the process is killed.

    Returns:
        str: The command's stdout.

    Raises:
        CalledProcessError: The command, or the engine, exited non-zero.
        TimeoutExpired: The command ran past ``timeout``.
        OSError: The engine binary could not be started.
    """
    mounts = list(dict.fromkeys(mounts))
    if engine is ContainerEngine.docker:
        argv = ["docker", "run", "--rm", "--entrypoint", "sh"]
        for mount in mounts:
            argv += ["-v", f"{mount}:{mount}"]
        argv += [image, "-c", cmd]
    else:
        argv = ["apptainer", "exec", "--cleanenv"]
        for mount in mounts:
            argv += ["--bind", mount]
        argv += [_image_uri(image), "sh", "-c", cmd]
    return check_output(argv, stdin=DEVNULL, timeout=timeout, env=_env(engine)).decode("utf-8")


def wrap(engine: ContainerEngine, image: str, cmd: str, *, mounts: list[str]) -> str:
    """The shell string that runs ``cmd`` inside ``image``, for pypiper to run.

    Each mount is visible at the same path inside, which the relative
    colocation symlinks in an output folder need. The current directory is
    mounted and used as the working directory, as Apptainer does by default,
    and Docker runs as the calling user so build outputs are not owned by root.
    Apptainer gets a clean environment, so host modules cannot leak in.

    Run the result with ``shell=True``: it is one command whose quoted
    argument may hold pipes and redirects.

    Args:
        engine: The engine to run the image with.
        image: The image, as a recipe's ``docker_image`` names it.
        cmd: The recipe command, interpreted by ``sh`` inside the image.
        mounts: Host paths to make visible at the same path inside.

    Returns:
        str: A shell command line.
    """
    mounts = list(dict.fromkeys(mounts))
    if engine is ContainerEngine.docker:
        cwd = os.getcwd()
        argv = ["docker", "run", "--rm", "--user", f"{os.getuid()}:{os.getgid()}"]
        for mount in dict.fromkeys([*mounts, cwd]):
            argv += ["-v", f"{mount}:{mount}"]
        argv += ["--workdir", cwd, "--entrypoint", "sh", image, "-c", cmd]
        return shlex.join(argv)
    argv = ["apptainer", "exec", "--cleanenv"]
    for mount in mounts:
        argv += ["--bind", mount]
    argv += [_image_uri(image), "sh", "-c", cmd]
    prefix = (
        f"APPTAINER_CACHEDIR={shlex.quote(str(config.apptainer_cachedir))} "
        if config.apptainer_cachedir
        else ""
    )
    return prefix + shlex.join(argv)


def image_digest(engine: ContainerEngine, image: str) -> str | None:
    """The registry digest of ``image`` as this host has it, for build provenance.

    Docker reports it with ``docker inspect``. Apptainer has no equivalent for a
    ``docker://`` image (``apptainer inspect`` needs a local image file), so
    under Apptainer this is None.
    """
    if engine is not ContainerEngine.docker:
        return None
    try:
        out = check_output(
            ["docker", "inspect", "--format", "{{index .RepoDigests 0}}", image],
            stdin=DEVNULL,
            stderr=DEVNULL,
            timeout=DOCKER_INFO_TIMEOUT,
        )
    except (OSError, CalledProcessError, TimeoutExpired):
        return None
    return out.decode("utf-8").strip() or None
