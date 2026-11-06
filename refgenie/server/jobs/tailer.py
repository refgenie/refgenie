"""Stream a build's pypiper log to the job by polling the log file.

Subprocess output reaches only pypiper's `<name>_log.md`, not our logger. The
file name is not known until the build resolves its asset name, so the tailer
globs for the newest match.
"""

import threading
import time
from contextlib import contextmanager
from pathlib import Path
from typing import TYPE_CHECKING, Iterator

from refgenie.logger import logger
from refgenie.utils.paths import get_build_dir

if TYPE_CHECKING:
    from refgenie.server.jobs.manager import JobContext

__all__ = ["build_log_tailer"]

#: How often to look for new lines.
POLL_SECONDS = 0.5
#: Emit at most this many lines per poll, so a chatty build cannot flood the ring.
MAX_LINES_PER_POLL = 20
#: Truncate any single line to this many characters.
MAX_LINE_CHARS = 2000


@contextmanager
def build_log_tailer(ctx: "JobContext", genome_name: str, asset_group_name: str) -> Iterator[None]:
    """Emit pypiper's log lines as `log` events while the block runs.

    Args:
        ctx: The job context to emit onto.
        genome_name: Genome the build is for.
        asset_group_name: Asset group being built.

    Yields:
        None. On exit the tailer drains once more, so the last lines of a
        failed build are not lost to the poll interval.
    """
    try:
        prefix = get_build_dir(
            genome_folder=Path(ctx.refgenie.genome_folder),
            genome_name=genome_name,
            asset_group_name=asset_group_name,
            # asset_name is not known yet; it may come from the recipe template.
        )
    except Exception:  # noqa: BLE001 - no log tail is not a failed build
        logger.debug("Could not resolve the build directory; build log will not be tailed.")
        yield
        return

    tailer = _Tailer(ctx, prefix)
    tailer.start()
    try:
        yield
    finally:
        tailer.stop()


class _Tailer:
    """A daemon thread that follows the newest `*_log.md` under `prefix`."""

    def __init__(self, ctx: "JobContext", prefix: Path):
        self._ctx = ctx
        self._prefix = prefix
        self._stop = threading.Event()
        self._thread = threading.Thread(
            target=self._loop, name=f"refgenie-tail-{ctx.job_id}", daemon=True
        )
        self._path: Path | None = None
        self._offset = 0
        self._started_at = time.time()

    def start(self) -> None:
        self._thread.start()

    def stop(self) -> None:
        self._stop.set()
        self._thread.join(timeout=POLL_SECONDS * 4)
        # One last read to catch the final lines.
        self._drain()

    def _loop(self) -> None:
        while not self._stop.wait(POLL_SECONDS):
            self._drain()

    def _drain(self) -> None:
        try:
            self._latch()
            if self._path is None:
                return
            with self._path.open("r", errors="replace") as handle:
                handle.seek(self._offset)
                lines = handle.readlines()
                self._offset = handle.tell()
        except OSError:
            return
        except Exception:  # noqa: BLE001 - a log tail must never fail a build
            logger.debug("Build log tail failed", exc_info=True)
            return

        lines = [line.rstrip("\n") for line in lines]
        lines = [line for line in lines if line.strip()]
        if not lines:
            return
        skipped = len(lines) - MAX_LINES_PER_POLL
        if skipped > 0:
            lines = lines[-MAX_LINES_PER_POLL:]
            self._emit(f"... {skipped} log line(s) omitted; full log at {self._path}")
        for line in lines:
            self._emit(line[:MAX_LINE_CHARS])

    def _emit(self, line: str) -> None:
        try:
            self._ctx.log(line, source="pipeline")
        except Exception:  # noqa: BLE001 - see above
            pass

    def _latch(self) -> None:
        """Find this build's pipeline log, once: the newest file touched since the job started."""
        if self._path is not None:
            return
        if not self._prefix.exists():
            return
        candidates = [
            path
            for path in self._prefix.rglob("*_log.md")
            if path.is_file() and path.stat().st_mtime >= self._started_at - 1
        ]
        if not candidates:
            return
        self._path = max(candidates, key=lambda path: path.stat().st_mtime)
        self._offset = 0
        self._ctx.set_log_file(str(self._path))
        self._ctx.phase("run")
