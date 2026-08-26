"""Shared rich console and progress-bar rendering."""

from rich.console import Console
from rich.progress import BarColumn, ProgressColumn, TextColumn, filesize
from rich.text import Text

CONSOLE = Console()


class _DownloadColumn(ProgressColumn):
    """Renders file size downloaded and total, e.g. '0.5/2.3 GB'."""

    @staticmethod
    def render(task):
        """Calculate common unit for completed and total."""
        completed = int(task.completed)
        if task.total is None:
            # Unknown total - just show completed bytes
            unit, suffix = filesize.pick_unit_and_suffix(
                completed,
                ["bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"],
                1024,
            )
            completed_ratio = completed / unit
            precision = 0 if unit == 1 else 1
            completed_str = f"{completed_ratio:,.{precision}f}"
            return Text(f"{completed_str} {suffix}", style="bright_white")

        total = int(task.total)
        unit, suffix = filesize.pick_unit_and_suffix(
            total, ["bytes", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"], 1024
        )
        completed_ratio = completed / unit
        total_ratio = total / unit
        precision = 0 if unit == 1 else 1
        completed_str = f"{completed_ratio:,.{precision}f}"
        total_str = f"{total_ratio:,.{precision}f}"
        download_status = f"{completed_str}/{total_str} {suffix}"
        return Text(download_status, style="bright_white")


class _TransferSpeedColumn(ProgressColumn):
    """Renders human readable transfer speed."""

    @staticmethod
    def render(task):
        """Show data transfer speed."""
        speed = task.speed
        if speed is None:
            return Text("?", style="bright_white")
        data_speed = filesize.decimal(int(speed))
        return Text(f"{data_speed}/s", style="bright_white")


RICH_PROGRESS_COLUMNS = [
    TextColumn("{task.fields[n]}", justify="right"),
    BarColumn(bar_width=None),
    "[magenta]{task.percentage:>3.1f}%",
    "•",
    _DownloadColumn(),
    "•",
    _TransferSpeedColumn(),
]
