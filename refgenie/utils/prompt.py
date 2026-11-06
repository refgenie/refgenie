"""
Confirmation prompts that never block a non-interactive caller.

Library code must not read stdin, since a script or server has no terminal.
Operations that want confirmation take ``confirm: Callable[[str], bool] | None``.
With none given, `resolve_confirmer` refuses. The CLI calls
`enable_interactive_prompts()` once at startup to make the default a real prompt;
nothing else should call it.
"""

from collections.abc import Callable

from refgenie.logger import logger

Confirmer = Callable[[str], bool]

_interactive = False


def enable_interactive_prompts() -> None:
    """Make the default confirmer prompt on the terminal. Called by the CLI only."""
    global _interactive
    _interactive = True


def deny(message: str) -> bool:
    """Refuse without asking. The default for a caller with no terminal."""
    logger.info(f"Not confirmed (non-interactive): {message}")
    return False


def ask(message: str, default: bool = True) -> bool:
    """Prompt on the terminal; refuse cleanly if there is no input to read."""
    from rich.prompt import Confirm

    try:
        return Confirm.ask(message, default=default)
    except EOFError:
        logger.info(f"No input available (non-interactive); declining: {message}")
        return False


def resolve_confirmer(confirm: Confirmer | None, default: bool = True) -> Confirmer:
    """
    Pick the confirmer to use for one prompt.

    Args:
        confirm: The caller's confirmer, if it supplied one.
        default: The answer to pre-select when prompting interactively.

    Returns:
        Confirmer: The caller's confirmer, an interactive prompt if the CLI
            enabled them, or `deny`.
    """
    if confirm is not None:
        return confirm
    if _interactive:
        return lambda message: ask(message, default=default)
    return deny
