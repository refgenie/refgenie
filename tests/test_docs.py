"""Every `refgenie ...` command shown in the docs must parse."""

import re
import shlex
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
DOC_FILES = [
    ROOT / "README.md",
    ROOT / "docs/site/SKILL.md",
    ROOT / "frontend/public/SKILL.md",
    *sorted((ROOT / "docs/refgenie").rglob("*.md")),
]
SHELL_TAGS = {"bash", "console", "sh", "shell"}
FENCE = re.compile(r"^\s*```\s*(\w*)")
OPERATORS = set("|&;()<>") | {"||", "&&", ";;", "|&"}
ASSIGNMENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*=")


def _blocks(path):
    """Yield (first_line_number, [(line_number, text), ...]) per shell block."""
    block, tag, start = None, "", 0
    for n, line in enumerate(path.read_text().splitlines(), 1):
        m = FENCE.match(line)
        if m and block is None:
            block, tag, start = [], m.group(1).lower(), n + 1
        elif m and block is not None:
            if tag in SHELL_TAGS:
                yield start, block
            block = None
        elif block is not None:
            block.append((n, line))


def _join_continuations(lines):
    out, buf, first = [], "", None
    for n, text in lines:
        first = n if first is None else first
        if text.rstrip().endswith("\\"):
            buf += text.rstrip()[:-1] + " "
            continue
        out.append((first, buf + text))
        buf, first = "", None
    if buf:
        out.append((first, buf))
    return out


def _segments(text):
    """Split a shell line into command segments (lists of tokens)."""
    text = re.sub(r"<[A-Za-z_][^<>\s]*>", "PLACEHOLDER", text)
    lex = shlex.shlex(text, posix=True, punctuation_chars=True)
    lex.whitespace_split = True
    lex.commenters = "#"
    try:
        tokens = list(lex)
    except ValueError:
        return []
    segments, cur = [], []
    for tok in tokens:
        if tok in OPERATORS or (set(tok) <= set("|&;()<>") and tok):
            segments.append(cur)
            cur = []
        else:
            cur.append(tok)
    segments.append(cur)
    result = []
    for seg in segments:
        while seg and ASSIGNMENT.match(seg[0]):
            seg = seg[1:]
        if seg and seg[0] == "refgenie":
            result.append(seg)
    return result


def documented_commands():
    """Return (path, line, argv) for each runnable refgenie command."""
    found = []
    for path in DOC_FILES:
        rel = path.relative_to(ROOT).as_posix()
        for _, block in _blocks(path):
            lines = _join_continuations(block)
            if any(t.lstrip().startswith("$ ") for _, t in lines):
                lines = [(n, t.lstrip()[2:]) for n, t in lines if t.lstrip().startswith("$ ")]
            for n, text in lines:
                # `$(refgenie ...)` is its own command.
                texts = re.findall(r"\$\((refgenie[^)]*)\)", text) + [
                    re.sub(r"\$\((refgenie[^)]*)\)", "", text)
                ]
                for t in texts:
                    for seg in _segments(t):
                        found.append((rel, n, seg[1:]))
    return found


COMMANDS = documented_commands()


def test_enough_commands_found():
    assert len(COMMANDS) >= 200
    assert any(
        p == "docs/site/SKILL.md" and argv[:2] == ["alias", "set"] for p, _, argv in COMMANDS
    )


@pytest.mark.parametrize("path,line,argv", COMMANDS, ids=[f"{p}:{n}" for p, n, _ in COMMANDS])
def test_documented_command_parses(path, line, argv, capsys):
    from pydantic import ValidationError

    from refgenie.cli.main import main
    from tests.support.cli import parse

    cmd = "refgenie " + shlex.join(argv)
    if argv == ["--version"]:
        with pytest.raises(SystemExit) as exc:
            main(test_args=argv)
        assert exc.value.code == 0
    elif "--help" in argv or "-h" in argv:
        with pytest.raises(SystemExit) as exc:
            parse(argv)
        assert exc.value.code == 0
    else:
        try:
            parse(argv)
        except (SystemExit, ValidationError) as e:
            err = capsys.readouterr().err.strip().splitlines()
            last = err[-1] if err else str(e)
            pytest.fail(f"{path}:{line}: `{cmd}` does not parse: {last}")


def test_no_local_paths_in_docs():
    files = [ROOT / "frontend/public/SKILL.md"]
    for ext in ("md", "mdx", "yml"):
        files += (ROOT / "docs").rglob(f"*.{ext}")
    bad = [f.relative_to(ROOT).as_posix() for f in files if "file://" in f.read_text()]
    assert not bad, f"file:// links in {bad}"
