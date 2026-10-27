#!/usr/bin/env python3
"""Wheel/sdist content guard for the packaged web UI bundle.

Dependency-free (stdlib only: zipfile, tarfile, argparse, json), usable both
locally (``./tests/scripts/test-package.sh``) and from any CI job that builds
a distribution artifact.

uv_build ships every file under the module root -- including whatever is
sitting in ``refgenie/server/webui/`` at build time, ``.gitignore``
notwithstanding. That is exactly the behavior release packaging relies on
(CI builds the frontend, then ``uv build`` sweeps it into the wheel), and
exactly the footgun this script exists to catch: a maintainer's stale local
build silently riding into a release, or a build that ran before the frontend
did, shipping no UI at all.

``uv build`` builds the wheel *from* the sdist, so the sdist is checked with
the same rules -- an sdist missing assets is the upstream cause of a wheel
missing them, and the sdist is itself an installable, Node-free artifact.

Usage:
    python tests/scripts/check-wheel-assets.py dist/*.whl dist/*.tar.gz
    python tests/scripts/check-wheel-assets.py dist/*.whl dist/*.tar.gz --require-fresh
    python tests/scripts/check-wheel-assets.py dist/*.whl --max-bytes 8000000
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import tarfile
import zipfile

#: Where the frontend build lands inside the package (refgenie/server/const.py
#: has no equivalent constant importable without the 'dash' extras, so this
#: script -- deliberately dependency-free -- restates the path literally).
WEBUI_PREFIX = "refgenie/server/webui/"
INDEX_PATH = WEBUI_PREFIX + "index.html"
BUILD_INFO_PATH = WEBUI_PREFIX + "build-info.json"
ASSETS_PREFIX = WEBUI_PREFIX + "_app/"

DEFAULT_MAX_BYTES = 8 * 1024 * 1024


#: Contents are only read eagerly for these two members -- everything else in
#: a wheel/sdist is irrelevant to this check and, for a tarball, can be large
#: enough that reading it into memory unconditionally is wasteful.
_READ_CONTENTS_FOR = {INDEX_PATH, BUILD_INFO_PATH}


class ArchiveMember:
    """A (name, size, optional pre-read content) view over one archive entry."""

    def __init__(self, name: str, size: int, data: "bytes | None"):
        self.name = name
        self.size = size
        self._data = data

    def read(self) -> bytes:
        if self._data is None:
            raise ValueError(f"{self.name}: contents were not captured (not in _READ_CONTENTS_FOR)")
        return self._data


def _iter_zip_members(path: str) -> list[ArchiveMember]:
    members = []
    with zipfile.ZipFile(path) as zf:
        for info in zf.infolist():
            if info.is_dir():
                continue
            data = zf.read(info.filename) if info.filename in _READ_CONTENTS_FOR else None
            members.append(ArchiveMember(info.filename, info.file_size, data))
    return members


def _iter_tar_members(path: str) -> list[ArchiveMember]:
    members = []
    with tarfile.open(path) as tf:
        # sdists nest everything under "<name>-<version>/"; strip the first
        # path component so members compare against the same prefixes as a
        # wheel's ("refgenie/server/webui/...", not
        # "refgenie-1.0.0a1/refgenie/server/webui/...").
        for member in tf.getmembers():
            if not member.isfile():
                continue
            parts = member.name.split("/", 1)
            relname = parts[1] if len(parts) == 2 else parts[0]
            data = None
            if relname in _READ_CONTENTS_FOR:
                extracted = tf.extractfile(member)
                data = extracted.read() if extracted is not None else b""
            members.append(ArchiveMember(relname, member.size, data))
    return members


def _iter_members(path: str) -> list[ArchiveMember]:
    if path.endswith(".whl") or path.endswith(".zip"):
        return _iter_zip_members(path)
    if path.endswith(".tar.gz") or path.endswith(".tgz"):
        return _iter_tar_members(path)
    raise ValueError(f"Don't know how to read archive members from: {path}")


def check_archive(path: str, *, require_fresh: bool, max_bytes: int) -> list[str]:
    """Return a list of failure messages for `path`; empty means it passed."""
    errors: list[str] = []
    members = _iter_members(path)
    by_name = {m.name: m for m in members}

    index_member = by_name.get(INDEX_PATH)
    if index_member is None:
        errors.append(f"{path}: missing {INDEX_PATH}")
    elif index_member.size == 0:
        errors.append(f"{path}: {INDEX_PATH} is empty")

    asset_js_members = [
        m for name, m in by_name.items() if name.startswith(ASSETS_PREFIX) and name.endswith(".js")
    ]
    if not asset_js_members:
        errors.append(f"{path}: no {ASSETS_PREFIX}*.js entry found (bundle is empty or unbuilt)")

    if index_member is not None and index_member.size > 0:
        index_text = index_member.read().decode("utf-8", errors="replace")
        if "_app/" not in index_text:
            errors.append(
                f"{path}: {INDEX_PATH} references no path under _app/ "
                "(index emitted against a wiped assets dir?)"
            )

    build_info_member = by_name.get(BUILD_INFO_PATH)
    build_info: dict | None = None
    if build_info_member is None:
        errors.append(f"{path}: missing {BUILD_INFO_PATH}")
    else:
        try:
            build_info = json.loads(build_info_member.read().decode("utf-8"))
        except (UnicodeDecodeError, json.JSONDecodeError) as exc:
            errors.append(f"{path}: {BUILD_INFO_PATH} does not parse as JSON: {exc}")

    if require_fresh:
        expected_sha = os.environ.get("GITHUB_SHA")
        if not expected_sha:
            errors.append("--require-fresh was passed but $GITHUB_SHA is not set")
        elif build_info is not None:
            actual_sha = build_info.get("commit")
            if actual_sha != expected_sha:
                errors.append(
                    f"{path}: {BUILD_INFO_PATH} commit {actual_sha!r} != "
                    f"$GITHUB_SHA {expected_sha!r} (stale build)"
                )
            if build_info.get("dirty", True) is not False:
                errors.append(
                    f"{path}: {BUILD_INFO_PATH} reports dirty={build_info.get('dirty')!r}; "
                    "a release must be built from a clean frontend checkout"
                )

    if max_bytes:
        webui_total = sum(m.size for name, m in by_name.items() if name.startswith(WEBUI_PREFIX))
        if webui_total > max_bytes:
            errors.append(
                f"{path}: {WEBUI_PREFIX} totals {webui_total:,} bytes, over the "
                f"{max_bytes:,}-byte budget (a source map, a font pack or a stray "
                "fixture riding into every `pip install refgenie`?)"
            )

    return errors


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifacts", nargs="+", help="wheel(s) and/or sdist(s) to check")
    parser.add_argument(
        "--require-fresh",
        action="store_true",
        help="also require build-info.json's commit == $GITHUB_SHA and dirty == false",
    )
    parser.add_argument(
        "--max-bytes",
        type=int,
        default=DEFAULT_MAX_BYTES,
        help=f"fail if the uncompressed webui/ payload exceeds this (default {DEFAULT_MAX_BYTES:,})",
    )
    args = parser.parse_args(argv)

    all_errors: list[str] = []
    for artifact in args.artifacts:
        if not os.path.isfile(artifact):
            all_errors.append(f"{artifact}: no such file")
            continue
        errors = check_archive(artifact, require_fresh=args.require_fresh, max_bytes=args.max_bytes)
        if errors:
            all_errors.extend(errors)
        else:
            print(f"OK: {artifact}")

    if all_errors:
        print("\ncheck-wheel-assets: FAILED", file=sys.stderr)
        for error in all_errors:
            print(f"  - {error}", file=sys.stderr)
        return 1

    print("check-wheel-assets: all artifacts carry a fresh, non-empty web UI bundle")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
