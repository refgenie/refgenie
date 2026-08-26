"""Tests for refgenie.utils.build.get_dir_digest, which computes the asset
primary key.

The digest must be determined by both file *content* and file *names*
(aligner-index filenames are load-bearing and seek keys resolve to specific
filenames), and must be independent of the builder's locale, or the same content
built on two machines gets two digests.
"""

import hashlib
import locale
import subprocess
import sys

import pytest

from refgenie.utils.build import DIGEST_SCHEME, get_dir_digest


def _make_dir(root, files: dict[str, bytes]):
    """Create a directory tree from a {relative path: content} mapping."""
    for rel, content in files.items():
        path = root / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    return root


# Collation of these three differs between C and en_US.UTF-8: C orders them
# B.txt, _c.txt, a.txt (by byte); en_US.UTF-8 orders them a.txt, B.txt, _c.txt.
_COLLATION_SENSITIVE = {"a.txt": b"one", "B.txt": b"two", "_c.txt": b"three"}


def _digest_under_locale(path, loc: str) -> str:
    """Compute a digest in a subprocess with LC_ALL set to `loc`."""
    result = subprocess.run(
        [
            sys.executable,
            "-c",
            "import sys; from refgenie.utils.build import get_dir_digest; "
            "print(get_dir_digest(sys.argv[1]))",
            str(path),
        ],
        capture_output=True,
        text=True,
        env={"LC_ALL": loc, "LANG": loc, "PATH": "/usr/bin:/bin"},
        check=True,
    )
    return result.stdout.strip()


class TestDirDigest:
    """
    refgenie.utils.build.get_dir_digest, which computes the asset primary key.

    The digest must be determined by both file *content* and file *names*
    (aligner-index filenames are load-bearing and seek keys resolve to specific
    filenames), and must be independent of the builder's locale, or the same
    content built on two machines gets two digests.
    """

    # --- Correctness of the core behavior --------------------------------

    def test_rename_changes_digest(self, tmp_path):
        """
        Identical bytes under a different filename must not be digest-identical.
        A bwa index with correct content under wrong names is broken.
        """
        a = _make_dir(tmp_path / "a", {"index.fa": b"abc"})
        b = _make_dir(tmp_path / "b", {"WRONGNAME.fa": b"abc"})
        assert get_dir_digest(a) != get_dir_digest(b)

    def test_relocation_into_subdirectory_changes_digest(self, tmp_path):
        """Moving a file into a subdirectory must change the digest."""
        a = _make_dir(tmp_path / "a", {"index.fa": b"abc"})
        c = _make_dir(tmp_path / "c", {"sub/index.fa": b"abc"})
        assert get_dir_digest(a) != get_dir_digest(c)

    def test_content_change_changes_digest(self, tmp_path):
        """The original property still holds: content is hashed."""
        a = _make_dir(tmp_path / "a", {"index.fa": b"abc"})
        b = _make_dir(tmp_path / "b", {"index.fa": b"xyz"})
        assert get_dir_digest(a) != get_dir_digest(b)

    def test_identical_layout_and_content_agree(self, tmp_path):
        """Two independently-built directories with the same layout must agree."""
        files = {"index.fa": b"abc", "sub/index.fa.fai": b"def", "n.bt2": b"ghi"}
        a = _make_dir(tmp_path / "a", dict(files))
        b = _make_dir(tmp_path / "b", dict(files))
        assert get_dir_digest(a) == get_dir_digest(b)

    def test_swapping_content_between_filenames_changes_digest(self, tmp_path):
        """
        Binding is per-(name, content) pair, so exchanging content between two
        names changes the digest even though the multiset of bytes is unchanged.
        """
        a = _make_dir(tmp_path / "a", {"one.fa": b"AAA", "two.fa": b"BBB"})
        b = _make_dir(tmp_path / "b", {"one.fa": b"BBB", "two.fa": b"AAA"})
        assert get_dir_digest(a) != get_dir_digest(b)

    # --- Preimage framing ------------------------------------------------

    def test_preimage_framing_is_length_prefixed(self, tmp_path):
        """
        Pin the exact preimage construction: scheme tag, then per entry, a 4-byte
        big-endian path length, the UTF-8 path, and the raw 32-byte content
        digest, ordered byte-wise by encoded path.
        """
        d = _make_dir(tmp_path / "d", {"b.fa": b"BBB", "a.fa": b"AAA"})

        expected = hashlib.sha256()
        expected.update(DIGEST_SCHEME)
        for rel, content in [(b"a.fa", b"AAA"), (b"b.fa", b"BBB")]:
            expected.update(len(rel).to_bytes(4, "big"))
            expected.update(rel)
            expected.update(hashlib.sha256(content).digest())

        assert get_dir_digest(d) == expected.hexdigest()

    # --- Locale independence ---------------------------------------------

    def test_digest_is_locale_independent(self, tmp_path):
        """
        The same directory must digest identically under C and under a UTF-8
        locale. The old shell implementation used `sort`, whose collation is
        locale-dependent, so identical content built on two machines with
        different locale settings produced different digests.
        """
        try:
            locale.setlocale(locale.LC_COLLATE, "en_US.UTF-8")
        except locale.Error:
            pytest.skip("en_US.UTF-8 locale not available on this runner")
        finally:
            locale.setlocale(locale.LC_COLLATE, "C")

        d = _make_dir(tmp_path / "d", dict(_COLLATION_SENSITIVE))
        assert _digest_under_locale(d, "C") == _digest_under_locale(d, "en_US.UTF-8")

    # --- Documented properties -------------------------------------------

    def test_symlinks_do_not_affect_digest(self, tmp_path):
        """
        Colocation symlinks are external references to parent assets, not this
        asset's content, and are stripped from archives. An asset with one and an
        asset without must hash identically.
        """
        external = _make_dir(tmp_path / "external", {"parent.fa": b"parent content"})

        plain = _make_dir(tmp_path / "plain", {"index.fa": b"abc"})
        linked = _make_dir(tmp_path / "linked", {"index.fa": b"abc"})
        (linked / "parent.fa").symlink_to(external / "parent.fa")

        assert get_dir_digest(plain) == get_dir_digest(linked)

    def test_empty_directories_are_invisible(self, tmp_path):
        """
        Only files are hashed, not directory structure, so all content-free
        assets share one digest. Accepted, and pinned here so it stays a
        documented property rather than an accident.
        """
        a = tmp_path / "a"
        a.mkdir()
        b = tmp_path / "b"
        b.mkdir()
        assert get_dir_digest(a) == get_dir_digest(b)

        # An empty subdirectory adds nothing.
        c = _make_dir(tmp_path / "c", {"index.fa": b"abc"})
        (c / "emptysub").mkdir()
        d = _make_dir(tmp_path / "d", {"index.fa": b"abc"})
        assert get_dir_digest(c) == get_dir_digest(d)

    # --- Exclusions ------------------------------------------------------

    def test_macos_metadata_is_always_excluded(self, tmp_path):
        """`._*` files are excluded without the caller asking."""
        a = _make_dir(tmp_path / "a", {"index.fa": b"abc"})
        b = _make_dir(tmp_path / "b", {"index.fa": b"abc", "._index.fa": b"junk"})
        assert get_dir_digest(a) == get_dir_digest(b)

    def test_inherent_globs_match_full_relative_path(self, tmp_path):
        """
        Globs match the normalized relative path, not the basename, so
        directory-scoped patterns work.
        """
        d = _make_dir(
            tmp_path / "d",
            {"index.fa": b"abc", "_refgenie_build/log.txt": b"noise"},
        )
        plain = _make_dir(tmp_path / "plain", {"index.fa": b"abc"})

        assert get_dir_digest(d) != get_dir_digest(plain)
        assert get_dir_digest(d, inherent=["*", "!_refgenie_build/*"]) == get_dir_digest(plain)

    # --- Recipe-declarable inherent set ----------------------------

    def test_default_is_include_all(self, tmp_path):
        """
        Omitting the declaration must hash everything, and must produce exactly
        the digest the hard-coded behavior produced. This pins that adding the
        feature did not silently re-key every asset already in a catalog.
        """
        files = {"index.fa": b"abc", "scratch/tmp.dat": b"junk", "build.log": b"noise"}
        a = _make_dir(tmp_path / "a", files)
        b = _make_dir(tmp_path / "b", files)
        assert get_dir_digest(a) == get_dir_digest(b, inherent=["*"])

    def test_declared_exclusion_drops_incidental_files(self, tmp_path):
        """
        The common shape: everything, minus junk. An asset that differs only in
        its declared-incidental files is the same asset.
        """
        quiet = _make_dir(tmp_path / "quiet", {"index.fa": b"abc"})
        noisy = _make_dir(tmp_path / "noisy", {"index.fa": b"abc", "build.log": b"whatever"})
        inherent = ["*", "!*.log"]

        assert get_dir_digest(quiet, inherent=inherent) == get_dir_digest(noisy, inherent=inherent)
        # ... and without the declaration the log is load-bearing, as it should be.
        assert get_dir_digest(quiet) != get_dir_digest(noisy)

    def test_pure_include_list_needs_no_leading_exclude(self, tmp_path):
        """
        The "only these files" shape, for scratch-heavy assets whose incidental
        outputs cannot be enumerated. A path no pattern matches is not inherent.
        """
        d = _make_dir(
            tmp_path / "d",
            {"index.1.bt2": b"a", "index.2.bt2": b"b", "tmp-4712.scratch": b"noise"},
        )
        plain = _make_dir(tmp_path / "plain", {"index.1.bt2": b"a", "index.2.bt2": b"b"})
        assert get_dir_digest(d, inherent=["index.*"]) == get_dir_digest(plain, inherent=["index.*"])

    def test_last_match_wins(self, tmp_path):
        """
        Ordering decides, gitignore-style — so a later pattern can re-include
        something an earlier one dropped.
        """
        d = _make_dir(tmp_path / "d", {"index.fa": b"abc", "keep.log": b"L", "drop.log": b"D"})

        # drop every log, then take one back
        kept = get_dir_digest(d, inherent=["*", "!*.log", "keep.log"])
        only_drop_gone = _make_dir(tmp_path / "e", {"index.fa": b"abc", "keep.log": b"L"})
        assert kept == get_dir_digest(only_drop_gone)

        # reversing the order drops it again
        assert get_dir_digest(d, inherent=["*", "keep.log", "!*.log"]) != kept

    def test_macos_metadata_cannot_be_re_included(self, tmp_path):
        """
        ALWAYS_EXCLUDED is appended after the recipe's list, so it has the last
        word even when a recipe explicitly asks for those files.
        """
        a = _make_dir(tmp_path / "a", {"index.fa": b"abc"})
        b = _make_dir(tmp_path / "b", {"index.fa": b"abc", "._index.fa": b"junk"})
        inherent = ["*", "._*"]
        assert get_dir_digest(a, inherent=inherent) == get_dir_digest(b, inherent=inherent)

    def test_empty_declaration_hashes_nothing(self, tmp_path):
        """
        An empty list declares nothing inherent, which is distinct from declaring
        nothing at all (None -> include everything). Pinned so the two do not get
        conflated by a falsy check.
        """
        d = _make_dir(tmp_path / "d", {"index.fa": b"abc"})
        empty = tmp_path / "empty"
        empty.mkdir()
        assert get_dir_digest(d, inherent=[]) == get_dir_digest(empty)
        assert get_dir_digest(d, inherent=None) != get_dir_digest(empty)

    # --- Failure behavior ------------------------------------------------

    def test_missing_path_raises(self, tmp_path):
        """
        A digest that cannot be computed must raise, not return None -- the value
        is written straight into the asset primary key.
        """
        with pytest.raises(Exception):
            get_dir_digest(tmp_path / "does_not_exist")

    def test_unreadable_file_raises(self, tmp_path):
        """An unreadable file must raise rather than being silently skipped."""
        if hasattr(sys, "getwindowsversion"):
            pytest.skip("POSIX permission semantics required")
        d = _make_dir(tmp_path / "d", {"index.fa": b"abc"})
        (d / "index.fa").chmod(0o000)
        try:
            if (d / "index.fa").read_bytes():  # root ignores the mode bits
                pytest.skip("running as root; permission bits not enforced")
        except PermissionError:
            pass
        with pytest.raises(PermissionError):
            get_dir_digest(d)
