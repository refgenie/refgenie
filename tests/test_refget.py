"""Tests for refgenie.refget module (digest functions)."""

from __future__ import annotations

import hashlib

from refgenie.refget import trunc512_digest


class TestTrunc512Digest:
    """Tests for trunc512_digest."""

    def test_deterministic(self):
        """Same input always produces same output."""
        result1 = trunc512_digest("ACGT")
        result2 = trunc512_digest("ACGT")
        assert result1 == result2

    def test_different_sequences_differ(self):
        assert trunc512_digest("ACGT") != trunc512_digest("TGCA")

    def test_default_offset_24_produces_48_hex_chars(self):
        """Default offset=24 means 24 bytes = 48 hex characters."""
        result = trunc512_digest("ACGT")
        assert len(result) == 48

    def test_custom_offset(self):
        """Offset=8 means 8 bytes = 16 hex chars."""
        result = trunc512_digest("ACGT", offset=8)
        assert len(result) == 16

    def test_known_value(self):
        """Verify against manual SHA-512 computation."""
        seq = "ACGT"
        digest = hashlib.sha512(seq.encode()).digest()
        import binascii

        expected = binascii.hexlify(digest[:24]).decode()
        assert trunc512_digest(seq) == expected

    def test_empty_string(self):
        result = trunc512_digest("")
        assert len(result) == 48
        assert isinstance(result, str)
