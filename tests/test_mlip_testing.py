"""Test ml_peg import."""

from __future__ import annotations

from ml_peg import __version__


def test_import():
    """Test ml_peg import."""
    assert __version__ is not None
