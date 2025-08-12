"""Utility functions for running calculations."""

from __future__ import annotations

import contextlib
import os
from pathlib import Path


@contextlib.contextmanager
def chdir(path: Path):
    """
    Change working directory and return to previous on exit.

    Parameters
    ----------
    path
        Path to temporarily change to.
    """
    prev_cwd = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(prev_cwd)
