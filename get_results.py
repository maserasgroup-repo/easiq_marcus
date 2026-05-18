#!/usr/bin/env python3
"""Compatibility wrapper for the installed `get_results` command."""

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from easiq_marcus.get_results import main


if __name__ == "__main__":
    raise SystemExit(main())
