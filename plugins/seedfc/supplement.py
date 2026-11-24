"""Compatibility entry point for the standalone SeedFC supplement."""
from pathlib import Path
import sys

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from plugins.function.seedfc_supplement import main, run_outputs

if __name__ == "__main__":
    main()
