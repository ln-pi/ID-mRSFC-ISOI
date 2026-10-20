"""Compatibility entry point; prefer python run.py --stage statistics."""
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from srcs.utils.config import load_context
from srcs.workflow.statistics import run

if __name__ == "__main__":
    run(load_context())
