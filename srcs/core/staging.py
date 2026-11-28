"""Same-volume staging directories with inherited workspace permissions."""
from contextlib import contextmanager
from pathlib import Path
from shutil import rmtree
from uuid import uuid4


@contextmanager
def staging_directory(parent, label):
    parent = Path(parent)
    parent.mkdir(parents=True, exist_ok=True)
    staging = parent / f".{label}-{uuid4().hex}.tmp"
    staging.mkdir()
    try:
        yield staging
    finally:
        rmtree(staging)
