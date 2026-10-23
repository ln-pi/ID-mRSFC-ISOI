"""Code and input identities used by numerical products."""
import hashlib
from importlib.metadata import version
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[2]


def sha256(path):
    hasher = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            hasher.update(block)
    return hasher.hexdigest()


def code_identity():
    def git(*args):
        result = subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True)
        return result.stdout.strip() if result.returncode == 0 else "unavailable"
    files = [ROOT / "run.py", *sorted((ROOT / "srcs").rglob("*.py"))]
    hashes = {p.relative_to(ROOT).as_posix(): sha256(p) for p in files}
    return dict(commit=git("rev-parse", "HEAD"), working_tree=git("status", "--porcelain"),
                source_sha256=hashes, packages={p: version(p) for p in ("numpy", "pandas", "scipy", "statsmodels")})
