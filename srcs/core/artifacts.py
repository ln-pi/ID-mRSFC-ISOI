import hashlib
import json
import os
from pathlib import Path
from uuid import uuid4

from .contracts import ProductRef


def sha256(path):
    with Path(path).open("rb") as stream:
        digest = hashlib.sha256()
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
        return digest.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pending = path.with_name(f".{path.name}.{uuid4().hex}.tmp")
    try:
        pending.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True,
                                      indent=2, allow_nan=False), encoding="utf-8")
        os.replace(pending, path)
    finally:
        pending.unlink(missing_ok=True)


def artifact_path(path, anchor):
    """Encode an artifact path relative to the directory containing its manifest."""
    path = Path(path).resolve()
    try:
        return Path(os.path.relpath(path, Path(anchor).resolve())).as_posix()
    except ValueError:  # Different Windows drives/UNC shares cannot be relative.
        return str(path)


def resolve_artifact(value, anchor):
    """Resolve a manifest reference; legacy absolute paths need no rebasing."""
    path = Path(value)
    return (path if path.is_absolute() else Path(anchor) / path).resolve()


def artifact_ref(path, anchor, **fields):
    """Build a checksummed manifest entry with optional identifying fields."""
    return dict(fields, path=artifact_path(path, anchor), sha256=sha256(path))


def artifact_entries(entries, anchor):
    """Yield reference records with absolute Path values from map or list forms."""
    if isinstance(entries, dict):
        entries = (dict(path=path, sha256=checksum) for path, checksum in entries.items())
    for item in entries:
        yield dict(item, path=resolve_artifact(item["path"], anchor))


def publish(path, name, schema, files, *, artifact_root=None):
    path = Path(path).resolve()
    root = Path(artifact_root).resolve() if artifact_root is not None else None
    entries = {}
    for item in files:
        item = Path(item).resolve()
        value = artifact_path(item, path.parent) if root is None or item.is_relative_to(root) else str(item)
        entries[value] = sha256(item)
    write_json(path, {"name": name, "schema": schema, "files": entries})
    return ProductRef(name, schema, path, sha256(path), {"files": len(entries)})
