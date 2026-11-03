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


def publish(path, name, schema, files):
    path = Path(path)
    entries = {str(Path(item).resolve()): sha256(item) for item in files}
    write_json(path, {"name": name, "schema": schema, "files": entries})
    return ProductRef(name, schema, path, sha256(path),{"files":len(entries)})
