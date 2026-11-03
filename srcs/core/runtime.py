from pathlib import Path

from .artifacts import sha256
from .identity import digest


def cache_key(code, options, inputs):
    return digest({"code": code, "options": options, "inputs": inputs})


def reusable(record, key):
    if record.get("status") != "complete" or record.get("cache", {}).get("key") != key:
        return False
    outputs = record["cache"].get("outputs", {})
    return bool(outputs) and all(Path(path).is_file() and sha256(path) == checksum
                                 for path, checksum in outputs.items())


def select_steps(steps, *, stage=None, start=None):
    if stage and start:
        raise ValueError("Select a single step or a starting point")
    first = stage or start or steps[0]
    if first not in steps:
        raise ValueError(f"Unknown step {first}; available: {', '.join(steps)}")
    return (stage,) if stage else tuple(steps[steps.index(first):])
