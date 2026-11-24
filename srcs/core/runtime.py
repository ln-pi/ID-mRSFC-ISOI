from pathlib import Path
from contextlib import contextmanager
import os
from functools import wraps

from .artifacts import sha256
from .identity import digest


def numerical_stage(name, option_keys):
    def bind(action):
        @wraps(action)
        def run(ctx, *, dry_run=False):
            if dry_run:
                return action(ctx,dry_run=True)
            return ctx.execute_stage(name,action)
        run.option_keys = tuple(option_keys)
        return run
    return bind


@contextmanager
def exclusive_run(path):
    """Serialize writers for one output; the OS releases the lock on process exit."""
    path = Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("a+b") as stream:
        stream.seek(0,2)
        if stream.tell()==0:
            stream.write(b"0")
            stream.flush()
        stream.seek(0)
        try:
            if os.name=="nt":
                import msvcrt
                msvcrt.locking(stream.fileno(),msvcrt.LK_NBLCK,1)
            else:
                import fcntl
                fcntl.flock(stream,fcntl.LOCK_EX | fcntl.LOCK_NB)
        except OSError as exc:
            raise RuntimeError(f"Output in use: {path.parent}") from exc
        yield


def cache_key(code, options, inputs):
    return digest({"code": code, "options": options, "inputs": inputs})


def reusable(record, key):
    if record.get("status") != "complete" or record.get("cache", {}).get("key") != key:
        return False
    if not record.get("products"):
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
