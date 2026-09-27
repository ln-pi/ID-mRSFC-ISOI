"""Small, atomic result writers for interrupted analysis runs."""
import json
import os
import pickle
import tempfile

import numpy as np


def _atomic(path, writer, binary=False):
    folder = os.path.dirname(os.fspath(path)) or "."
    os.makedirs(folder, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(mode="wb" if binary else "w",
        dir=folder, delete=False, **({} if binary else {"encoding": "utf-8"}))
    try:
        with handle:
            writer(handle)
        os.replace(handle.name, path)
    finally:
        if os.path.exists(handle.name):
            os.unlink(handle.name)


def _json_value(value):
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    return value


def save_json(path, value):
    _atomic(path, lambda stream: json.dump(_json_value(value), stream,
        ensure_ascii=False, indent=2, allow_nan=False))


def save_pickle(path, value):
    _atomic(path, lambda stream: pickle.dump(value, stream, protocol=pickle.HIGHEST_PROTOCOL), binary=True)


def save_table(path, frame, index=True):
    _atomic(path, lambda stream: frame.to_csv(stream, index=index))
