"""Small, atomic result writers for interrupted analysis runs."""
import json
import os
import pickle
import tempfile
import time
import logging

import numpy as np


def _atomic(path, writer, binary=False):
    path = os.path.abspath(os.fspath(path))
    folder = os.path.dirname(path)
    os.makedirs(folder, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(mode="wb" if binary else "w",
        dir=folder, delete=False, **({} if binary else {"encoding": "utf-8"}))
    try:
        with handle:
            writer(handle)
        for attempt in range(6):
            try:
                os.replace(handle.name, path)
                break
            except PermissionError:
                if attempt == 5:
                    raise
                time.sleep(0.05 * (2 ** attempt))
    finally:
        if os.path.exists(handle.name):
            try:
                os.unlink(handle.name)
            except OSError:
                logging.getLogger("project_logger").warning(
                    "Could not remove temporary output file: %s", handle.name)


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
