"""Validate preprocessing outputs without loading complete recordings into RAM."""
import json
import os
from pathlib import Path
import time

import numpy as np
from PIL import Image
import tifffile

from srcs.utils.results import save_json


def tiff_shape(path):
    """Check every stored byte range and decode the first/last page."""
    path = Path(path)
    size = path.stat().st_size
    with tifffile.TiffFile(path) as tif:
        shapes = []
        for page in tif.pages:
            if not page.dataoffsets or any(offset <= 0 or count <= 0 or offset + count > size
                                          for offset, count in zip(page.dataoffsets, page.databytecounts)):
                raise ValueError(f"Incomplete TIFF data: {path}")
            shapes.append(page.shape)
        if not shapes or len(set(shapes)) != 1:
            raise ValueError(f"Inconsistent TIFF pages: {path}")
        tif.pages[0].asarray()
        tif.pages[-1].asarray()
        if len(shapes[0]) == 2:
            return (len(shapes), *shapes[0])
        if len(shapes) == 1 and len(shapes[0]) == 3:
            return tuple(shapes[0])
        raise ValueError(f"Expected grayscale TIFF frames: {path}")


def npy_shape(path):
    data = np.load(path, mmap_mode="r", allow_pickle=False)
    try:
        if data.ndim != 3 or any(n <= 0 for n in data.shape):
            raise ValueError(f"Expected nonempty (time, height, width) array: {path}")
        if data.offset + data.nbytes > Path(path).stat().st_size:
            raise ValueError(f"Incomplete NPY: {path}")
        return tuple(data.shape)
    finally:
        del data


def image_shape(path):
    with Image.open(path) as image:
        image.load()
        return image.height, image.width


def registered_files(folder):
    folder = Path(folder)
    return sorted(p for p in folder.glob("*") if p.suffix.lower() in (".tif", ".tiff"))


def check_tiff(path, expected=None):
    shape = tiff_shape(path)
    return expected is None or shape == tuple(expected)


def check_registration(folder, original):
    paths = registered_files(folder)
    if not paths:
        return False
    source = tiff_shape(original)
    shapes = [tiff_shape(p) for p in paths]
    return all(s[1:] == source[1:] for s in shapes) and sum(s[0] for s in shapes) == source[0]


def check_alignment(array, preview, expected):
    return npy_shape(array) == tuple(expected) and image_shape(preview) == tuple(expected[1:])


def check_filtered(array, affine, spatial_shape):
    original = npy_shape(affine)
    expected = ((original[0] + 5) // 6, *spatial_shape)
    return npy_shape(array) == expected


def signatures(paths):
    result = []
    for path in paths:
        path = Path(path)
        stat = path.stat()
        if not path.is_file() or stat.st_size <= 0:
            raise ValueError(f"Missing or empty file: {path}")
        result.append([str(path.resolve()), stat.st_size, stat.st_mtime_ns])
    return result


class StepRunner:
    def __init__(self, folder, logger, resume=True):
        self.folder = Path(folder)
        self.logger = logger
        self.resume = resume
        self.changed = False

    def run(self, name, action, outputs, validate, dependencies=(), parameters=None):
        """Existing legacy outputs may be adopted only after validation."""
        marker = self.folder / f".{name}.json"
        inputs = signatures(dependencies)
        parameters = parameters or {}
        paths = lambda: list(outputs() if callable(outputs) else outputs)
        state = None
        if marker.exists():
            try:
                state = json.loads(marker.read_text(encoding="utf-8"))
            except (ValueError, OSError):
                state = {"status": "invalid"}
        valid = False
        if self.resume and not self.changed:
            try:
                valid = bool(validate())
                actual = signatures(paths())
                if state is not None:
                    valid = valid and state.get("status") == "complete" and state.get("inputs") == inputs and state.get("parameters") == parameters and state.get("outputs") == actual
                elif inputs and actual:
                    valid = valid and min(p[2] for p in actual) >= max(p[2] for p in inputs)
            except (OSError, ValueError, EOFError, RuntimeError, TypeError, tifffile.TiffFileError):
                valid = False
        if valid:
            if state is None:
                save_json(marker, {"status": "complete", "inputs": inputs,
                    "parameters": parameters, "outputs": signatures(paths()), "adopted": True})
            self.logger.info(f"{name}: valid outputs found; skipping.")
            return False
        self.changed = True
        save_json(marker, {"status": "running", "inputs": inputs, "parameters": parameters})
        self.logger.info(f"{name}: running (missing, incomplete or changed outputs/inputs).")
        action()
        if not validate():
            raise RuntimeError(f"{name} finished without complete valid outputs.")
        save_json(marker, {"status": "complete", "inputs": inputs,
            "parameters": parameters, "outputs": signatures(paths())})
        return True


def preserve_suite2p(folder):
    """Keep an incomplete previous run for inspection instead of mixing chunks."""
    folder = Path(folder).resolve()
    if folder.exists():
        destination = folder.with_name(f"suite2p.previous-{time.time_ns()}")
        if folder.name != "suite2p" or destination.parent != folder.parent:
            raise ValueError("Unexpected Suite2p output directory")
        folder.rename(destination)
