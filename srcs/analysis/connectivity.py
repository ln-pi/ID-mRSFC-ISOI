"""Validated per-record FC reuse strictly within one configuration."""
import hashlib
import inspect
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from srcs.utils.dataset import resolve_recording
from srcs.utils.results import _atomic, save_pickle
import glob
import os
from srcs.utils.logger import setup_logger

logger = setup_logger()


def signature(path):
    path = Path(path)
    stat = path.stat()
    return [str(path.resolve()), stat.st_size, stat.st_mtime_ns]


def valid_fc(frame, names):
    if not isinstance(frame, pd.DataFrame):
        return False
    if list(frame.index) != names or list(frame.columns) != names:
        return False
    try:
        values = frame.to_numpy(dtype=float)
        return (len(names) >= 2 and np.isfinite(values).all()
                and np.all(np.abs(values) <= 1 + 1e-10)
                and np.allclose(values, values.T, atol=1e-10)
                and np.allclose(np.diag(values), 1))
    except (TypeError, ValueError):
        return False


class CorrelationCache:
    def __init__(self, ctx, atlas):
        self.ctx, self.atlas = ctx, atlas
        self.options = ctx.section("correlation")
        self.masks = sorted(p for p in Path(ctx.mask_root).glob("*.png")
                            if p.stem.casefold() != "chamber_mask")
        if not self.masks:
            raise FileNotFoundError(f"No ROI masks in {ctx.mask_root}")
        shapes = []
        mask_signatures = []
        for path in self.masks:
            with Image.open(path) as img:
                shapes.append(np.asarray(img).shape)
            mask_signatures.append([str(path.resolve()), hashlib.sha256(path.read_bytes()).hexdigest()])
        if any(len(s) != 2 or s != shapes[0] for s in shapes):
            raise ValueError("ROI masks must be 2-D and have matching dimensions")
        self.shape = shapes[0]
        names = ([p.stem + "_L" for p in self.masks]
                 + [p.stem + "_R" for p in reversed(self.masks)])
        self.names = [n for n in names if not should_exclude(n, atlas)]
        self.base = {"atlas": atlas, "masks": mask_signatures, "mode": "symmetric",
                     "algorithm": hashlib.sha256((inspect.getsource(extract_roi)
                         + inspect.getsource(should_exclude)
                         + Path(__file__).read_text(encoding="utf-8")).encode()).hexdigest()}

    def process(self, week, mouse, run):
        ctx = self.ctx
        source = resolve_recording(ctx, week, mouse, run)
        key = dict(self.base, config_tag=ctx.analysis_tag, source=signature(source))  # Missing source must never be hidden by a cache.
        stem = f"{week}_{mouse}_{run}_FC_Matrix_{'_'.join(self.atlas)}_excluded"
        folder = Path(ctx.correlation_path) / "FC_Matrix"
        folder.mkdir(parents=True, exist_ok=True)
        cached = folder / (stem + ".cache.pkl")
        numeric = folder / (stem + ".npy")
        frame = None
        reused = False
        resume = self.options.get("resume", True)
        if resume and cached.exists():
            try:
                candidate = pd.read_pickle(cached)
                if valid_fc(candidate, self.names) and candidate.attrs.get("correlation_signature") == key:
                    frame = candidate
            except Exception as exc:
                ctx.logger.warning("Invalid FC cache %s; recomputing: %s", cached, exc)
        if frame is not None:
            reused = True
            ctx.logger.info("Reusing valid FC: %s", stem)
        else:
            matrix = np.load(source)
            series, names = extract_roi(matrix, mask_root=ctx.mask_root,
                                        atlas_exclude=self.atlas, mode="symmetric")
            frame = pd.DataFrame(np.corrcoef(series.T), index=names, columns=names)
            if not valid_fc(frame, self.names):
                raise ValueError("FC contains non-finite values or insufficient/invalid ROI signals")
        frame.attrs["correlation_signature"] = key
        # Cache is committed before plotting, so a plot failure does not lose computed FC.
        if not reused or not cached.exists():
            save_pickle(cached, frame)
        values = frame.to_numpy()
        numeric_ok = False
        if resume and numeric.exists():
            try:
                numeric_ok = np.array_equal(np.load(numeric, allow_pickle=False), values)
            except Exception:
                pass
        if not numeric_ok:
            _atomic(numeric, lambda f: np.save(f, values, allow_pickle=False), binary=True)
        return frame, "reused" if reused else "success"



def should_exclude(name, exclude_list):
    if exclude_list is None:
        return False

    for exclude in exclude_list:
        if isinstance(exclude, dict):
            for side, area in exclude.items():
                if area in name:
                    if side == "LR":
                        return True
                    elif side == "L" and name.endswith("_L"):
                        return True
                    elif side == "R" and name.endswith("_R"):
                        return True
        elif isinstance(exclude, str):
            if "_LR" in exclude:
                area = exclude.replace("_LR", "")
                if area in name:
                    return True
            elif "_L" in exclude:
                area = exclude.replace("_L", "")
                if area in name and name.endswith("_L"):
                    return True
            elif "_R" in exclude:
                area = exclude.replace("_R", "")
                if area in name and name.endswith("_R"):
                    return True
            else:
                if exclude in name:
                    return True

    return False

def extract_roi(matrix, mask_root, atlas_exclude=None, fmt="png", mode="interleaved"):

    if matrix.ndim != 3:
        logger.raise_(ValueError("Input matrix must be 3-D (Time, H, W)"))

    if mask_root is None:
        logger.raise_(IOError("FATAL: no mask_root provided!"))

    mask_files = [
        p for p in sorted(glob.glob(os.path.join(mask_root, f"*.{fmt}")))
        if os.path.splitext(os.path.basename(p))[0].casefold() != "chamber_mask"
    ]
    if not mask_files:
        raise FileNotFoundError(f"No ROI masks found in {mask_root} (excluding Chamber_Mask).")

    pairs = []

    for file in mask_files:
        name = os.path.basename(file).split(".")[0]

        if fmt == "png":
            from PIL import Image
            arr = np.array(Image.open(file))
        elif fmt == "npy":
            arr = np.load(file)
        else:
            logger.raise_(ValueError(f"Unsupported format: {fmt}"))

        if arr.shape != matrix.shape[1:]:
            raise ValueError(f"ROI mask {file} shape {arr.shape} does not match data shape {matrix.shape[1:]}.")
        mask_L = arr > 0
        mask_R = np.fliplr(mask_L)

        pairs.append((name, mask_L, mask_R))

    if mode == "interleaved":
        items = [(f"{n}_L", mL) for n, mL, _ in pairs] + [(f"{n}_R", mR) for n, _, mR in pairs]
        items = [x for pair in zip(items[:len(pairs)], items[len(pairs):]) for x in pair]
    elif mode == "symmetric":
        items = [(f"{n}_L", mL) for n, mL, _ in pairs] + [(f"{n}_R", mR) for n, _, mR in reversed(pairs)]
    else:
        logger.raise_(ValueError(f"Unknown mode: {mode}"))

    if atlas_exclude is not None:
        items = [(n, m) for n, m in items if not should_exclude(n, atlas_exclude)]

    roi_ts = []
    names = []

    for name, mask in items:
        time_series = matrix[:, mask]

        valid_voxel = ~np.isnan(time_series).all(axis=0)
        time_series = time_series[:, valid_voxel]

        if time_series.shape[1] == 0:
            roi_mean = np.full(matrix.shape[0], np.nan)
        else:
            roi_mean = np.nanmean(time_series, axis=1)

        roi_ts.append(roi_mean)
        names.append(name)

    roi_time_series = np.stack(roi_ts, axis=1)

    return roi_time_series, names
