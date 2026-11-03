"""Validated per-record FC reuse strictly within one configuration."""
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from srcs.utils.dataset import resolve_recording
from srcs.utils.storage import _atomic, save_pickle
import glob
import os
from srcs.utils.logger import setup_logger

logger = setup_logger()


def signature(path):
    path = Path(path)
    from srcs.core.artifacts import sha256
    return [str(path.resolve()), sha256(path)]


from srcs.analysis.connectivity import valid_fc, should_exclude, roi_signals


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
                     "algorithm": ctx.data_identity["stages"]["fc"]["code"]}

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

    return roi_signals(matrix, pairs, atlas_exclude, mode)
