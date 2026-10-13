"""Validated per-record FC reuse strictly within one configuration."""
import hashlib
import inspect
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from srcs.utils.dataset import resolve_recording
from srcs.utils.io import _atomic, save_json, save_pickle
from srcs.preprocessing.signals import extract_roi, should_exclude
from srcs.visualization.plot_series import save_fc_heatmap


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
        # Numerical FC is independent of visualization; report owns all plots.
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

    def plot(self, frame, stem, resume):
        ctx = self.ctx
        options, style = ctx.render_options(), ctx.viz_style("fc_heatmap")
        if not options.get("save", True) and not options.get("show", False):
            return
        folder = Path(ctx.output_dir("Correlation", "Figures"))
        png, marker = folder / (stem + ".png"), folder / (stem + ".plot.json")
        identity = {"values": hashlib.sha256(frame.to_numpy().tobytes()).hexdigest(),
                    "labels": list(frame.index), "options": options, "style": style}
        complete = False
        if resume and png.exists() and marker.exists():
            try:
                complete = json.loads(marker.read_text(encoding="utf-8")) == identity
                with Image.open(png) as img:
                    img.verify()
            except Exception:
                complete = False
        if complete and not options.get("show", False):
            return
        save_fc_heatmap(frame.to_numpy(), names=list(frame.index), name=stem,
                        save_path=str(folder), config=style, **options)
        if options.get("save", True):
            with Image.open(png) as img:
                img.verify()
            save_json(marker, identity)
