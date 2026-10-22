"""Labeled FC edge store and per-Run pattern similarities.

All axes follow configuration order. No cohort selection or Run averaging occurs
when building the store. NPZ files contain only numeric/string arrays (no pickle).
"""
from dataclasses import dataclass
import json

import numpy as np
import pandas as pd

from srcs.utils.results import _atomic

EDGE_CLIP = 0.999999  # P1 edge-level transform, distinct from similarity boundaries
SCHEMA = "fc_edges_v1"
SIM_SCHEMA = "pairwise_similarity_v1"


def save_arrays(path, metadata, **arrays):
    _atomic(path, lambda stream: np.savez_compressed(
        stream, metadata=np.asarray(json.dumps(metadata, ensure_ascii=False)), **arrays), binary=True)


def load_arrays(path, schema):
    with np.load(path, allow_pickle=False) as data:
        metadata = json.loads(str(data["metadata"]))
        if metadata.get("schema") != schema:
            raise ValueError(f"Unsupported schema: {path}")
        arrays = {k: data[k] for k in data.files if k != "metadata"}
    return metadata, arrays


def pattern_pearson(first, second):
    """Pairwise Pearson over complete edge vectors; undefined patterns stay NaN."""
    a, b = np.asarray(first, float), np.asarray(second, float)
    if a.ndim != 2 or b.ndim != 2 or a.shape[1] != b.shape[1] or a.shape[1] < 3:
        raise ValueError("Pattern arrays require matching edge axes with at least three edges")
    def normalize(x):
        finite = np.isfinite(x).all(axis=1)
        # A rounded mean of a constant vector can leave a nonzero residual.
        # Test variation before centering so it cannot become a spurious pattern.
        varied = np.ptp(np.where(finite[:, None], x, 0), axis=1) > 0
        centered = np.where(finite[:, None], x, 0)
        centered = centered - centered.mean(axis=1, keepdims=True)
        scale = np.max(np.abs(centered), axis=1, keepdims=True)
        centered = np.divide(centered, scale, out=np.zeros_like(centered), where=scale != 0)
        length = np.linalg.norm(centered, axis=1, keepdims=True)
        valid = finite & varied & (length[:, 0] > 0)
        return np.divide(centered, length, out=np.zeros_like(centered), where=length != 0), valid
    x, valid_a = normalize(a)
    y, valid_b = normalize(b)
    result = np.clip(x @ y.T, -1, 1)
    result[~(valid_a[:, None] & valid_b[None, :])] = np.nan
    return result


def similarity_z(values):
    """Numerical finite boundary only; no display cap. NaNs are preserved."""
    bound = np.nextafter(1.0, 0.0)
    return np.arctanh(np.clip(values, -bound, bound))


@dataclass
class EdgeStore:
    metadata: dict
    F_r: np.ndarray
    F_z: np.ndarray
    valid_record: np.ndarray
    invalid_reason: np.ndarray
    edge_index: np.ndarray

    @property
    def mice(self):
        return self.metadata["mouse_ids"]

    @property
    def weeks(self):
        return self.metadata["week_ids"]

    @property
    def runs(self):
        return self.metadata["run_ids"]

    @property
    def roi_labels(self):
        return self.metadata["roi_labels"]

    @classmethod
    def from_tables(cls, tables, mice, weeks, runs, atlas, record_audit=None):
        # Fix the ROI/edge basis using the first valid labeled FC in config order.
        from srcs.analysis.connectivity import valid_fc
        labels = None
        for week in weeks:
            table = tables.get(week)
            if table is None:
                continue
            if (not isinstance(table, pd.DataFrame) or list(table.index) != list(mice)
                    or list(table.columns) != list(runs) or not table.index.is_unique
                    or not table.columns.is_unique):
                raise ValueError(f"FC table axes differ from configuration: {week}")
            for mouse in mice:
                for run in runs:
                    cell = table.at[mouse, run]
                    if (isinstance(cell, pd.DataFrame) and cell.index.is_unique
                            and len(cell) >= 3 and valid_fc(cell, list(cell.index))):
                        if labels is None:
                            labels = list(cell.index)
        if labels is None:
            raise ValueError(f"No valid labeled FC matrices for {atlas}")
        edge_index = np.asarray(np.triu_indices(len(labels), k=1))
        shape = (len(mice), len(weeks), len(runs))
        raw = np.full((*shape, edge_index.shape[1]), np.nan)
        valid = np.zeros(shape, bool)
        reasons = np.full(shape, "missing_record", dtype="U256")
        for u, week in enumerate(weeks):
            table = tables.get(week)
            if table is None:
                reasons[:, u, :] = "missing_week_table"
                continue
            for i, mouse in enumerate(mice):
                for r, run in enumerate(runs):
                    cell = table.at[mouse, run]
                    if record_audit is not None:
                        row = record_audit[(week, mouse, run, atlas)]
                        if row["status"] not in ("success", "reused"):
                            reasons[i, u, r] = f"audit_{row['status']}: {row.get('reason', '')}"
                            continue
                    if cell is None or (pd.api.types.is_scalar(cell) and pd.isna(cell)):
                        continue
                    if not isinstance(cell, pd.DataFrame) or not valid_fc(cell, labels):
                        reasons[i, u, r] = "invalid_fc_or_roi_order"
                        continue
                    raw[i, u, r] = cell.to_numpy(dtype=float)[tuple(edge_index)]
                    valid[i, u, r], reasons[i, u, r] = True, ""
        meta = dict(schema=SCHEMA, atlas=atlas, axes=["mouse", "week", "run", "edge"],
                    shape=list(raw.shape), mouse_ids=list(mice), week_ids=list(weeks),
                    run_ids=list(runs), roi_labels=labels, edge_order="numpy.triu_indices(k=1)",
                    edge_transform=dict(method="atanh", clip=[-EDGE_CLIP, EDGE_CLIP]),
                    missing_policy="retain_configured_axes_nan_no_imputation",
                    algorithm="fc-edges-1")
        return cls(meta, raw, np.arctanh(np.clip(raw, -EDGE_CLIP, EDGE_CLIP)),
                   valid, reasons, edge_index)

    def save(self, path):
        save_arrays(path, self.metadata, F_r=self.F_r, F_z=self.F_z,
                    valid_record=self.valid_record, invalid_reason=self.invalid_reason,
                    edge_index=self.edge_index)

    @classmethod
    def load(cls, path):
        meta, arrays = load_arrays(path, SCHEMA)
        obj = cls(meta, **arrays)
        shape = (len(obj.mice), len(obj.weeks), len(obj.runs), len(obj.edge_index[0]))
        if (obj.F_r.shape != shape or obj.F_z.shape != shape
                or obj.valid_record.shape != shape[:3] or obj.invalid_reason.shape != shape[:3]
                or not np.array_equal(obj.edge_index, np.asarray(np.triu_indices(len(obj.roi_labels), 1)))):
            raise ValueError("Edge store dimensions/labels are inconsistent")
        return obj

    def subset(self, mice, weeks, runs, scale="z"):
        indices = [self.mice.index(m) for m in mice], [self.weeks.index(w) for w in weeks], [self.runs.index(r) for r in runs]
        if scale not in ("r", "z"):
            raise ValueError("Edge scale must be r or z")
        values = self.F_z if scale == "z" else self.F_r
        return values[np.ix_(*indices, range(values.shape[-1]))]

    def average_halves(self, mice, weeks, runs, split_axis, groups):
        """Two mouse x edge mean-F_z arrays; never reconstruct nested matrices."""
        from srcs.analysis.selection import validate_groups
        if split_axis not in ("week", "run"):
            raise ValueError("split_axis must be week or run")
        values = self.subset(mice, weeks, runs)
        if not np.isfinite(values).all():
            raise ValueError("Combined selection requires complete FC records")
        axis = 1 if split_axis == "week" else 2
        validate_groups(groups, values.shape[axis])
        return tuple(np.take(values, group, axis=axis).mean(axis=(1, 2)) for group in groups)


def pairwise_tensor(store, run_id):
    """S[i,j,u,v]; same-record r=1 is retained but excluded from statistics/z."""
    r = store.runs.index(run_id)
    n, w, _, e = store.F_z.shape
    vectors = store.F_z[:, :, r].reshape(n * w, e)
    flat = pattern_pearson(vectors, vectors)
    self_indices = np.arange(n * w)
    defined = np.isfinite(np.diag(flat))
    flat[self_indices[defined], self_indices[defined]] = 1.0
    self_mask = np.eye(n * w, dtype=bool).reshape(n, w, n, w).transpose(0, 2, 1, 3)
    raw = flat.reshape(n, w, n, w).transpose(0, 2, 1, 3)
    valid = np.isfinite(raw)
    stat_mask = valid & ~self_mask
    z = similarity_z(raw)
    z[~stat_mask] = np.nan
    reasons = np.zeros(raw.shape, np.uint8)
    reasons[~valid] = 1
    reasons[self_mask & valid] = 2
    pattern_valid = np.isfinite(np.diag(flat)).reshape(n, w)
    pattern_reason = np.full((n, w), "", dtype="U256")
    pattern_reason[~store.valid_record[:, :, r]] = store.invalid_reason[:, :, r][~store.valid_record[:, :, r]]
    pattern_reason[store.valid_record[:, :, r] & ~pattern_valid] = "constant_edge_vector"
    meta = dict(schema=SIM_SCHEMA, atlas=store.metadata["atlas"], run_id=run_id,
                axes=["mouse_i", "mouse_j", "week_u", "week_v"], shape=list(raw.shape),
                mouse_ids=store.mice, week_ids=store.weeks,
                input_scale="F_z", similarity="Pearson_over_edges",
                similarity_z_boundary="clip to +/- nextafter(1,0), no display cap; same-record z=NaN",
                invalid_reason_codes={"0": "defined", "1": "undefined_endpoint_see_pattern_reason",
                                      "2": "same_record_excluded_from_statistics"},
                algorithm="pairwise-1")
    return meta, dict(S_r=raw, S_z=z, valid_similarity=valid, statistical_mask=stat_mask,
                      invalid_reason=reasons, pattern_reason=pattern_reason)
