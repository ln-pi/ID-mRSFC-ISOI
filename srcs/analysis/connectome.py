"""Labeled FC edge store and per-Run pattern similarities.

All axes follow configuration order. No cohort selection or Run averaging occurs
when building the store. NPZ files contain only numeric/string arrays (no pickle).
"""

import numpy as np


EDGE_CLIP = 0.999999  # P1 edge-level transform, distinct from similarity boundaries
SCHEMA = "fc_edges_v1"
SIM_SCHEMA = "pairwise_similarity_v1"


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
