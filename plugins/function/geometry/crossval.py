import numpy as np
from sklearn.covariance import LedoitWolf


def crossvalidated(first, second, mode="identity"):
    first, second = np.asarray(first, float), np.asarray(second, float)
    if first.ndim != 2 or first.shape != second.shape or not first.shape[1]:
        raise ValueError("Expected two aligned mouse-by-edge matrices")
    if mode not in ("identity", "shrinkage", "animal_holdout"):
        raise ValueError("Unknown noise weighting")
    n, e = first.shape
    valid = np.isfinite(first).all(axis=1) & np.isfinite(second).all(axis=1)
    noise = (first - second) / np.sqrt(2)
    out = np.full((n, n), np.nan)
    counts = np.zeros((n, n), dtype=int)
    common = None
    if mode == "shrinkage" and valid.sum() >= 2:
        common = LedoitWolf().fit(noise[valid]).precision_
    for i in np.flatnonzero(valid):
        for j in np.flatnonzero(valid & (np.arange(n) < i)):
            training = valid.copy()
            weight = common
            if mode == "animal_holdout":
                training[[i, j]] = False
                if training.sum() >= 2:
                    weight = LedoitWolf().fit(noise[training]).precision_
            counts[i,j] = counts[j,i] = training.sum() if mode != "identity" else 0
            if weight is None and mode != "identity":
                continue
            delta = first[i] - first[j]
            if weight is not None:
                delta = delta @ weight
            out[i,j] = out[j,i] = delta @ (second[i]-second[j]) / e
    return out, counts
