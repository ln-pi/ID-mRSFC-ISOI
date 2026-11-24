import numpy as np


def distances(query, reference):
    query, reference = np.asarray(query, float), np.asarray(reference, float)
    if query.ndim != 2 or reference.ndim != 2 or query.shape[1] != reference.shape[1]:
        raise ValueError("Expected matching mouse-by-edge matrices")
    rms = np.sqrt(np.mean((query[:, None] - reference[None, :]) ** 2, axis=-1))
    a = query - query.mean(axis=1, keepdims=True)
    b = reference - reference.mean(axis=1, keepdims=True)
    denominator = np.outer(np.linalg.norm(a, axis=1), np.linalg.norm(b, axis=1))
    with np.errstate(invalid="ignore", divide="ignore"):
        correlation = np.clip(a @ b.T / denominator, -1, 1)
    return rms, 1 - correlation


def dispersion(values):
    if len(values) < 2:
        return np.nan
    return np.var(values, axis=0, ddof=1).mean()
