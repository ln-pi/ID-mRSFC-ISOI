"""Pure paired fingerprint and exploratory edge calculations."""
import numpy as np
from srcs.analysis.statistics import paired_effect_summary,individual_identifiability


def edge_vector(matrix, fisher=True):
    a = np.asarray(matrix, dtype=float)
    if a.ndim != 2 or a.shape[0] != a.shape[1] or a.shape[0] < 3:
        raise ValueError('At least three aligned ROIs are required')
    v = a[np.triu_indices(len(a), 1)]
    if not np.isfinite(a).all() or not np.allclose(a, a.T):
        raise ValueError('Matrix must be finite and symmetric')
    if np.any(np.abs(v) > 1 + 1e-10):
        raise ValueError('FC values outside [-1, 1]')
    return np.arctanh(np.clip(v, -0.999999, 0.999999)) if fisher else v


def pattern_similarity(first, second):
    a, b = np.asarray(first), np.asarray(second)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Aligned finite edge vectors required')
    if np.ptp(a) == 0 or np.ptp(b) == 0:
        raise ValueError('Constant edge vectors have undefined similarity')
    return float(np.corrcoef(a, b)[0, 1])


def accuracy(similarity):
    s = np.asarray(similarity)
    n = len(s)
    if s.shape != (n, n) or n < 2:
        raise ValueError('At least two paired mice required')
    other = ~np.eye(n, dtype=bool)
    row = np.diag(s) > np.max(np.where(other, s, -np.inf), axis=1)
    col = np.diag(s) > np.max(np.where(other, s, -np.inf), axis=0)
    return dict(accuracy=float(np.mean(row & col)), row_accuracy=float(row.mean()),
                column_accuracy=float(col.mean()))


def identify(first, second, mice, n_iter=0, seed=42):
    """Inputs are mouse x Fisher-z edge vectors; no second edge transform."""
    a, b = np.asarray(first, float), np.asarray(second, float)
    if a.shape != b.shape or a.ndim != 2 or len(a) != len(mice) or len(mice) < 2:
        raise ValueError('Matching arrays and at least two paired mice required')
    raw = np.array([[pattern_similarity(x, y) for y in b] for x in a])
    sim = np.arctanh(np.clip(raw, -np.tanh(4), np.tanh(4)))
    own = np.diag(sim)
    other = np.array([np.mean(np.r_[np.delete(sim[i], i), np.delete(sim[:, i], i)])
                      for i in range(len(mice))])
    effect = paired_effect_summary(own, other)
    out = dict(mice=list(mice), raw_similarity=raw, similarity_z=sim, self_z=own,
               other_z=other, advantage_z=own-other, effect=effect, **accuracy(sim))
    out['individual_metrics']=individual_identifiability(sim)
    out['relative_effect']=1-np.divide(other,own,out=np.full_like(own,np.nan),where=own!=0)
    if isinstance(n_iter, bool) or not isinstance(n_iter, int) or n_iter < 0:
        raise ValueError('n_iter must be a nonnegative integer')
    if n_iter:
        rng = np.random.default_rng(seed)
        null = np.array([accuracy(sim[:, rng.permutation(len(mice))])['accuracy']
                         for _ in range(n_iter)])
        out.update(permutation_accuracy=null,
                   p_permutation=float((1 + np.sum(null >= out['accuracy']))/(n_iter+1)))
    return out


def edge_contribution(first, second):
    """Original Phi/DP convention on raw FC edges, fixed explicit pairs."""
    a, b = np.asarray(first, float), np.asarray(second, float)
    if a.shape != b.shape or a.ndim != 2 or len(a) < 2 or not np.isfinite([a,b]).all():
        raise ValueError('At least two finite paired mice required for Edge')
    phi = np.mean(a*b, axis=0)
    dp = np.zeros(a.shape[1])
    for i in range(len(a)):
        own = a[i]*b[i]
        count = sum((a[i]*b[j] > own).astype(int) + (a[j]*b[i] > own).astype(int)
                    for j in range(len(a)) if j != i)
        dp -= np.log(np.clip(count/(2*(len(a)-1)), 1e-10, None))
    return dict(phi=phi, dp=dp, effective_repeat=1, interpretation='exploratory identification contribution')
