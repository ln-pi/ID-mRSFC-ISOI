"""Shared descriptive conventions; no group-level resampling is performed here."""
import numpy as np


def individual_identifiability(similarity):
    """Per-mouse identity margin and midrank win fraction in [0, 1].

    Rows and columns must represent the same mice in the same order, with
    independent measurements on the two axes. A tie contributes half a win.
    These are descriptive measures, not independent mouse-level test results.
    """
    s = np.asarray(similarity, dtype=float)
    if s.ndim != 2 or s.shape[0] != s.shape[1] or len(s) < 2:
        raise ValueError('Similarity must be a square matrix with at least two mice.')
    if not np.isfinite(s).all():
        raise ValueError('Similarity must contain only finite values.')
    n = len(s)
    mask = ~np.eye(n, dtype=bool)
    row_others = s[mask].reshape(n, n-1)
    col_others = s.T[mask].reshape(n, n-1)
    own = np.diag(s)
    other = (row_others.mean(axis=1) + col_others.mean(axis=1)) / 2
    def percentile(competitors):
        return ((own[:, None] > competitors) + .5 * (own[:, None] == competitors)).mean(axis=1)
    row, col = percentile(row_others), percentile(col_others)
    return dict(individual_similarity=own.copy(), other_similarity=other,
                self_match_advantage=own-other,
                self_match_percentile=(row+col)/2,
                row_self_match_percentile=row, col_self_match_percentile=col)


def paired_effect_summary(first, second):
    """Signed Cohen's dz = mean(first-second) / sample SD(first-second).

    Pairing is positional and must be aligned by the caller. Nonfinite values
    fail explicitly rather than silently changing the cohort. Undefined dz
    (fewer than two pairs or zero difference variance) is represented by NaN.
    cohen_d is retained as a compatibility alias, explicitly meaning dz.
    """
    a, b = np.asarray(first, dtype=float), np.asarray(second, dtype=float)
    if a.ndim != 1 or a.shape != b.shape:
        raise ValueError('Paired effects require matching one-dimensional arrays.')
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError('Paired effects require finite, aligned observations.')
    diff = a - b
    n = len(diff)
    mean = float(diff.mean()) if n else np.nan
    sd = float(diff.std(ddof=1)) if n >= 2 else np.nan
    status = 'insufficient_pairs' if n < 2 else 'zero_difference_sd' if sd == 0 else 'complete'
    dz = mean / sd if status == 'complete' else np.nan
    return dict(cohen_d=dz, cohen_dz=dz, effect_size_method='paired_cohen_dz',
                effect_size_ddof=1, effect_size_direction='first_minus_second',
                effect_size_status=status, n_pairs=n, mean_difference=mean,
                sd_difference=sd)


def monte_carlo_pvalue_greater(observed, null_values):
    """Inclusive upper-tail randomized p=(b+1)/(B+1), not an exact test."""
    null = np.asarray(null_values, dtype=float)
    if null.ndim != 1 or not len(null) or not np.isfinite(null).all() or not np.isfinite(observed):
        raise ValueError('Monte Carlo p requires a finite observation and nonempty finite 1-D null sample.')
    count = int(np.count_nonzero(null >= observed))
    n = len(null)
    return dict(p_perm=(count + 1) / (n + 1), n_extreme=count, n_iter=n,
                p_perm_min=1 / (n + 1), alternative='greater',
                p_perm_method='monte_carlo_plus_one', ties='greater_or_equal')
