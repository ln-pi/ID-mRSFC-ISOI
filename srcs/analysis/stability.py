"""Paired-Run repeatability and within-mouse Delta FC persistence."""
from itertools import combinations

import numpy as np

from srcs.analysis.connectome import pattern_pearson
from srcs.analysis.magnitude import magnitude

SCHEMA = "stability_v1"


def run_repeatability(store, weeks):
    rows, audits = [], []
    for week in weeks:
        u = store.weeks.index(week)
        for r, s in combinations(range(len(store.runs)), 2):
            ids = np.flatnonzero(store.valid_record[:, u, r] & store.valid_record[:, u, s])
            mask = np.zeros(len(store.mice), bool)
            mask[ids] = True
            residuals = [magnitude(store.F_z[:, u, k], mask)["residual_group_loo"] for k in (r, s)]
            audits.append(dict(week=week, run_1=store.runs[r], run_2=store.runs[s],
                               template_mice=[store.mice[i] for i in ids]))
            for label, first, second in (("raw_fc", store.F_z[:, u, r], store.F_z[:, u, s]),
                                         ("delta_fc_common_template", *residuals)):
                for i in ids:
                    a, b = first[i], second[i]
                    finite = np.isfinite(a).all() and np.isfinite(b).all()
                    corr = pattern_pearson(a[None], b[None])[0, 0] if finite else np.nan
                    rows.append(dict(mouse=store.mice[i], week=week, run_1=store.runs[r], run_2=store.runs[s],
                        representation=label, n_template_mice=len(ids), pattern_r=corr,
                        pair_mean=float(np.mean((a+b)/2)), mean_signed_difference=float(np.mean(a-b)),
                        rms_difference=float(np.sqrt(np.mean((a-b)**2))),
                        status="complete" if np.isfinite(corr) else "undefined_pattern_or_template"))
    return rows, audits


def residual_persistence(residuals):
    """N x W x V x E -> N x V x W x W; same-record diagonal is not evidence."""
    n, w, v, _ = residuals.shape
    result = np.full((n, v, w, w), np.nan)
    for i in range(n):
        for r in range(v):
            result[i, r] = pattern_pearson(residuals[i, :, r], residuals[i, :, r])
            np.fill_diagonal(result[i, r], np.nan)
    return result
