"""Same-day Run repeatability, distinct from group differences or ICC."""
from itertools import combinations

import numpy as np

from srcs.analysis.connectome import pattern_pearson


def same_day_reliability(store):
    rows = []
    for u, week in enumerate(store.weeks):
        for r1, r2 in combinations(range(len(store.runs)), 2):
            for i, mouse in enumerate(store.mice):
                if not (store.valid_record[i, u, r1] and store.valid_record[i, u, r2]):
                    continue
                a, b = store.F_z[i, u, r1], store.F_z[i, u, r2]
                corr = pattern_pearson(a[None], b[None])[0, 0]
                rows.append(dict(atlas=store.metadata["atlas"], week=week, mouse=mouse,
                    run_1=store.runs[r1], run_2=store.runs[r2], pattern_r=corr,
                    pattern_status="complete" if np.isfinite(corr) else "constant_edge_vector",
                    rms_difference=float(np.sqrt(np.mean((a-b)**2))),
                    mean_signed_difference=float(np.mean(a-b))))
    return rows
