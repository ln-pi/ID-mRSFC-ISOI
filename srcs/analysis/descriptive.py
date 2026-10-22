"""Run-level amplitude/topology and same-day reliability without side effects."""
from itertools import combinations
import numpy as np
from srcs.analysis.connectome import pattern_pearson


def describe(store, reference_week="03w"):
    reliability, deviations = [], []
    baseline, atlas = reference_week, store.metadata["atlas"]
    z = store.F_z
    reference_index = store.weeks.index(baseline) if baseline in store.weeks else None
    variance = np.full(z.shape[1:], np.nan)
    delta = np.full_like(z, np.nan)
    for u, week in enumerate(store.weeks):
        for r, run_id in enumerate(store.runs):
            indices = np.flatnonzero(store.valid_record[:, u, r])
            if len(indices) >= 2:
                variance[u, r] = z[indices, u, r].var(axis=0, ddof=1)
            for i in indices:
                others = indices[indices != i]
                residual = (z[i, u, r] - z[others, u, r].mean(axis=0)) if len(others) else np.full(z.shape[-1], np.nan)
                has_ref = reference_index is not None and store.valid_record[i, reference_index, r]
                change = z[i, u, r] - z[i, reference_index, r] if has_ref else np.full(z.shape[-1], np.nan)
                delta[i, u, r] = change
                deviations.append(dict(atlas=atlas, mouse=store.mice[i], week=week, run=run_id,
                    n_reference_mice=len(others), reference_mice="|".join(store.mice[j] for j in others),
                    rms_from_group_loo=float(np.sqrt(np.mean(residual ** 2))),
                    mean_signed_group_deviation=float(residual.mean()),
                    reference_week=baseline, reference_status="available" if has_ref else "missing",
                    rms_from_reference=float(np.sqrt(np.mean(change ** 2))),
                    mean_signed_reference_change=float(np.mean(change))))
        for r1, r2 in combinations(range(len(store.runs)), 2):
            for i, mouse in enumerate(store.mice):
                if not (store.valid_record[i, u, r1] and store.valid_record[i, u, r2]):
                    continue
                a, b = z[i, u, r1], z[i, u, r2]
                corr = pattern_pearson(a[None], b[None])[0, 0]
                reliability.append(dict(atlas=atlas, week=week, mouse=mouse,
                    run_1=store.runs[r1], run_2=store.runs[r2], pattern_r=corr,
                    pattern_status="complete" if np.isfinite(corr) else "constant_edge_vector",
                    rms_difference=float(np.sqrt(np.mean((a-b)**2))),
                    mean_signed_difference=float(np.mean(a-b))))
    return dict(reliability=reliability, deviations=deviations,
                edge_variance=variance, delta_from_reference=delta)
