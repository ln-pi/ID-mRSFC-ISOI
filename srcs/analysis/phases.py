"""Stage summaries on in-memory FC patterns and already calculated metrics."""
from itertools import product
import numpy as np
import pandas as pd
from srcs.analysis.connectome import pattern_pearson


def definitions(weeks, config):
    phases = config.get("phases", {})
    names = config.get("summary_phases", [p for p in ("baseline", "early", "middle", "late") if p in phases])
    if len(names) != len(set(names)) or any(p not in phases for p in names):
        raise ValueError("summary_phases must name unique configured phases")
    result = {p: [weeks.index(w) for w in phases[p]] for p in names}
    if any(not indices or len(indices) != len(set(indices)) for indices in result.values()):
        raise ValueError("Phase weeks must be nonempty and unique")
    flattened = [i for indices in result.values() for i in indices]
    if len(flattened) != len(set(flattened)):
        raise ValueError("Summary phases must be disjoint")
    return result


def similarities(store, config):
    phases = definitions(store.weeks, config)
    cases, matrices, counts = [], [], []
    n = len(store.mice)
    window = sorted({i for indices in phases.values() for i in indices})
    for cohort, a, b, first, second in product(("available_phase", "fixed_window"),
            range(len(store.runs)), range(len(store.runs)), phases, phases):
        u, v = phases[first], phases[second]
        keep_a = store.valid_record[:, u, a].all(axis=1)
        keep_b = store.valid_record[:, v, b].all(axis=1)
        if cohort == "fixed_window":
            fixed = store.valid_record[:, window][:, :, sorted({a, b})].all(axis=(1, 2))
            keep_a &= fixed
            keep_b &= fixed
        x = store.F_z[:, u, a].mean(axis=1)
        y = store.F_z[:, v, b].mean(axis=1)
        averaged_fc = pattern_pearson(x, y)
        support = np.zeros((n, n), dtype=int)
        total = np.zeros((n, n))
        for i, j in product(u, v):
            panel = pattern_pearson(store.F_z[:, i, a], store.F_z[:, j, b])
            # Exclude the identical acquisition, preserving valid between-animal comparisons.
            if a == b and i == j:
                np.fill_diagonal(panel, np.nan)
            good = np.isfinite(panel)
            total += np.where(good, panel, 0)
            support += good
        mean_r = np.divide(total, support, out=np.full((n, n), np.nan), where=support > 0)
        mask = keep_a[:, None] & keep_b[None, :]
        for method, matrix in (("correlate_mean_fc", averaged_fc), ("mean_correlations", mean_r)):
            matrix = matrix.copy()
            matrix[~mask] = np.nan
            if method == "correlate_mean_fc" and a == b and first == second:
                np.fill_diagonal(matrix, np.nan)
            cases.append(dict(cohort=cohort, query_run=store.runs[a], reference_run=store.runs[b],
                              phase_u=first, phase_v=second, method=method,
                              weeks_u=[store.weeks[i] for i in u], weeks_v=[store.weeks[i] for i in v]))
            matrices.append(matrix)
            count = support.copy() if method == "mean_correlations" else np.full((n,n),len(u)*len(v))
            count[~np.isfinite(matrix)] = 0
            counts.append(count)
    shape = (0, n, n)
    return cases, np.stack(matrices) if matrices else np.empty(shape), np.stack(counts) if counts else np.empty(shape,dtype=int)


def weekly_metric_means(frame, phases, group_keys, metrics):
    rows = []
    if frame.empty:
        return pd.DataFrame(columns=[*group_keys,"phase","n_requested_weeks","n_valid_weeks",*metrics])
    for labels, part in frame.groupby(group_keys, sort=False, dropna=False):
        if not isinstance(labels, tuple):
            labels = (labels,)
        for phase, weeks in phases.items():
            selected = part[part.week.isin(weeks)]
            row = dict(zip(group_keys, labels), phase=phase, n_requested_weeks=len(weeks))
            for metric in metrics:
                finite = selected.loc[np.isfinite(selected[metric]), metric]
                row[metric] = finite.mean()
                row[metric+"_n_weeks"] = len(finite)
            row["n_valid_weeks"] = int(np.isfinite(selected[metrics[0]]).sum())
            row["complete_phase"] = row["n_valid_weeks"] == len(weeks)
            rows.append(row)
    return pd.DataFrame(rows)
