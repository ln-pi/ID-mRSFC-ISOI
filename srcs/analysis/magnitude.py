"""FC amplitude and pattern distances in a single, explicitly selected cohort."""
import numpy as np

from srcs.analysis.connectome import pattern_pearson

SCHEMA = "descriptive_magnitude_topology_v2"


def magnitude(values, valid):
    """Return full mouse axes; undefined statistics remain NaN (never imputed)."""
    values = np.asarray(values, float)
    valid = np.asarray(valid, bool)
    if values.ndim != 2 or valid.shape != values.shape[:1] or values.shape[1] < 3:
        raise ValueError("Expected mouse x edge vectors with at least three edges")
    if not np.isfinite(values[valid]).all():
        raise ValueError("Valid FC vectors must have finite edges")
    n, e = values.shape
    indices = np.flatnonzero(valid)
    count = len(indices)
    distances = {key: np.full((n, n), np.nan) for key in ("euclidean", "rms", "correlation")}
    residual = np.full_like(values, np.nan)
    mean, variance = np.full(e, np.nan), np.full(e, np.nan)
    dispersion = np.nan
    if count:
        x = values[indices]
        mean = x.mean(axis=0)
        # Direct differences avoid cancellation of near-identical FC vectors.
        euclidean = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=-1)
        corr_distance = 1 - pattern_pearson(x, x)
        varied = np.ptp(x, axis=1) > 0
        corr_distance[np.diag_indices(count)] = np.where(varied, 0., np.nan)
        for key, array in (("euclidean", euclidean), ("rms", euclidean / np.sqrt(e)),
                           ("correlation", corr_distance)):
            distances[key][np.ix_(indices, indices)] = array
        if count >= 2:
            variance = x.var(axis=0, ddof=1)
            dispersion = float(np.mean((x - mean) ** 2))
            for k, i in enumerate(indices):
                residual[i] = x[k] - np.delete(x, k, axis=0).mean(axis=0)
    return dict(distances=distances, group_mean=mean, edge_variance=variance,
                residual_group_loo=residual, dispersion=dispersion, n_mice=count)


def select_views(store, config):
    """Run means require all configured Runs, and are computed before distances."""
    weeks = config.get("weeks", store.weeks)
    if (not isinstance(weeks, list) or not weeks or len(set(weeks)) != len(weeks)
            or not set(weeks) <= set(store.weeks)):
        raise ValueError("descriptive.weeks must be unique configured week IDs")
    # Respect configuration order rather than arbitrary request order.
    weeks = [w for w in store.weeks if w in weeks]
    views = config.get("run_views", ["per_run", "mean_fc"])
    if (not isinstance(views, list) or not views or len(set(views)) != len(views)
            or not set(views) <= {"per_run", "mean_fc"}):
        raise ValueError("descriptive.run_views must select per_run and/or mean_fc")
    wi = [store.weeks.index(w) for w in weeks]
    values, masks, labels, definitions = [], [], [], []
    if "per_run" in views:
        for r, run in enumerate(store.runs):
            values.append(store.F_z[:, wi, r, :])
            masks.append(store.valid_record[:, wi, r])
            labels.append(run)
            definitions.append(dict(id=run, strategy="per_run", run_ids=[run]))
    if "mean_fc" in views:
        if "mean_fc" in labels:
            raise ValueError("Run ID mean_fc is reserved when using the averaged FC view")
        valid = store.valid_record[:, wi, :].all(axis=2)
        averaged = store.F_z[:, wi, :, :].mean(axis=2)
        averaged[~valid] = np.nan
        values.append(averaged)
        masks.append(valid)
        labels.append("mean_fc")
        definitions.append(dict(id="mean_fc", strategy="mean_F_z_before_metrics",
                                run_ids=store.runs, missing="all_runs_required"))
    return weeks, labels, np.stack(values, axis=2), np.stack(masks, axis=2), definitions
