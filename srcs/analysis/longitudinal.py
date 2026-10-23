"""Animal-block trajectory uncertainty and network-level random slope models."""
import warnings

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SCHEMA = "longitudinal_v1"


def linear_slope(time, values):
    time, values = np.asarray(time, float), np.asarray(values, float)
    keep = np.isfinite(values) & np.isfinite(time)
    if keep.sum() < 3 or np.ptp(time[keep]) <= 0:
        return np.nan
    t = time[keep]-time[keep].mean()
    return float(np.dot(t, values[keep])/np.dot(t, t))


def dispersion_trajectory(values, valid, time, n_bootstrap=200, seed=42):
    """Resample animal trajectories; recompute each week's reference in each draw."""
    values, valid, time = np.asarray(values, float), np.asarray(valid, bool), np.asarray(time, float)
    n, w, _ = values.shape
    if valid.shape != (n, w) or time.shape != (w,) or not np.isfinite(values[valid]).all():
        raise ValueError("Invalid trajectory axes or valid FC vectors")
    if not np.isfinite(time).all() or np.any(np.diff(time) <= 0):
        raise ValueError("Time positions must be finite and strictly increasing")
    def calculate(ids):
        out = np.full(w, np.nan)
        for u in range(w):
            x = values[ids[valid[ids, u]], u]
            if len(x) >= 2:
                out[u] = np.mean((x-x.mean(axis=0))**2)
        return out
    # Do not bootstrap animals absent from the entire selected analysis.
    ids = np.flatnonzero(valid.any(axis=1))
    observed = calculate(ids)
    rng = np.random.default_rng(seed)
    boot = np.full((n_bootstrap, w), np.nan)
    for b in range(n_bootstrap):
        if len(ids):
            boot[b] = calculate(rng.choice(ids, size=len(ids), replace=True))
    ci = np.full((2, w), np.nan)
    for u in range(w):
        finite = boot[:, u][np.isfinite(boot[:, u])]
        if len(finite) >= 20 and valid[:, u].sum() >= 3:
            ci[:, u] = np.quantile(finite, [.025, .975])
    slope = linear_slope(time, observed)
    slopes = np.array([linear_slope(time, y) for y in boot])
    finite = slopes[np.isfinite(slopes)]
    limits = np.quantile(finite, [.025, .975]) if len(finite) >= 20 and len(ids) >= 3 else [np.nan, np.nan]
    direction = ("increasing_dispersion" if limits[0] > 0 else
                 "decreasing_dispersion" if limits[1] < 0 else "inconclusive")
    return dict(dispersion=observed, dispersion_ci=ci, slope=slope, slope_ci=np.array(limits),
                direction=direction, n_valid_mice=valid.sum(axis=0), bootstrap_valid=len(finite))


def rank_persistence(scores, valid):
    w = scores.shape[1]
    rows = []
    for u in range(w):
        for v in range(u+1, w):
            ids = np.flatnonzero(valid[:, u] & valid[:, v] & np.isfinite(scores[:, u]) & np.isfinite(scores[:, v]))
            defined = len(ids) >= 3 and np.ptp(scores[ids, u]) > 0 and np.ptp(scores[ids, v]) > 0
            value = float(spearmanr(scores[ids, u], scores[ids, v]).statistic) if defined else np.nan
            rows.append(dict(u=u, v=v, n_mice=len(ids), spearman_r=value,
                             status="complete" if defined else "insufficient_or_constant_ranks"))
    return rows


def random_slopes(frame, maxiter=300):
    """Gaussian REML: fixed linear time+Run; mouse intercept/slope; mouse-week VC.

    Input columns: mouse, week, run, time, value. Missing observations are not
    imputed. All reported slopes use the original numeric time unit.
    """
    from statsmodels.regression.mixed_linear_model import MixedLM
    frame = frame.copy()
    frame = frame[np.isfinite(frame.value) & np.isfinite(frame.time)]
    summary = dict(n_mice=frame.mouse.nunique(), n_weeks=frame.week.nunique(), n_observations=len(frame),
                   model="value ~ time + Run + (intercept,time|mouse) + (1|mouse:week)",
                   status="unsupported", converged=False)
    visits = frame.groupby("mouse").week.nunique()
    repeated = frame.groupby(["mouse", "week"]).run.nunique()
    if (len(visits) < 4 or (visits < 3).any() or not len(repeated) or (repeated >= 2).sum() < 4
            or np.ptp(frame.value) <= 1e-12):
        summary["reason"] = "need_4_mice_each_3_weeks_and_repeated_Runs_with_nonconstant_values"
        return summary, [], []
    center = float(frame.time.mean())
    scale = float(frame.time.std(ddof=0))
    frame["time_scaled"] = (frame.time-center)/scale
    messages, attempts = [], []
    fitted = None
    for method in ("lbfgs", "powell"):
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                model = MixedLM.from_formula("value ~ time_scaled + C(run)", groups="mouse",
                    re_formula="~time_scaled", vc_formula={"mouse_week": "0+C(week)"}, data=frame)
                candidate = model.fit(reml=True, method=method, maxiter=maxiter, disp=False)
            messages.extend(str(w.message) for w in caught)
            attempts.append(dict(method=method, converged=bool(candidate.converged)))
            fitted = candidate
            if candidate.converged:
                break
        except (ValueError, np.linalg.LinAlgError, FloatingPointError) as exc:
            attempts.append(dict(method=method, error=str(exc)))
    summary.update(time_center=center, time_scale=scale, attempts=attempts, warnings=sorted(set(messages)))
    if fitted is None:
        summary.update(status="failed", reason="all_optimizers_failed")
        return summary, [], []
    cov = np.asarray(fitted.cov_re)
    variance = np.diag(cov)
    boundary = bool(np.min(np.linalg.eigvalsh(cov)) <= 1e-6*max(float(fitted.scale), float(variance.max()), 1e-12))
    hessian_warning = any("Hessian" in s or "singular" in s for s in messages)
    summary.update(converged=bool(fitted.converged), boundary=boundary, hessian_warning=hessian_warning,
        status=("nonconverged" if not fitted.converged else "boundary" if boundary else "converged"),
        fixed_intercept=float(fitted.fe_params["Intercept"]),
        fixed_slope=float(fitted.fe_params["time_scaled"]/scale),
        mouse_intercept_variance=float(cov[0, 0]), mouse_slope_variance=float(cov[1, 1]/scale**2),
        intercept_slope_covariance=float(cov[0, 1]/scale),
        mouse_week_variance=float(fitted.vcomp[0]), error_variance=float(fitted.scale),
        reml_log_likelihood=float(fitted.llf), residual_rmse=np.nan)
    effects, diagnostics = [], []
    if not fitted.converged:
        return summary, effects, diagnostics
    try:
        random = fitted.random_effects
        residuals = np.asarray(fitted.resid)
        predictions = np.asarray(fitted.fittedvalues)
        summary["residual_rmse"] = float(np.sqrt(np.mean(residuals**2)))
        for mouse, effect in random.items():
            effects.append(dict(mouse=mouse, random_intercept=float(effect.iloc[0]),
                random_slope=float(effect["time_scaled"]/scale),
                total_intercept=summary["fixed_intercept"]+float(effect.iloc[0]),
                total_slope=summary["fixed_slope"]+float(effect["time_scaled"]/scale)))
        for k, row in enumerate(frame.itertuples()):
            diagnostics.append(dict(mouse=row.mouse, week=row.week, run=row.run, time=row.time,
                observed=row.value, fitted=predictions[k], residual=residuals[k]))
    except (ValueError, np.linalg.LinAlgError) as exc:
        summary.update(blup_status="unavailable", blup_error=str(exc))
    return summary, effects, diagnostics
