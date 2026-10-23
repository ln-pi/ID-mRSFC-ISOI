"""Balanced categorical-week/Run adjusted nested REML, fitted to all edges.

Model: y_itr = week_t + run_r + u_i + v_it + error_itr.
Orthogonal restricted contrasts have variances error, error+R*interaction,
and error+R*interaction+R*W*mouse. Ordered pooling is the constrained REML
optimum, including zero-variance boundaries; it is not moment clipping.
"""
import numpy as np


def balanced_reml(y):
    y = np.asarray(y, float)
    if y.ndim != 4 or not np.isfinite(y).all():
        raise ValueError("Balanced REML requires finite mouse x week x Run x edge data")
    n, w, r, e = y.shape
    if n < 3 or w < 2 or r < 2:
        raise ValueError("Balanced REML requires >=3 mice, >=2 weeks, >=2 Runs")
    grand = y.mean(axis=(0, 1, 2))
    cell = y.mean(axis=2)
    mouse = cell.mean(axis=1)
    week = cell.mean(axis=0)
    run = y.mean(axis=(0, 1))
    df = np.array([(n*w-1)*(r-1), (n-1)*(w-1), n-1], float)
    ss = np.stack([
        np.sum((y-cell[:, :, None, :]-run[None, None, :, :]+grand)**2, axis=(0, 1, 2)),
        r*np.sum((cell-mouse[:, None, :]-week[None, :, :]+grand)**2, axis=(0, 1)),
        w*r*np.sum((mouse-grand)**2, axis=0)])
    ms = ss / df[:, None]
    eigenvalues = np.empty_like(ms)
    # Only three strata. Weighted pool-adjacent-violators enforces 0<=l0<=l1<=l2.
    for j in range(e):
        pools = []
        for k in range(3):
            pools.append(([k], df[k], ss[k, j]))
            while len(pools) > 1 and pools[-2][2]/pools[-2][1] > pools[-1][2]/pools[-1][1]:
                right, left = pools.pop(), pools.pop()
                pools.append((left[0]+right[0], left[1]+right[1], left[2]+right[2]))
        for ids, weight, total in pools:
            eigenvalues[ids, j] = total/weight
    error = eigenvalues[0]
    interaction = (eigenvalues[1]-error)/r
    stable = (eigenvalues[2]-eigenvalues[1])/(w*r)
    total = stable+interaction+error
    divide = lambda a, b: np.divide(a, b, out=np.full_like(a, np.nan), where=b > 0)
    # BLUPs are conditional predictions, not independent observed FC components.
    shrink = divide(w*r*stable, eigenvalues[2])
    u = (mouse-grand)*np.nan_to_num(shrink)
    v_shrink = divide(interaction, interaction+error/r)
    v = (cell-week[None, :, :]-u[:, None, :])*np.nan_to_num(v_shrink)
    return dict(mouse_variance=stable, mouse_week_variance=interaction, error_variance=error,
        stable_icc=divide(stable, total), same_week_icc=divide(stable+interaction, total),
        mean_run_stable_icc=divide(stable, stable+interaction+error/r),
        stable_mouse_blup=u, mouse_week_blup=v, week_mean=week, run_effect=run-grand,
        mean_squares=ms, restricted_eigenvalues=eigenvalues, degrees_of_freedom=df,
        boundary=(stable == 0) | (interaction == 0) | (error == 0),
        status=np.where(total <= 0, "constant", np.where((stable == 0) | (interaction == 0) | (error == 0),
                                                       "boundary", "interior")))


def variance_intervals(y, n_bootstrap=200, seed=42):
    """Resample complete animals; retain each animal's weeks and Runs together."""
    rng = np.random.default_rng(seed)
    names = ("mouse_variance", "mouse_week_variance", "error_variance", "stable_icc")
    samples = {key: [] for key in names}
    for _ in range(n_bootstrap):
        fit = balanced_reml(y[rng.integers(len(y), size=len(y))])
        for key in names:
            samples[key].append(fit[key])
    output = {}
    for key, values in samples.items():
        a = np.asarray(values)
        ci = np.full((2, y.shape[-1]), np.nan)
        for j in range(y.shape[-1]):
            finite = a[:, j][np.isfinite(a[:, j])] if len(a) else np.array([])
            if len(finite) >= 20:
                ci[:, j] = np.quantile(finite, [.025, .975])
        output[key+"_ci"] = ci
        output[key+"_bootstrap_valid"] = np.isfinite(a).sum(axis=0) if len(a) else np.zeros(y.shape[-1], int)
    return output
