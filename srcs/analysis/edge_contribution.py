"""P1 Phi/DP from native mouse x week x Run x edge arrays."""
import numpy as np

from srcs.analysis.selection import validate_groups


def edge_contributions(values, split_axis, groups, opposite_ids):
    """Keep opposite-axis samples separate and average raw FC on the split axis."""
    values = np.asarray(values, float)
    if values.ndim != 4 or values.shape[0] < 2 or not np.isfinite(values).all():
        raise ValueError("Edge contribution requires complete (mouse, week, run, edge) data")
    if split_axis not in ("week", "run"):
        raise ValueError("split_axis must be week or run")
    axis = 1 if split_axis == "week" else 2
    validate_groups(groups, values.shape[axis])
    first, second = [np.take(values, group, axis=axis).mean(axis=axis) for group in groups]
    if first.shape[1] != len(opposite_ids):
        raise ValueError("Opposite-axis labels do not match edge array")
    phi, dp = {}, {}
    n, _, edge_count = first.shape
    for k, sample in enumerate(opposite_ids):
        a, b = first[:, k], second[:, k]
        phi[sample] = np.mean(a * b, axis=0)
        total = np.zeros(edge_count)
        for i in range(n):
            own = a[i] * b[i]
            count = np.zeros(edge_count)
            for j in range(n):
                if i != j:
                    count += a[i] * b[j] > own
                    count += a[j] * b[i] > own
            total += -np.log(np.clip(count / (2 * (n-1)), 1e-10, None))
        dp[sample] = total
    return phi, dp


def edge_matrix(vector, edge_index, roi_count):
    result = np.zeros((roi_count, roi_count))
    result[tuple(edge_index)] = vector
    return result + result.T
