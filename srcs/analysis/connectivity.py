"""FC validation and extraction on in-memory arrays."""
import numpy as np
import pandas as pd


def valid_fc(frame, names):
    if not isinstance(frame, pd.DataFrame):
        return False
    if list(frame.index) != names or list(frame.columns) != names:
        return False
    try:
        values = frame.to_numpy(dtype=float)
        return (len(names) >= 2 and np.isfinite(values).all()
                and np.all(np.abs(values) <= 1 + 1e-10)
                and np.allclose(values, values.T, atol=1e-10)
                and np.allclose(np.diag(values), 1))
    except (TypeError, ValueError):
        return False


def should_exclude(name, exclude_list):
    if exclude_list is None:
        return False

    for exclude in exclude_list:
        if isinstance(exclude, dict):
            for side, area in exclude.items():
                if area in name:
                    if side == "LR":
                        return True
                    elif side == "L" and name.endswith("_L"):
                        return True
                    elif side == "R" and name.endswith("_R"):
                        return True
        elif isinstance(exclude, str):
            if "_LR" in exclude:
                area = exclude.replace("_LR", "")
                if area in name:
                    return True
            elif "_L" in exclude:
                area = exclude.replace("_L", "")
                if area in name and name.endswith("_L"):
                    return True
            elif "_R" in exclude:
                area = exclude.replace("_R", "")
                if area in name and name.endswith("_R"):
                    return True
            else:
                if exclude in name:
                    return True

    return False

def roi_signals(matrix, pairs, atlas_exclude, mode):
    if mode == "interleaved":
        items = [(f"{n}_L", mL) for n, mL, _ in pairs] + [(f"{n}_R", mR) for n, _, mR in pairs]
        items = [x for pair in zip(items[:len(pairs)], items[len(pairs):]) for x in pair]
    elif mode == "symmetric":
        items = [(f"{n}_L", mL) for n, mL, _ in pairs] + [(f"{n}_R", mR) for n, _, mR in reversed(pairs)]
    else:
        raise ValueError(f"Unknown mode: {mode}")

    if atlas_exclude is not None:
        items = [(n, m) for n, m in items if not should_exclude(n, atlas_exclude)]

    roi_ts = []
    names = []

    for name, mask in items:
        time_series = matrix[:, mask]

        valid_voxel = ~np.isnan(time_series).all(axis=0)
        time_series = time_series[:, valid_voxel]

        if time_series.shape[1] == 0:
            roi_mean = np.full(matrix.shape[0], np.nan)
        else:
            roi_mean = np.nanmean(time_series, axis=1)

        roi_ts.append(roi_mean)
        names.append(name)

    roi_time_series = np.stack(roi_ts, axis=1)

    return roi_time_series, names
