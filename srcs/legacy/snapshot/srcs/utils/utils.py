import argparse
import glob
import numpy as np
import os
import pandas as pd
import pickle

from srcs.utils.logger import setup_logger


logger = setup_logger()


def parse_args():
    parser = argparse.ArgumentParser(description="Initial config inputting...")
    parser.add_argument(
        "-c", "--config",
        default="configs/OIS.yaml",
        help="Path to the config yaml file"
    )
    return parser.parse_args()


def normalize_data(data, mask=None, axis=-1, eps=1e-8, mask_fill=np.nan):
    """
    Z-score normalize data along a given axis.

    For OIS stacks arranged as H x W x T, the default normalizes each pixel's
    time series. If mask is provided, pixels outside mask are filled with
    mask_fill, which defaults to NaN.
    """
    data = np.asarray(data, dtype=np.float64)

    if mask is None:
        mean = np.nanmean(data, axis=axis, keepdims=True)
        std = np.nanstd(data, axis=axis, keepdims=True)
        normalized = (data - mean) / (std + eps)
        return np.nan_to_num(normalized, nan=0.0, posinf=0.0, neginf=0.0)

    mask = np.asarray(mask, dtype=bool)
    if mask.shape != data.shape[:mask.ndim]:
        logger.raise_(ValueError(
            f"Mask shape {mask.shape} is incompatible with data shape {data.shape}"
        ))

    normalized = np.full(data.shape, mask_fill, dtype=np.float64)
    valid_data = data[mask]
    mean = np.nanmean(valid_data, axis=axis, keepdims=True)
    std = np.nanstd(valid_data, axis=axis, keepdims=True)
    valid_normalized = (valid_data - mean) / (std + eps)
    normalized[mask] = np.nan_to_num(valid_normalized, nan=0.0, posinf=0.0, neginf=0.0)

    return normalized


# ========== Extract ROIs needed ==========
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

def load_chamber_mask(mask_root, expected_shape=None):
    """Load the shared chamber mask, allowing filename case differences."""
    from PIL import Image
    files = [
        p for p in glob.glob(os.path.join(mask_root, "*"))
        if os.path.basename(p).casefold() == "chamber_mask.png" and os.path.isfile(p)
    ]
    if len(files) != 1:
        raise FileNotFoundError(f"Expected one Chamber_Mask.png in {mask_root}; found {len(files)}.")
    with Image.open(files[0]) as image:
        mask = np.asarray(image.convert("L")) > 0
        if "A" in image.getbands():
            mask &= np.asarray(image.getchannel("A")) > 0
    if expected_shape is not None and mask.shape != tuple(expected_shape):
        raise ValueError(f"Chamber mask shape {mask.shape} does not match data shape {tuple(expected_shape)}.")
    if not mask.any():
        raise ValueError("Chamber mask contains no valid pixels.")
    return mask


def extract_roi(matrix, mask_root, atlas_exclude=None, fmt="png", mode="interleaved"):

    if matrix.ndim != 3:
        logger.raise_(ValueError("Input matrix must be 3-D (Time, H, W)"))

    if mask_root is None:
        logger.raise_(IOError("FATAL: no mask_root provided!"))
    
    mask_files = [
        p for p in sorted(glob.glob(os.path.join(mask_root, f"*.{fmt}")))
        if os.path.splitext(os.path.basename(p))[0].casefold() != "chamber_mask"
    ]
    if not mask_files:
        raise FileNotFoundError(f"No ROI masks found in {mask_root} (excluding Chamber_Mask).")

    pairs = []

    for file in mask_files:
        name = os.path.basename(file).split(".")[0]

        if fmt == "png":
            from PIL import Image
            arr = np.array(Image.open(file))
        elif fmt == "npy":
            arr = np.load(file)
        else:
            logger.raise_(ValueError(f"Unsupported format: {fmt}"))
        
        if arr.shape != matrix.shape[1:]:
            raise ValueError(f"ROI mask {file} shape {arr.shape} does not match data shape {matrix.shape[1:]}.")
        mask_L = arr > 0
        mask_R = np.fliplr(mask_L)

        pairs.append((name, mask_L, mask_R))
    
    if mode == "interleaved":
        items = [(f"{n}_L", mL) for n, mL, _ in pairs] + [(f"{n}_R", mR) for n, _, mR in pairs]
        items = [x for pair in zip(items[:len(pairs)], items[len(pairs):]) for x in pair]
    elif mode == "symmetric":
        items = [(f"{n}_L", mL) for n, mL, _ in pairs] + [(f"{n}_R", mR) for n, _, mR in reversed(pairs)]
    else:
        logger.raise_(ValueError(f"Unknown mode: {mode}"))
    
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


# ========== Average over matrices ==========
def load_pkl(path):
    with open(path, "rb") as f:
        return pickle.load(f)

def extrack_matrix(cell):
    if isinstance(cell, pd.DataFrame):
        return cell.values
    elif isinstance(cell, np.ndarray):
        return cell
    else:
        logger.raise_(TypeError(f"Unsupport cell type: {type(cell)}"))

def fisher_z(r):
    r = np.clip(r, -0.999999, 0.999999)
    return np.arctanh(r)

def average_matrices(m1, m2):
    z1 = fisher_z(extrack_matrix(m1))
    z2 = fisher_z(extrack_matrix(m2))
    return (z1 + z2) / 2

def split_group_indices(session_count, group_mode="sequential", repeat_session=None, seed=None):
    """Return zero-based groups; repeat_session is a one-based position."""
    if isinstance(session_count, bool) or not isinstance(session_count, (int, np.integer)) or session_count < 2:
        raise ValueError("At least two sessions are required.")
    if group_mode not in ("sequential", "random", 12, 13, 14):
        raise ValueError("group_mode must be sequential, random, or legacy 12/13/14.")
    if repeat_session is not None:
        if group_mode != "sequential" or session_count % 2 == 0:
            raise ValueError("repeat_session requires sequential mode and an odd session count.")
        if isinstance(repeat_session, bool) or not isinstance(repeat_session, (int, np.integer)) or not 1 <= repeat_session <= session_count:
            raise ValueError("repeat_session must be a one-based session position.")
    if group_mode in (12, 13, 14):
        if session_count != 4:
            raise ValueError("Legacy 12/13/14 require four sessions; use sequential or random.")
        return {12: ([0, 1], [2, 3]), 13: ([0, 2], [1, 3]), 14: ([0, 3], [1, 2])}[group_mode]
    order = list(range(session_count))
    if group_mode == "random":
        order = np.random.default_rng(seed).permutation(session_count).tolist()
    half = session_count // 2
    if session_count % 2 == 0:
        return order[:half], order[half:]
    shared = order[half] if repeat_session is None else int(repeat_session) - 1
    remaining = [i for i in order if i != shared]
    first, second = remaining[:half] + [shared], [shared] + remaining[half:]
    if group_mode == "sequential":
        first.sort()
        second.sort()
    return first, second


split_session_indices = split_group_indices  # Compatibility for external callers.


def build_half_dataframe(session_paths, trial_cnt=0, group_mode="sequential",
                         repeat_session=None, seed=None, mouse_ids=None, trial_ids=None):
    """Average Fisher-z matrices over equal session halves and selected trials.

    Input order defines session order. trial_cnt=0 uses all trials.
    repeat_session is a one-based position for odd sequential splits.
    Outputs remain in Fisher-z space, matching the original implementation.
    """
    session_paths = list(session_paths)
    groups = split_group_indices(len(session_paths), group_mode, repeat_session, seed)
    sessions = [p if isinstance(p, pd.DataFrame) else load_pkl(p) for p in session_paths]
    if not all(isinstance(frame, pd.DataFrame) for frame in sessions):
        raise ValueError("Each session must contain a mouse-by-trial DataFrame.")
    if mouse_ids is not None or trial_ids is not None:
        sessions = [
            frame.loc[mouse_ids if mouse_ids is not None else frame.index,
                      trial_ids if trial_ids is not None else frame.columns]
            for frame in sessions
        ]
    reference = sessions[0]
    for frame in sessions:
        if not frame.index.is_unique or not frame.columns.is_unique:
            raise ValueError("Mouse and trial labels must be unique.")
        if not reference.index.equals(frame.index) or not reference.columns.equals(frame.columns):
            raise ValueError("All sessions must have identical mouse and trial labels in the same order.")
    if reference.empty:
        raise ValueError("Sessions must contain at least one mouse and one trial.")
    if isinstance(trial_cnt, bool) or not isinstance(trial_cnt, (int, np.integer)) or not 0 <= trial_cnt <= len(reference.columns):
        raise ValueError("trial_cnt must be between 0 (all trials) and the available trial count.")
    trials = reference.columns[:trial_cnt] if trial_cnt else reference.columns
    df_out = pd.DataFrame(index=["1st", "2nd"], columns=reference.index, dtype=object)
    df_out.attrs["session_groups"] = {
        name: [sessions[i].attrs.get("axis_id", str(i)) if isinstance(session_paths[i], pd.DataFrame)
               else str(session_paths[i]) for i in indices]
        for name, indices in zip(df_out.index, groups)
    }
    for mouse in reference.index:
        first_cell = reference.loc[mouse, trials[0]]
        shape = extrack_matrix(first_cell).shape
        if len(shape) != 2 or shape[0] != shape[1]:
            raise ValueError("Connectivity matrices must be square.")
        labels = first_cell.index if isinstance(first_cell, pd.DataFrame) else pd.RangeIndex(shape[0])
        for name, indices in zip(df_out.index, groups):
            matrices = []
            for i in indices:
                for trial in trials:
                    cell = sessions[i].loc[mouse, trial]
                    matrix = extrack_matrix(cell)
                    if matrix.shape != shape:
                        raise ValueError("Connectivity matrices must have matching shapes.")
                    if isinstance(cell, pd.DataFrame) and (not cell.index.equals(labels) or not cell.columns.equals(labels)):
                        raise ValueError("Connectivity matrices must have matching ROI labels and order.")
                    matrices.append(fisher_z(matrix))
            df_out.at[name, mouse] = pd.DataFrame(np.mean(np.stack(matrices), axis=0), index=labels, columns=labels)
    return df_out


# ========== DataFrame -> Array ==========
def df_to_array(df):

    ind = df.index.tolist()
    cols = df.columns.tolist()

    n_ind = len(ind)
    n_col = len(cols)
    
    n_reg = 0
    labels = None

    for r in ind:
        for c in cols:
            cell = df.loc[r, c]

            if isinstance(cell, pd.DataFrame):
                labels = cell.index.tolist()
                mat_1st = cell.values
            elif isinstance(cell, np.ndarray):
                labels = list(range(cell.shape[0]))
                mat_1st = cell
            else:
                continue

            n_reg = mat_1st.shape[0]
            break

        if n_reg is not None:
            break
    
    if n_reg is None:
        logger.raise_(ValueError("No valid data connectivity matrices found."))
    
    tensor = np.full(
        (n_ind, n_col, n_reg, n_reg),
        np.nan,
        dtype=float
    )

    for i, r in enumerate(ind):
        for j, c in enumerate(cols):
            cell = df.loc[r, c]

            if isinstance(cell, pd.DataFrame):
                tensor[i, j] = cell.values
            elif isinstance(cell, np.ndarray):
                tensor[i, j] = cell
    
    return tensor, ind, cols, labels
