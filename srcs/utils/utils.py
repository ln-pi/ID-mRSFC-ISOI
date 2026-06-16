import argparse
import glob
import numpy as np
import os
import pandas as pd
import pickle

from srcs.utils.logger import setup_logger


logger = setup_logger("logs/default.log")


def parse_args():
    parser = argparse.ArgumentParser(description="Initial config inputting...")
    parser.add_argument(
        "-c", "--config",
        default="configs/OIS.yaml",
        help="Path to the config yaml file"
    )
    return parser.parse_args()


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

def extract_roi(matrix, mask_root, atlas_exclude=None, fmt="png", mode="interleaved"):

    if matrix.ndim != 3:
        logger.raise_(ValueError("Input matrix must be 3-D (Time, H, W)"))

    if mask_root is None:
        logger.raise_(IOError("FATAL: no mask_root provided!"))
    
    mask_files = sorted(glob.glob(os.path.join(mask_root, f"*.{fmt}")))

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

def build_half_dataframe(session_paths, trial_cnt=0, group_mode=14):
    s1, s2, s3, s4 = [load_pkl(p) for p in session_paths]

    # check index & columns
    assert s1.index.equals(s2.index) and s1.index.equals(s3.index) and s1.index.equals(s4.index)
    assert s1.columns.equals(s2.columns) and s1.columns.equals(s3.columns) and s1.columns.equals(s4.columns)

    mouse_ids = s1.index
    trials = s1.columns

    if trial_cnt > 0 and trial_cnt <= len(trials):
        trials = trials[:trial_cnt]
    else:
        logger.raise_(ValueError("Trial extracting out of range!"))
    
    df_out = pd.DataFrame(index=["1st", "2nd"], columns=mouse_ids, dtype=object)

    for mouse in mouse_ids:

        mats_1st = []
        mats_2nd = []
        roi_labels=None

        for trial in trials:

            m1 = s1.loc[mouse, trial]
            m2 = s2.loc[mouse, trial]
            m3 = s3.loc[mouse, trial]
            m4 = s4.loc[mouse, trial]

            if roi_labels is None:
                if isinstance(m1, pd.DataFrame):
                    roi_labels = m1.index.tolist()
                else:
                    roi_labels = list(range(extrack_matrix(m1).shape[0]))
            
            if group_mode == 12:
                mats_1st.append(average_matrices(m1, m2))
                mats_2nd.append(average_matrices(m3, m4))
            elif group_mode == 13:
                mats_1st.append(average_matrices(m1, m3))
                mats_2nd.append(average_matrices(m2, m4))
            elif group_mode == 14:
                mats_1st.append(average_matrices(m1, m4))
                mats_2nd.append(average_matrices(m2, m3))
        
        avg_1st = np.mean(np.stack(mats_1st), axis=0)
        avg_2nd = np.mean(np.stack(mats_2nd), axis=0)

        df_out.at["1st", mouse] = pd.DataFrame(avg_1st, index=roi_labels, columns=roi_labels)
        df_out.at["2nd", mouse] = pd.DataFrame(avg_2nd, index=roi_labels, columns=roi_labels)

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
