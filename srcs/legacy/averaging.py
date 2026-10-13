import argparse
import glob
import numpy as np
import os
import pandas as pd
import pickle

from srcs.utils.logger import setup_logger


logger = setup_logger()






# ========== Extract ROIs needed ==========





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
