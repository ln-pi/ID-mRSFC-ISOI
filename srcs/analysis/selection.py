"""Deterministic disjoint observation splits."""
import numpy as np


def validate_groups(groups, count):
    if (len(groups) != 2 or not all(len(g) for g in groups)
            or sorted(list(groups[0]) + list(groups[1])) != list(range(count))):
        raise ValueError("Halves must partition the selected axis without shared observations")


def split_group_indices(count, group_mode="sequential", seed=None):
    """Disjoint halves; odd counts use floor/ceil sizes without shared weeks."""
    if isinstance(count, bool) or not isinstance(count, (int, np.integer)) or count < 2:
        raise ValueError("At least two observations are required")
    if group_mode not in ("sequential", "random"):
        raise ValueError("group_mode must be sequential or random; numeric legacy modes were removed")
    order = list(range(count))
    if group_mode == "random":
        order = np.random.default_rng(seed).permutation(count).tolist()
    half = count // 2
    return order[:half], order[half:]
