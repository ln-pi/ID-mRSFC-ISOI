"""Explicit edge memberships and conserved regional contributions of FC differences."""
from itertools import combinations_with_replacement
import re

import numpy as np


ATLAS_FAMILIES = {"AUD", "LAT", "MED", "PF", "SM", "VIS"}


def edge_layout(roi_labels, edge_index, network_mapping=None):
    """Built-in atlas families are anatomical groupings, not inferred networks."""
    if network_mapping is not None and (
            not isinstance(network_mapping, dict) or not set(roi_labels) <= set(network_mapping)
            or any(not isinstance(network_mapping[r], str) or not network_mapping[r].strip()
                   for r in roi_labels)):
        raise ValueError("descriptive.network_mapping must map every included ROI to a nonempty network")
    families, hemispheres = [], []
    for label in roi_labels:
        match = re.fullmatch(r"([A-Z]+)-.+_([LR])", label)
        families.append(match[1] if match and match[1] in ATLAS_FAMILIES else None)
        hemispheres.append(match[2] if match else None)
    networks = ([network_mapping[r] for r in roi_labels] if network_mapping is not None else families)
    network_kind = "user_defined" if network_mapping is not None else "atlas_anatomical_families"
    groups = []
    a, b = edge_index
    for i, roi in enumerate(roi_labels):
        groups.append(dict(level="roi", name=roi, indices=np.flatnonzero((a == i) | (b == i)),
                           allocation=.5))
    if all(g is not None for g in networks):
        for first, second in combinations_with_replacement(sorted(set(networks)), 2):
            indices = [e for e, (i, j) in enumerate(zip(a, b))
                       if tuple(sorted((networks[i], networks[j]))) == (first, second)]
            if indices:
                groups.append(dict(level="network_pair", name=f"{first} | {second}",
                                   indices=np.asarray(indices), allocation=1.))
    if all(h is not None for h in hemispheres):
        for first, second in combinations_with_replacement(sorted(set(hemispheres)), 2):
            indices = [e for e, (i, j) in enumerate(zip(a, b))
                       if tuple(sorted((hemispheres[i], hemispheres[j]))) == (first, second)]
            if indices:
                groups.append(dict(level="hemisphere_pair", name=f"{first} | {second}",
                                   indices=np.asarray(indices), allocation=1.))
    edges = [dict(edge=e, roi_1=roi_labels[i], roi_2=roi_labels[j],
                  network_1=networks[i], network_2=networks[j],
                  hemisphere_1=hemispheres[i], hemisphere_2=hemispheres[j])
             for e, (i, j) in enumerate(zip(a, b))]
    return groups, edges, dict(network_kind=network_kind,
                              networks_available=all(g is not None for g in networks),
                              roi_network_mapping=dict(zip(roi_labels, networks)))


def regional_summaries(vector, groups, quantity):
    """ROI endpoints get half each; network/hemisphere pairs own unique edges."""
    if quantity not in ("residual_group_loo", "edge_variance"):
        raise ValueError("Unknown topology quantity")
    vector = np.asarray(vector, float)
    if not np.isfinite(vector).all():
        return []
    energy = vector ** 2 if quantity == "residual_group_loo" else vector
    total = energy.sum()
    rows = []
    for group in groups:
        ids, allocation = group["indices"], group["allocation"]
        x = vector[ids]
        allocated = float(energy[ids].sum() * allocation)
        row = dict(level=group["level"], region=group["name"], quantity=quantity,
                   n_edges=len(ids), mean_value=float(x.mean()),
                   allocated_energy=allocated,
                   contribution=allocated / total if total > 0 else np.nan)
        if quantity == "residual_group_loo":
            row.update(rms=float(np.sqrt(np.mean(x ** 2))),
                       positive_fraction=float(np.mean(x > 0)), negative_fraction=float(np.mean(x < 0)),
                       positive_rms=float(np.sqrt(np.mean(np.maximum(x, 0) ** 2))),
                       negative_rms=float(np.sqrt(np.mean(np.minimum(x, 0) ** 2))))
        rows.append(row)
    return rows


def edge_matrix(vector, edge_index, n_rois):
    """No self-connections are measured: retain a masked diagonal."""
    result = np.full((n_rois, n_rois), np.nan)
    a, b = edge_index
    result[a, b] = vector
    result[b, a] = vector
    return result
