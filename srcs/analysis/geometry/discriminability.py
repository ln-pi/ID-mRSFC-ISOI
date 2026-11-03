import numpy as np


def score(row, subject, others):
    own = row[subject]
    others = np.asarray(others, dtype=int)
    others = others[np.isfinite(row[others])]
    if not np.isfinite(own) or not len(others):
        return dict(self_distance=own, q=np.nan, margin=np.nan, nearest_distance=np.nan,
                    nearest_indices=[], ties=0, n_candidates=len(others), status="undefined")
    foreign = row[others]
    nearest = foreign.min()
    return dict(self_distance=own, q=np.mean((own < foreign) + 0.5 * (own == foreign)),
                margin=nearest-own, nearest_distance=nearest,
                nearest_indices=others[foreign == nearest].tolist(), ties=int((own == foreign).sum()),
                n_candidates=len(others), status="ok")
