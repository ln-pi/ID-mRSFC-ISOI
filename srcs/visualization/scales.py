"""Shared color normalization for the finite values actually displayed."""
from dataclasses import dataclass
from math import ceil, floor, log10

from matplotlib import colormaps
import numpy as np

from srcs.utils.visualization import color_policy, settings


def palette(viz, role):
    return colormaps[viz["colormaps"][role]].with_extremes(
        bad=viz["missing_color"], under=viz["under_color"], over=viz["over_color"])


def step(span):
    target = span / 5
    unit = 10 ** floor(log10(target))
    return next(factor * unit for factor in (1, 2, 2.5, 5, 10) if factor * unit >= target)


@dataclass(frozen=True)
class Scale:
    cmap: object
    vmin: float
    vmax: float
    mode: str
    role: str
    observed: tuple
    finite: int
    missing: int
    under: int
    over: int

    @property
    def image(self):
        return dict(cmap=self.cmap, vmin=self.vmin, vmax=self.vmax)

    @property
    def extend(self):
        if self.under and self.over:
            return "both"
        return "min" if self.under else "max" if self.over else "neither"

    def record(self, **selection):
        return dict(**selection, cmap=self.cmap.name, color_mode=self.role,
                    mode=self.mode, center=0 if self.role == "diverging" else None,
                    observed=list(self.observed), limits=[self.vmin, self.vmax],
                    finite=self.finite, missing=self.missing, under=self.under, over=self.over)


def resolve(values, viz, metric, *, anchor_zero=False, bounds=None, fixed=None):
    viz = settings() if viz is None else viz
    role, policy = color_policy(viz, metric)
    raw = np.asarray(values, dtype=float).ravel()
    finite = raw[np.isfinite(raw)]
    if not len(finite):
        return None
    observed = float(finite.min()), float(finite.max())
    low, high = observed
    mode = policy["mode"]
    if fixed is not None:
        mode = "fixed"
        low, high = fixed
    elif mode == "fixed":
        low, high = policy["limits"]
    else:
        if role == "diverging":
            limit = max(abs(low), abs(high))
            if limit:
                tick = step(2 * limit)
                limit = ceil(limit / tick) * tick
            else:
                limit = 1.
            low, high = -limit, limit
        else:
            if anchor_zero:
                low, high = min(0., low), max(0., high)
            if low == high:
                if anchor_zero and low == 0:
                    high = 1.
                else:
                    pad = abs(low) * .05 or .05
                    low, high = low - pad, high + pad
            tick = step(high - low)
            low, high = min(low, floor(low / tick) * tick), max(high, ceil(high / tick) * tick)
        if bounds is not None:
            low, high = max(low, bounds[0]), min(high, bounds[1])
    if not low < high or (role == "diverging" and not low == -high):
        raise ValueError(f"Invalid {role} color limits for {metric}: {(low, high)}")
    return Scale(palette(viz, role), float(low), float(high), mode, role, observed,
                 len(finite), raw.size - len(finite), int((finite < low).sum()), int((finite > high).sum()))
