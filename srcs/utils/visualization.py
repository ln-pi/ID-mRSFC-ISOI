"""Resolve project visualization settings, including for saved historical outputs."""
from collections.abc import Mapping
from copy import deepcopy
from functools import lru_cache
from math import isfinite
from pathlib import Path

import yaml


@lru_cache(maxsize=1)
def defaults():
    path = Path(__file__).resolve().parents[2] / "configs/default.yaml"
    return yaml.safe_load(path.read_text(encoding="utf-8"))["visualization"]


def merge(base, overrides):
    result = deepcopy(base)
    for key, value in overrides.items():
        if isinstance(value, Mapping) and isinstance(result.get(key), Mapping):
            result[key] = merge(result[key], value)
            if key == "color_scale" and value.get("mode") == "auto" and "limits" not in value:
                result[key].pop("limits", None)
        else:
            result[key] = deepcopy(value)
    return result


def settings(config=None):
    return merge(defaults(), (config or {}).get("visualization", {}))


def metric_options(viz, metric):
    return viz["metrics"].get(metric, {})


def color_policy(viz, metric, role=None):
    options = metric_options(viz, metric)
    role = role or options["color_mode"]
    override = options.get("color_scale", {})
    scale = merge(viz["color_scale"], override)
    if override.get("mode") == "auto" and "limits" not in override:
        scale.pop("limits", None)
    return role, scale


def validate(viz):
    if viz["time_start"] not in ("start", "baseline"):
        raise ValueError("visualization.time_start must be start or baseline")
    maps = viz["colormaps"]
    if set(maps) != {"sequential", "diverging"} or any(
            not isinstance(value, str) or not value for value in maps.values()):
        raise ValueError("visualization.colormaps requires sequential and diverging names")
    from matplotlib import colormaps
    from matplotlib.colors import is_color_like

    for name in (*maps.values(), viz["graph_palette"], viz["alignment"]["atlas_palette"]):
        if name not in colormaps:
            raise ValueError(f"Unknown Matplotlib color map: {name}")
    for key in ("missing_color", "under_color", "over_color"):
        if not is_color_like(viz[key]):
            raise ValueError(f"Invalid visualization.{key}")
    if viz["color_scale"]["scope"] != "batch" or viz["color_scale"]["mode"] not in ("auto", "fixed"):
        raise ValueError("Invalid visualization.color_scale scope or mode")
    if viz["color_scale"]["mode"] == "auto" and "limits" in viz["color_scale"]:
        raise ValueError("Auto visualization.color_scale cannot specify limits")
    if len(viz["series_colors"]) < 2 or not viz["view_markers"]:
        raise ValueError("Provide at least two series colors and one acquisition view marker")
    if any(not is_color_like(color) for color in viz["series_colors"]):
        raise ValueError("Invalid visualization.series_colors")
    cohorts = viz.get("display_cohorts")
    if cohorts is not None and (not isinstance(cohorts, list) or not cohorts
                               or not set(cohorts) <= {"available_week", "fixed_window"}):
        raise ValueError("Invalid visualization.display_cohorts")
    for name, options in viz["metrics"].items():
        role, scale = color_policy(viz, name)
        if role not in maps:
            raise ValueError(f"Unknown color mode for {name}: {role}")
        if scale["scope"] != "batch" or scale["mode"] not in ("auto", "fixed"):
            raise ValueError(f"Invalid color scale for {name}")
        if scale["mode"] == "fixed":
            limits = scale.get("limits")
            if (not isinstance(limits, (list, tuple)) or len(limits) != 2
                    or any(isinstance(x, bool) or not isinstance(x, (int, float))
                           or not isfinite(x) for x in limits) or limits[0] >= limits[1]):
                raise ValueError(f"Fixed color scale for {name} requires finite, increasing limits")
            if role == "diverging" and not limits[0] == -limits[1]:
                raise ValueError(f"Diverging fixed limits for {name} must be symmetric around zero")
        elif "limits" in options.get("color_scale", {}):
            raise ValueError(f"Auto color scale cannot specify limits for {name}")


def displayed_weeks(weeks, viz, reference):
    weeks = list(weeks)
    if viz["time_start"] == "start":
        return weeks
    if reference not in weeks:
        raise ValueError(f"Plot baseline {reference!r} is absent from saved weeks")
    return weeks[weeks.index(reference):]


def filter_cohorts(frame, viz):
    cohorts = viz.get("display_cohorts")
    return frame if cohorts is None or "cohort" not in frame else frame[frame.cohort.isin(cohorts)]


def saved_config(folder):
    import json
    path = Path(folder) / "metadata/analysis_config.json"
    return json.loads(path.read_text(encoding="utf-8"))["config"]
