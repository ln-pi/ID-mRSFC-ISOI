import os
import hashlib
import json
from dataclasses import dataclass
from typing import Any

import yaml
from srcs.utils.output_paths import DIRECTORIES

from srcs.utils.logger import setup_logger
from srcs.legacy.averaging import split_group_indices


@dataclass
class RunContext:
    args: Any
    config: dict[str, Any]
    logger: Any

    @property
    def paths(self):
        return self.config["paths"]

    @property
    def variables(self):
        return self.config["variables"]

    @property
    def raw_root(self):
        return self.paths["raw"]

    @property
    def processed_root(self):
        return self.paths["processed"]

    @property
    def output_root(self):
        return self.paths["output"]

    @property
    def mask_root(self):
        return self.paths["mask"]

    @property
    def dataset(self):
        return self.variables["dataset"]

    @property
    def weeks(self):
        return self._ids("Week_ids")

    def _ids(self, name):
        values = self.variables.get(name)
        if not isinstance(values, list) or not values or any(
            not isinstance(v, str) or not v or v in (".", "..") or "/" in v or "\\" in v
            for v in values
        ) or len(set(values)) != len(values):
            raise ValueError(f"variables.{name} requires unique, nonempty directory IDs")
        return values

    @property
    def mice(self):
        return self._ids("Mouse_ids")

    @property
    def runs(self):
        return self._ids("Run_ids")

    @property
    def split_axis(self):
        value = self.section("analysis").get("split_axis", "run")
        if value not in ("week", "run"):
            raise ValueError("analysis.split_axis must be week or run")
        return value

    @property
    def selection_mode(self):
        value = self.section("analysis").get("selection_mode", "cumulative")
        if value not in ("single", "cumulative"):
            raise ValueError("analysis.selection_mode must be single or cumulative")
        return value

    @property
    def atlases(self):
        return self.variables["Atlas"]

    @property
    def processed_path(self):
        return os.path.join(self.processed_root, self.dataset)

    @property
    def output_path(self):
        return os.path.join(self.output_root, self.analysis_tag, self.dataset)

    @property
    def session_seed(self):
        seed = self.section("ID").get("session_seed", 42)
        if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
            raise ValueError("ID.session_seed must be a nonnegative integer so 03/04/05 share one output directory.")
        return seed

    @property
    def group_mode(self):
        return self.section("ID").get("group_mode", "random")

    def group_indices(self, weeks=None, runs=None):
        ids = (self.weeks if weeks is None else weeks) if self.split_axis == "week" else (self.runs if runs is None else runs)
        repeat = self.section("ID").get("repeat_id")
        if isinstance(repeat, str):
            if repeat not in ids:
                raise ValueError(f"Unknown ID.repeat_id on selected {self.split_axis} axis: {repeat}")
            repeat = ids.index(repeat) + 1
        return split_group_indices(len(ids), self.group_mode, repeat, self.session_seed)

    def split_groups(self, weeks=None, runs=None):
        ids = (self.weeks if weeks is None else weeks) if self.split_axis == "week" else (self.runs if runs is None else runs)
        return [[ids[i] for i in group] for group in self.group_indices(weeks, runs)]

    @property
    def analysis_tag(self):
        edge_seed = self.section("Edge").get("seed", 42)
        if isinstance(edge_seed, bool) or not isinstance(edge_seed, int) or edge_seed < 0:
            raise ValueError("Edge.seed must be a nonnegative integer for reproducible analysis.")
        identity = self.config  # Deterministic isolation by the complete parsed configuration.
        digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:10]
        repeat_id = self.section("ID").get("repeat_id")
        repeat = f"_repeat{repeat_id}" if repeat_id is not None else ""
        return f"seed_{self.group_mode}_{self.session_seed}{repeat}_edge{edge_seed}_{digest}"

    @property
    def correlation_path(self):
        return os.path.join(self.output_path, "Correlation")

    @property
    def statistics_file(self):
        # Legacy supplement-only configurations still resolve to one shared name.
        name = self.section("analysis").get(
            "statistics_file", self.section("supplement").get("statistics_file", "statistics.pkl"))
        if not isinstance(name, str) or not name or os.path.basename(name) != name:
            raise ValueError("analysis.statistics_file must be a filename without directories.")
        return name

    @property
    def statistics_path(self):
        return os.path.join(self.output_dir('Identification', create=False), self.statistics_file)

    def save_analysis_metadata(self, selection=None):
        from srcs.utils.io import save_json
        save_json(os.path.join(self.output_path, "metadata", "analysis_config.json"), {
            "config": self.config, "split_axis": self.split_axis,
            "selection_mode": self.selection_mode, "session_seed": self.session_seed,
            "week_order": self.weeks, "run_order": self.runs,
            "output_path": self.output_path, "selection": selection,
        })

    def section(self, name, default=None):
        return self.config.get(name, default or {})

    def viz(self, name=None, default=None):
        visualization = self.config.get("visualization", {})
        if name is None:
            return visualization
        return visualization.get(name, default or {})

    def viz_style(self, *names, default=None):
        style = dict(default or {})
        for name in names:
            style.update(self.viz(name))
        return style

    def render_options(self):
        defaults = self.viz("defaults")
        return {
            "show": defaults.get("show", False),
            "save": defaults.get("save", True),
            "log": defaults.get("log", False),
        }

    def figure_dpi(self):
        return self.viz("defaults").get("dpi", 180)

    def render_dpi(self):
        return self.viz("defaults").get("render_dpi", self.figure_dpi())

    def output_dir(self, *parts, create=True):
        if parts:
            parts = (DIRECTORIES.get(parts[0], parts[0]), *parts[1:])
        path = os.path.join(self.output_path, *parts)
        if create:
            os.makedirs(path, exist_ok=True)
        return path

    def processed_dir(self, *parts, create=True):
        path = os.path.join(self.processed_path, *parts)
        if create:
            os.makedirs(path, exist_ok=True)
        return path


def load_config(path):
    try:
        with open(path, "r", encoding="utf-8") as f:
            return yaml.safe_load(f)
    except Exception as e:
        setup_logger().error(f"Failed to load config from {path}: {e}")
        return None


def load_context(config_path="configs/default.yaml", create_dirs=False):
    from types import SimpleNamespace
    config = load_config(config_path)
    if not isinstance(config, dict):
        raise ValueError(f"Invalid configuration: {config_path}")
    ctx = RunContext(SimpleNamespace(config=str(config_path)), config, setup_logger())
    _ = ctx.weeks, ctx.runs, ctx.mice, ctx.session_seed, ctx.analysis_tag
    if create_dirs:
        os.makedirs(ctx.output_path, exist_ok=True)
        ctx.logger = setup_logger(config["paths"]["log"])
    return ctx
