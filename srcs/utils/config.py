import os
import hashlib
import json
from dataclasses import dataclass
from typing import Any

import yaml

from srcs.utils.logger import setup_logger
from srcs.utils.utils import parse_args, split_session_indices


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
    def sessions(self):
        return self.variables["Session_ids"]

    @property
    def mice(self):
        return self.variables["Mouse_ids"]

    @property
    def trials(self):
        return self.variables["Trial_ids"]

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

    @property
    def repeat_session(self):
        repeat = self.section("ID").get("repeat_session")
        if isinstance(repeat, str):
            if repeat not in self.sessions:
                raise ValueError(f"Unknown repeat_session: {repeat}")
            repeat = self.sessions.index(repeat) + 1
        return repeat

    @property
    def session_groups(self):
        groups = split_session_indices(len(self.sessions), self.group_mode,
                                       self.repeat_session, self.session_seed)
        return [[self.sessions[i] for i in group] for group in groups]

    @property
    def analysis_tag(self):
        edge_seed = self.section("Edge").get("seed", 42)
        if isinstance(edge_seed, bool) or not isinstance(edge_seed, int) or edge_seed < 0:
            raise ValueError("Edge.seed must be a nonnegative integer for reproducible analysis.")
        identity = {
            "variables": self.variables,
            "raw": self.raw_root,
            "processed": self.processed_root,
            "groups": self.session_groups,
            "ID": self.section("ID"),
            "Edge": self.section("Edge"),
        }
        digest = hashlib.sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:10]
        repeat = f"_repeat{self.repeat_session}" if self.repeat_session is not None else ""
        return f"seed_{self.group_mode}_{self.session_seed}{repeat}_edge{edge_seed}_{digest}"

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
        return os.path.join(self.output_path, self.statistics_file)

    def save_analysis_metadata(self):
        from srcs.utils.results import save_json
        save_json(os.path.join(self.output_path, "analysis_config.json"), {
            "config": self.config, "session_groups": self.session_groups,
            "session_seed": self.session_seed, "output_path": self.output_path,
            "trial_order": self.trials,
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


def load_context(config_path=None, create_dirs=True):
    args = parse_args()
    if config_path is not None:
        args.config = config_path

    config = load_config(args.config)
    if config is None:
        raise ValueError("Config is empty (None). Exiting.")

    logger = setup_logger(config["paths"]["log"])
    ctx = RunContext(args=args, config=config, logger=logger)

    if create_dirs:
        os.makedirs(ctx.processed_path, exist_ok=True)
        os.makedirs(ctx.output_path, exist_ok=True)

    return ctx
