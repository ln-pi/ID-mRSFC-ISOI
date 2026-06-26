import os
from dataclasses import dataclass
from typing import Any

import yaml

from srcs.utils.logger import setup_logger
from srcs.utils.utils import parse_args


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
        return os.path.join(self.output_root, self.dataset)

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
