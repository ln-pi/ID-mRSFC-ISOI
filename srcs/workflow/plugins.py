"""Discover source and function plugins without importing numerical implementations."""
from collections.abc import Mapping
from dataclasses import dataclass
from importlib import import_module
from importlib.util import find_spec
from pathlib import Path
import re

import yaml


@dataclass(frozen=True)
class Plugin:
    spec: dict
    options: dict
    directory: Path

    @property
    def module(self):
        return f"plugins.{self.spec['kind']}.{self.spec['id']}"

    def entry(self, entry):
        module, name = entry.split(":")
        return getattr(import_module(f"{self.module}{'.' + module if module else ''}"), name)

    def stage_module(self, definition):
        suffix = definition.get("module")
        return f"{self.module}.{suffix}" if suffix else self.module


def configured(config, kind):
    plugins = config.get("plugins", {})
    settings = plugins.get(kind, {})
    if not isinstance(settings, Mapping):
        raise ValueError(f"plugins.{kind} requires a mapping")
    if kind == "source":
        legacy = {key: value for key, value in plugins.items() if key not in ("source", "function")}
        conflict = legacy.keys() & settings.keys()
        if conflict:
            raise ValueError(f"Duplicate source plugin configuration: {sorted(conflict)}")
        return {**legacy, **settings}
    if "seedfc_supplement" in config:
        legacy = config["seedfc_supplement"]
        return {"seedfc_supplement": {"enabled": legacy.get("enabled", False),
                "options": {key: value for key, value in legacy.items() if key != "enabled"}}, **settings}
    return settings


def function_options(config, plugin_id, subsection=None):
    options = configured(config, "function").get(plugin_id, {}).get("options", {})
    if subsection is None:
        return dict(options)
    return {**config.get(subsection, {}), **options.get(subsection, {})}


def load_plugin(kind, name, root, options=None, *, optional=False, dependencies=True):
    if not re.fullmatch(r"[a-z][a-z0-9_]*", name):
        raise ValueError(f"Invalid plugin ID: {name}")
    directory = Path(root) / "plugins" / kind / name
    path = directory / "plugin.yaml"
    if not path.is_file() and kind == "function":
        try:
            spec = find_spec(f"plugins.{kind}.{name}")
        except ModuleNotFoundError:
            spec = None
        if spec and spec.submodule_search_locations:
            directory = Path(next(iter(spec.submodule_search_locations)))
            path = directory / "plugin.yaml"
    if not path.is_file():
        if optional:
            return None
        raise ValueError(f"Enabled {kind} plugin is unavailable: {name}")
    spec = yaml.safe_load(path.read_text(encoding="utf-8"))
    if spec.get("kind") != kind or spec.get("id") != name:
        raise ValueError(f"Invalid {kind} plugin declaration: {path}")
    if dependencies:
        for module in spec.get("dependencies", []):
            if find_spec(module) is None:
                raise ValueError(f"Plugin {name} requires {module}")
    return Plugin(spec, dict(options or {}), directory)


def enabled_plugins(config, kind, root):
    plugins = []
    for name, setting in configured(config, kind).items():
        if not isinstance(setting, Mapping) or type(setting.get("enabled", False)) is not bool:
            raise ValueError(f"plugins.{kind}.{name}.enabled must be a boolean")
        options = setting.get("options", {})
        if not isinstance(options, Mapping):
            raise ValueError(f"plugins.{kind}.{name}.options requires a mapping")
        if setting.get("enabled", False):
            plugins.append(load_plugin(kind, name, root, options))
    return plugins



def enable_functions(config, names):
    if not names:
        return config
    plugins = dict(config.get("plugins", {}))
    functions = dict(configured(config, "function"))
    for name in names:
        functions[name] = {**functions.get(name, {}), "enabled": True}
    plugins["function"] = functions
    return {**config, "plugins": plugins}
