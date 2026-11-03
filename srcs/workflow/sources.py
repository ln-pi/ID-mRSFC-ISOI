from pathlib import Path
import yaml

from srcs.core.identity import code_closure, digest


def active_source(config, root):
    enabled = []
    known = {}
    for path in sorted((Path(root)/"plugins").glob("*/plugin.yaml")):
        spec = yaml.safe_load(path.read_text(encoding="utf-8"))
        if spec["kind"] != "source" or spec["id"] != path.parent.name:
            raise ValueError(f"Invalid source plugin: {path}")
        known[spec["id"]] = spec
    for name,settings in config.get("plugins",{}).items():
        if not settings.get("enabled"):
            continue
        if name not in known:
            raise ValueError(f"Unknown source: {name}")
        enabled.append((known[name],settings["options"]))
    if len(enabled)>1:
        raise ValueError("Select one external source")
    if enabled:
        spec,options = enabled[0]
        if options.get("branch") not in ("GSR","noGSR") or not all(options.get(k) for k in ("source","series","partitions")):
            raise ValueError("Invalid SeedFC source options")
        if options["series"] != config["variables"].get("series"):
            raise ValueError("SeedFC series differs from configured series")
        return spec,options
    return None


def source_code(spec,root):
    folder = root/"plugins"/spec["id"]
    listed = sorted(p.relative_to(folder).as_posix() for p in folder.rglob("*.py"))
    if sorted(spec["numerics"]) != listed:
        raise ValueError("Source numerics must cover all Python files")
    return digest([spec,code_closure(root,"plugins."+spec["id"])])
