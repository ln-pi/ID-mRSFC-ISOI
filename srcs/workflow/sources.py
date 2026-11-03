from pathlib import Path
import yaml

from srcs.core.identity import code_closure, digest


def active_source(config, root):
    enabled = []
    known = {}
    for directory in sorted((Path(root)/"plugins").glob("*")):
        if not directory.is_dir() or directory.name=="__pycache__":
            continue
        path=directory/"plugin.yaml"
        if not path.is_file():
            raise ValueError(f"Missing source declaration: {directory}")
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
    closure={}
    for filename in listed:
        parts=list(Path(filename).with_suffix("").parts)
        if parts[-1]=="__init__":
            parts.pop()
        module=".".join(["plugins",spec["id"],*parts])
        closure.update(code_closure(root,module))
    for module in spec.get("kernel",[]):
        closure.update(code_closure(root,module))
    return digest([spec,closure])
