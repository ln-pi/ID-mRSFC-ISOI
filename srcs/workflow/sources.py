"""Resolve external FC sources and identify their numerical implementation."""
from srcs.core.identity import code_closure, digest
from srcs.workflow.plugins import enabled_plugins


def active_source(config, root):
    enabled = enabled_plugins(config, "source", root)
    if len(enabled) > 1:
        raise ValueError("Select one external source")
    if not enabled:
        return None
    source = enabled[0]
    source.entry(source.spec["validate_options"])(config, source.options)
    return source


def source_code(source, root):
    spec = source.spec
    numerical = spec["numerics"]
    listed = sorted(p.relative_to(source.directory).as_posix() for p in source.directory.rglob("*.py"))
    if sorted(numerical) != listed:
        raise ValueError("Source numerics must cover each implementation file once")
    closure = {}
    for filename in sorted(numerical):
        parts = filename.removesuffix(".py").split("/")
        if parts[-1] == "__init__":
            parts.pop()
        module = ".".join([source.module, *parts])
        closure.update(code_closure(root, module))
    for module in spec.get("kernel", []):
        closure.update(code_closure(root, module))
    return digest([spec, closure])
