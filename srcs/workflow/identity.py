import re

from srcs.core.identity import code_closure, digest
from srcs.workflow import MODULES


def data_identity(ctx, root):
    options = {"preprocess": ctx.section("preprocess"), "fc": ctx.section("correlation"), "edges": {}}
    return {
        "schema": "data_identity_v1", "schema_version": 1,
        "scope": {
            "dataset": ctx.dataset, "series": ctx.variables.get("series", ""),
            "branch": "GSR" if ctx.section("preprocess").get("gsr", True) else "noGSR",
            "source": {"kind": "native", "id": "native", "digest": None},
            "weeks": ctx.weeks, "mice": ctx.mice, "runs": ctx.runs,
            "partitions": [{"name": "_".join(a), "exclude": a} for a in ctx.atlases],
            "paths": {key: ctx.paths[key] for key in ("raw", "processed", "mask")},
        },
        "stages": {name: {"kind": "builtin", "code": digest(code_closure(root, MODULES[name])),
                          "options": value} for name, value in options.items()},
    }


def data_tag(identity, suffix=None):
    scope = identity["scope"]
    prefix = "-".join(filter(None, (scope["dataset"], scope["series"],
                                    scope["branch"] if scope["branch"] != "GSR" else "")))
    parts = [prefix, f'{scope["weeks"][0]}-{scope["weeks"][-1]}',
             f'{len(scope["mice"])}mice', f'{len(scope["runs"])}run', digest(identity)[:8]]
    if suffix is not None:
        parts.append(suffix)
    if any(not re.fullmatch(r"[A-Za-z0-9._-]+", part) for part in parts if part):
        raise ValueError("Invalid data tag component")
    return "_".join(filter(None, parts))
