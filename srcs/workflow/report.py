"""Report dispatch from saved products; no source FC or edge cache is required."""
from importlib import import_module
from pathlib import Path

from srcs.utils.provenance import code_identity, sha256
from srcs.utils.results import save_json


def report_inputs(ctx):
    root = Path(ctx.output_path)
    return dict(pairwise=root / "results/pairwise/selection.json",
                descriptive=root / "results/descriptive/v2/index.json",
                stability=root / "results/stability/index.json",
                longitudinal=root / "results/longitudinal/index.json",
                combined=Path(ctx.statistics_path),
                statistics=root / "results/statistics/individual_identifiability.csv")


def run(ctx):
    section = getattr(ctx.args, "report_section", "all")
    inputs = report_inputs(ctx)
    selected = list(inputs) if section == "all" else [section]
    rendered = []
    for name in selected:
        if not inputs[name].is_file():
            if section != "all":
                raise FileNotFoundError(f"Missing {name} numerical input: {inputs[name]}")
            continue
        result = import_module("srcs.visualization." + name).run(ctx)
        paths = result if isinstance(result, list) else [str(inputs[name])]
        rendered.append(dict(section=name, inputs=[dict(path=p, sha256=sha256(p)) for p in paths]))
    if not rendered:
        raise FileNotFoundError("No saved numerical products available for report")
    save_json(Path(ctx.output_dir("report")) / "index.json",
              dict(schema="unified_report_v1", sections=rendered, code=code_identity()))
    ctx.logger.info("Report sections: %s", ", ".join(r["section"] for r in rendered))
