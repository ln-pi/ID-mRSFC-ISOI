"""Report dispatch from saved products; no source FC or edge cache is required."""
from importlib import import_module
from pathlib import Path
import json

from srcs.utils.provenance import code_identity, sha256
from srcs.utils.results import save_json


def redraw_saved(folder, section, dry_run=False):
    from argparse import Namespace
    from srcs.utils.config import RunContext
    from srcs.utils.logger import setup_logger
    metadata = json.loads((folder/"metadata/analysis_config.json").read_text(encoding="utf-8"))
    from srcs.utils.config import PROJECT_ROOT
    from uuid import uuid4
    target = PROJECT_ROOT/"reports"/"redraw"/(folder.name+"_"+uuid4().hex[:8])
    ctx = RunContext(Namespace(saved_output=folder,report_section=section,redraw_root=target),metadata["config"],setup_logger())
    if not dry_run:
        run(ctx)
        print(f"Report saved under {target}")
    return 0


def report_inputs(ctx):
    root = Path(ctx.output_path)
    return dict(geometry=root / "results/geometry/v1/index.json",
                pairwise=root / "results/pairwise/selection.json",
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
    figures, skipped = [], []
    for name in selected:
        if not inputs[name].is_file():
            skipped.append(name)
            continue
        if name == "geometry":
            from srcs.visualization.geometry import render
            figures.extend(render(ctx))
            rendered.append(dict(section=name, inputs=[dict(path=str(inputs[name]),sha256=sha256(inputs[name]))]))
            continue
        result = import_module("srcs.visualization." + name).run(ctx)
        paths = result if isinstance(result, list) else [str(inputs[name])]
        rendered.append(dict(section=name, inputs=[dict(path=p, sha256=sha256(p)) for p in paths]))
    from srcs.visualization.geometry import publish_report
    for item in rendered:
        if item["section"] != "geometry":
            figures.extend(sorted(Path(ctx.output_dir("report",item["section"],create=False)).rglob("*.png")))
    publish_report(ctx,figures,rendered,skipped)
    save_json(Path(ctx.output_dir("report")) / "index.json",
              dict(schema="unified_report_v1", sections=rendered, code=code_identity()))
    ctx.logger.info("Report sections: %s", ", ".join(r["section"] for r in rendered))
