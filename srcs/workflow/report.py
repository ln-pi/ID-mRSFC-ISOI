"""Report dispatch from saved products; no source FC or edge cache is required."""
from importlib import import_module
from pathlib import Path
import json

from srcs.utils.provenance import code_identity, sha256
from srcs.utils.storage import save_json
from srcs.core.artifacts import publish
from srcs.core.contracts import StageResult
from srcs.core.identity import digest


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
                pairwise=Path(ctx.output_dir("results","pairwise",create=False))/"selection.json",
                descriptive=root / "results/descriptive/v2/index.json",
                stability=Path(ctx.output_dir("results","stability",create=False))/"index.json",
                longitudinal=Path(ctx.output_dir("results","longitudinal",create=False))/"index.json",
                combined=Path(ctx.statistics_path),
                statistics=Path(ctx.output_dir("results","statistics",create=False))/"individual_identifiability.csv")


def run(ctx):
    section = getattr(ctx.args, "report_section", "all")
    inputs = report_inputs(ctx)
    selected = list(inputs) if section == "all" else [section]
    rendered = []
    figures, skipped = [], []
    for name in selected:
        marker = Path(ctx.output_path)/"metadata/stages"/(name+".json")
        if marker.is_file():
            state=json.loads(marker.read_text(encoding="utf-8"))
            if state["status"]!="complete":
                skipped.append(name)
                continue
            for path,expected in state["cache"]["outputs"].items():
                if not Path(path).is_file() or sha256(path)!=expected:
                    raise ValueError(f"Saved product changed: {path}; rerun {name}")
        if not inputs[name].is_file():
            skipped.append(name)
            continue
        if name == "geometry":
            from srcs.visualization.geometry import render
            figures.extend(render(ctx))
            index=json.loads(inputs[name].read_text(encoding="utf-8"))
            rendered.append(dict(section=name, inputs=[dict(path=p,sha256=h) for p,h in
                {**index["files"],str(inputs[name]):sha256(inputs[name])}.items()]))
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
    folder=Path(ctx.output_dir("report"))
    ref=publish(folder/"products.json","unified_report_v1","unified_report_v1",
                [*figures,*(folder/name for name in ("report.html","report.pdf","manifest.json","index.json"))])
    if not getattr(ctx.args,"saved_output",None):
        save_json(Path(ctx.output_path)/"metadata/stages/report.json",
                  dict(schema="stage_record_v1",stage="report",status="complete",
                       options_digest=digest([ctx.viz(),ctx.section("report"),section]),
                       inputs=rendered,products=[dict(name=ref.name,schema=ref.schema,path=str(ref.path),sha256=ref.sha256)]))
    ctx.logger.info("Report sections: %s", ", ".join(r["section"] for r in rendered))
    return StageResult((ref,))
