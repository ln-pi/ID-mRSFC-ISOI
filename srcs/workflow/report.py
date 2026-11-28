"""Report dispatch from saved products; no source FC or edge cache is required."""
from importlib import import_module
from pathlib import Path
import json

from srcs.utils.provenance import code_identity, sha256
from srcs.utils.storage import save_json
from srcs.core.artifacts import artifact_entries, artifact_path, publish
from srcs.core.contracts import StageResult
from srcs.core.identity import digest
from srcs.workflow import REPORT_SECTIONS


def redraw_saved(folder, section, dry_run=False, visualization_config=None, functions=(), *, output_root=None):
    from argparse import Namespace
    from srcs.core.artifacts import write_json
    from srcs.utils.config import PROJECT_ROOT, RunContext, load_config
    from srcs.utils.logger import setup_logger, close_log_file
    from srcs.utils.visualization import merge, settings, validate
    from srcs.workflow.batch import timestamp
    from srcs.workflow.batch_layout import BatchLayout, dataset_label, saved_scopes
    from srcs.workflow.plugins import enable_functions

    folder = Path(folder).resolve()
    metadata = json.loads((folder / "metadata/analysis_config.json").read_text(encoding="utf-8"))
    config = metadata["config"]
    if visualization_config:
        override = load_config(PROJECT_ROOT / visualization_config)
        config = merge(config, {"visualization": override.get("visualization", override)})
    config = enable_functions(config, functions)
    validate(settings(config))
    parent = PROJECT_ROOT / (output_root or config["paths"]["output"])
    allocate = BatchLayout.preview if dry_run else BatchLayout.create
    scopes = list(saved_scopes([folder], [config]))
    # Validate the selected report section before allocating a batch.
    args = Namespace(saved_output=folder, report_section=section)
    ctx = RunContext(args, config, setup_logger())
    selected_sections(ctx, report_inputs(ctx))
    layout = allocate(parent, dataset_label(scopes), f"Redraw-{section}")
    target = layout.root / folder.name
    print(f"Batch output: {layout.root}")
    if dry_run:
        return 0
    ctx.args.redraw_root = target
    ctx.logger = setup_logger(str(layout.logger_path("batch")))
    summary = dict(kind="redraw", status="running", started=timestamp(), inputs=[str(folder)],
                   series=[dict(path=folder.name, report=f"{folder.name}/report", section=section)])
    layout.save_summary(summary)
    try:
        write_json(target / "metadata" / "redraw_source.json", {"source": str(folder), "config": config})
        run(ctx)
        summary["status"] = "complete"
        print(f"Report saved under {target / 'report'}")
        return 0
    except BaseException as exc:
        summary.update(status="interrupted" if isinstance(exc, KeyboardInterrupt) else "failed", reason=str(exc))
        ctx.logger.error("Redraw %s: %s", summary["status"], exc)
        raise
    finally:
        summary["finished"] = timestamp()
        layout.save_summary(summary)
        close_log_file()


def report_inputs(ctx):
    root = Path(ctx.output_path)
    inputs = dict(pairwise=Path(ctx.output_dir("results","pairwise",create=False))/"selection.json",
                descriptive=root / "results/descriptive/index.json",
                stability=Path(ctx.output_dir("results","stability",create=False))/"index.json",
                longitudinal=Path(ctx.output_dir("results","longitudinal",create=False))/"index.json",
                combined=Path(ctx.statistics_path),
                statistics=Path(ctx.output_dir("results","statistics",create=False))/"individual_identifiability.csv")
    inputs.update({name: root / spec["product"] for name, (_, spec) in report_plugins(ctx).items()})
    return inputs


def report_plugins(ctx):
    if hasattr(ctx.args, "report_plugins"):
        return ctx.args.report_plugins
    reports = dict(ctx.workflow.reports)
    name = getattr(ctx.args, "report_section", "all")
    if getattr(ctx.args, "saved_output", None) and name not in ("all", *REPORT_SECTIONS) and name not in reports:
        from srcs.utils.config import PROJECT_ROOT
        from srcs.workflow.plugins import load_plugin
        plugin = load_plugin("function", name, PROJECT_ROOT, optional=True)
        if plugin and name in plugin.spec.get("reports", {}):
            reports[name] = (plugin, plugin.spec["reports"][name])
    ctx.args.report_plugins = reports
    return reports


def selected_sections(ctx, inputs):
    section = getattr(ctx.args, "report_section", "all")
    if section == "all":
        return list(inputs)
    if section not in inputs:
        raise ValueError(f"Unknown or disabled report section {section}; available: {', '.join(inputs)}")
    return [section]


def run(ctx):
    section = getattr(ctx.args, "report_section", "all")
    inputs = report_inputs(ctx)
    selected = selected_sections(ctx, inputs)
    plugins = report_plugins(ctx)
    from inspect import getfile

    plugin_renderers = {name: sha256(getfile(plugin.entry(spec["render"])))
                        for name, (plugin, spec) in plugins.items() if name in selected}
    render_key = digest([ctx.viz(), ctx.section("report"), ctx.section("descriptive"), section,
                         {name: sha256(inputs[name]) for name in selected if inputs[name].is_file()},
                         plugin_renderers])
    if not getattr(ctx.args, "saved_output", None):
        previous = Path(ctx.output_dir("report", create=False)) / "index.json"
        if previous.is_file():
            previous_key = json.loads(previous.read_text(encoding="utf-8")).get("render_key")
            if previous_key != render_key:
                variants = Path(ctx.output_path) / "report_variants"
                for index in range(1, 10000):
                    suffix = "" if index == 1 else f"_{index:02d}"
                    variant = variants / f"{section}{suffix}"
                    marker = variant / "report/index.json"
                    if not variant.exists() or (marker.is_file() and
                            json.loads(marker.read_text(encoding="utf-8")).get("render_key") == render_key):
                        ctx.args.redraw_root = variant
                        break
                else:
                    raise ValueError(f"No free report variant for {section}")
    rendered = []
    figures, skipped = [], []
    for name in selected:
        marker = Path(ctx.output_path)/"metadata/stages"/(name+".json")
        if marker.is_file():
            state=json.loads(marker.read_text(encoding="utf-8"))
            if state["status"]!="complete":
                skipped.append(name)
                continue
            for item in artifact_entries(state["cache"]["outputs"], marker.parent):
                path = item["path"]
                if not path.is_file() or sha256(path) != item["sha256"]:
                    raise ValueError(f"Saved product changed: {path}; rerun {name}")
        if not inputs[name].is_file():
            skipped.append(name)
            continue
        if name in plugins:
            plugin, spec = plugins[name]
            images, paths = plugin.entry(spec["render"])(ctx, inputs[name])
            figures.extend(images)
        else:
            result = import_module("srcs.visualization." + name).run(ctx)
            paths = result if isinstance(result, list) else [str(inputs[name])]
            figures.extend(sorted(Path(ctx.output_dir("report", name, create=False)).rglob("*.png")))
        rendered.append(dict(section=name, inputs=[dict(path=str(p), sha256=sha256(p)) for p in paths]))
    from srcs.visualization.report import publish_report
    folder = Path(ctx.output_dir("report"))
    internal = not getattr(ctx.args, "saved_output", None)
    if internal:
        rendered = [dict(item, inputs=[dict(source, path=artifact_path(source["path"], folder))
                                       for source in item["inputs"]]) for item in rendered]
    publish_report(ctx, figures, rendered, skipped)
    save_json(folder / "index.json",
              dict(schema="unified_report", sections=rendered, render_key=render_key, code=code_identity()))
    ref=publish(folder/"products.json","unified_report","unified_report",
                [*figures,*(folder/name for name in ("report.html","report.pdf","manifest.json","index.json"))])
    if internal:
        stage_dir = Path(ctx.output_path) / "metadata/stages"
        stage_inputs = [dict(item, inputs=[dict(source, path=artifact_path(folder / source["path"], stage_dir))
                                          for source in item["inputs"]]) for item in rendered]
        save_json(stage_dir / "report.json",
                  dict(schema="stage_record", stage="report", status="complete", options_digest=render_key,
                       inputs=stage_inputs, products=[dict(name=ref.name, schema=ref.schema,
                           path=artifact_path(ref.path, stage_dir), sha256=ref.sha256)]))
    ctx.logger.info("Report sections: %s", ", ".join(r["section"] for r in rendered))
    return StageResult((ref,))
