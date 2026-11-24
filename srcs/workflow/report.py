"""Report dispatch from saved products; no source FC or edge cache is required."""
from importlib import import_module
from pathlib import Path
import json

from srcs.utils.provenance import code_identity, sha256
from srcs.utils.storage import save_json
from srcs.core.artifacts import publish
from srcs.core.contracts import StageResult
from srcs.core.identity import digest
from srcs.workflow import REPORT_SECTIONS


def redraw_saved(folder, section, dry_run=False, visualization_config=None, functions=()):
    from argparse import Namespace
    from srcs.utils.config import RunContext
    from srcs.utils.logger import setup_logger
    metadata = json.loads((folder/"metadata/analysis_config.json").read_text(encoding="utf-8"))
    from srcs.utils.config import PROJECT_ROOT, load_config
    from uuid import uuid4
    defaults = load_config(PROJECT_ROOT / "configs/default.yaml")
    output_root = (PROJECT_ROOT / defaults["paths"]["output"]).resolve()
    target = output_root / "_redraw" / (folder.name + "_" + uuid4().hex[:8])
    from srcs.utils.visualization import merge, settings, validate
    config = metadata["config"]
    if visualization_config:
        override = load_config(PROJECT_ROOT / visualization_config)
        config = merge(config, {"visualization": override.get("visualization", override)})
    from srcs.workflow.plugins import enable_functions
    config = enable_functions(config, functions)
    validate(settings(config))
    if not dry_run:
        logger = setup_logger(str(output_root / "_logs" / "redraw.log"))
        args = Namespace(saved_output=folder, report_section=section, redraw_root=target)
        run(RunContext(args, config, logger))
        print(f"Report saved under {target}")
    return 0


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
                ctx.args.redraw_root = Path(ctx.output_path) / "rendered_reports" / render_key
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
        if name in plugins:
            plugin, spec = plugins[name]
            images, paths = plugin.entry(spec["render"])(ctx, inputs[name])
            figures.extend(images)
        else:
            result = import_module("srcs.visualization." + name).run(ctx)
            paths = result if isinstance(result, list) else [str(inputs[name])]
            figures.extend(sorted(Path(ctx.output_dir("report", name, create=False)).rglob("*.png")))
        rendered.append(dict(section=name, inputs=[dict(path=p, sha256=sha256(p)) for p in paths]))
    from srcs.visualization.report import publish_report
    publish_report(ctx,figures,rendered,skipped)
    save_json(Path(ctx.output_dir("report")) / "index.json",
              dict(schema="unified_report", sections=rendered, render_key=render_key, code=code_identity()))
    folder=Path(ctx.output_dir("report"))
    ref=publish(folder/"products.json","unified_report","unified_report",
                [*figures,*(folder/name for name in ("report.html","report.pdf","manifest.json","index.json"))])
    if not getattr(ctx.args,"saved_output",None):
        save_json(Path(ctx.output_path)/"metadata/stages/report.json",
                  dict(schema="stage_record",stage="report",status="complete",
                       options_digest=render_key,
                       inputs=rendered,products=[dict(name=ref.name,schema=ref.schema,path=str(ref.path),sha256=ref.sha256)]))
    ctx.logger.info("Report sections: %s", ", ".join(r["section"] for r in rendered))
    return StageResult((ref,))
