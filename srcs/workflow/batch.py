"""Coordinate series workers and publish their products under one batch root."""
from contextlib import nullcontext
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

from srcs.core.identity import digest
from srcs.workflow.batch_layout import BatchLayout, dataset_label, series_label


def timestamp():
    return datetime.now(timezone.utc).isoformat()


def batch_functions(contexts):
    enabled = [{plugin.spec["id"]: plugin for plugin in ctx.workflow.functions} for ctx in contexts]
    for name, plugin in enabled[0].items():
        spec = plugin.spec.get("batch")
        if not spec or not all(name in item for item in enabled):
            continue
        source = spec.get("source")
        if source and not all(ctx.workflow.source and ctx.workflow.source.spec["id"] == source for ctx in contexts):
            continue
        yield name, plugin, spec


def is_an_pair(scopes):
    return (len(scopes) == 2 and {scope["series"] for scope in scopes} == {"Aseries", "Nseries"}
            and scopes[0]["branch"] == scopes[1]["branch"])


def run_batch(configs, *, project_root, stage=None, start=None, dry_run=False,
              import_source=None, report_section="all", output_root=None, resume_batch=None,
              execute=None, tag=None, enabled_functions=()):
    root = Path(project_root).resolve()
    from srcs.core.runtime import exclusive_run, select_steps
    from srcs.utils.config import load_context
    from srcs.workflow import REPORT_SECTIONS

    if import_source and (len(configs) != 1 or (stage or start) != "edges"):
        raise ValueError("FC import requires one configuration and an edges start")
    configs = [(root / config).resolve() for config in configs]
    contexts = [load_context(path, create_dirs=False, enabled_functions=enabled_functions) for path in configs]
    if import_source and contexts[0].workflow.source:
        raise ValueError("Cannot import FC while an external source is enabled")
    for context in contexts:
        context.args.tag = tag
    workflows = [context.workflow for context in contexts]
    if report_section not in ("all", *REPORT_SECTIONS) and not all(
            report_section in workflow.reports for workflow in workflows):
        raise ValueError(f"Unknown or disabled report section: {report_section}")
    plans = []
    for workflow in workflows:
        start_at = start or (workflow.source.spec["id"] if workflow.source else "fc")
        stages = select_steps(workflow.order, stage=stage, start=None if stage else start_at)
        workflow.validate_plan(stages)
        plans.append(stages)
    full_report = report_section == "all" and all("report" in plan for plan in plans)
    scopes = [context.data_identity["scope"] for context in contexts]
    can_compare = full_report and is_an_pair(scopes)
    functions = list(batch_functions(contexts)) if full_report and not stage and (len(contexts) == 1 or can_compare) else []
    names = [context.analysis_tag for context in contexts]
    if len(set(name.casefold() for name in names)) != len(names):
        raise ValueError("Two configurations resolve to the same series directory; use distinct configurations")
    configured_roots = {Path(context.output_root).resolve() for context in contexts}
    if resume_batch is not None and output_root is not None:
        raise ValueError("--resume-batch cannot be combined with --output-root")
    if output_root is None and resume_batch is None and len(configured_roots) != 1:
        raise ValueError("All configurations must share paths.output, or specify --output-root")
    parent = (root / output_root).resolve() if output_root is not None else Path(contexts[0].output_root)
    task = stage or (f"From-{start}" if start and start not in ("seedfc", "fc", "preprocess") else "Core")
    if not stage and any("geometry" in plan for plan in plans):
        task += "-Geometry"
    if tag and tag.casefold() not in task.casefold().split("-"):
        task += "-" + tag
    if resume_batch is not None:
        if import_source is not None:
            raise ValueError("Import FC into a fresh batch; do not replace its source during resume")
        layout = BatchLayout((root / resume_batch).resolve())
        if not layout.summary_path.is_file():
            raise FileNotFoundError(f"No batch summary: {layout.summary_path}")
    else:
        allocate = BatchLayout.preview if dry_run else BatchLayout.create
        layout = allocate(parent, dataset_label(scopes), task)
    outputs = [layout.root / name for name in names]
    for context, target in zip(contexts, outputs):
        context.args.series_output = target

    supplied = [series_label(scope) for scope in scopes]
    log_names = [f"{label}-{index:02d}" if supplied.count(label) > 1 else label
                 for index, label in enumerate(supplied, 1)]
    series = [dict(config=str(config), config_digest=digest(context.config), path=output.name,
                   log=layout.relative(layout.logger_path(label)), status="planned",
                   stages=[dict(stage=current, status="planned") for current in plan])
              for config, context, output, label, plan in zip(configs, contexts, outputs, log_names, plans)]
    batch = dict(status="planned", kind="pipeline", started=timestamp(), project_root=str(root),
                 python=sys.executable, series=series, report_section=report_section,
                 functions={name: dict(status="planned", inputs=names, results=[]) for name, _, _ in functions})
    batch["comparison"] = (dict(status="planned", path="compare") if can_compare else
                           dict(status="skipped", reason="Requires a full Aseries/Nseries report for the same processing branch"))
    execute = execute or subprocess.run
    lock = nullcontext() if dry_run else exclusive_run(layout.root / ".batch.lock")
    with lock:
        if resume_batch is not None:
            previous = json.loads(layout.summary_path.read_text(encoding="utf-8"))
            old = {item["path"]: item for item in previous.get("series", []) if "path" in item}
            if previous.get("schema") != "workflow_batch_v2" or previous.get("kind") != "pipeline" or any(
                    item["path"] not in old or item["config_digest"] != old[item["path"]].get("config_digest")
                    for item in series):
                raise ValueError("Resume requires the same effective configuration, series identity and tag as the saved batch")
            if not dry_run:
                history = layout.root / "metadata" / "previous_runs"
                copies = (history / f"attempt_{index:02d}.json" for index in range(1, 10000))
                from srcs.core.artifacts import write_json
                write_json(next(path for path in copies if not path.exists()), previous)
                updated = {item["path"]: item for item in series}
                batch["series"] = [updated.get(name, item) for name, item in old.items()]
                if not can_compare and previous.get("comparison", {}).get("status") == "complete":
                    batch["comparison"] = previous["comparison"]
                batch["functions"] = {**previous.get("functions", {}), **batch["functions"]}
        return _execute(layout, batch, root, series, outputs, contexts, functions, execute,
                        can_compare=can_compare, dry_run=dry_run, tag=tag, import_source=import_source,
                        enabled_functions=enabled_functions)


def _execute(layout, batch, root, series, outputs, contexts, functions, execute, *,
             can_compare, dry_run, tag, import_source, enabled_functions):
    from srcs.utils.logger import setup_logger, close_log_file

    print(f"Batch output: {layout.root}", flush=True)
    logger = setup_logger() if dry_run else setup_logger(str(layout.logger_path("batch")))
    section = batch["report_section"]
    active = None

    def save():
        if not dry_run:
            layout.save_summary(batch)

    def fail(status, reason):
        batch.update(status=status, reason=reason, finished=timestamp())
        if active is not None:
            active.update(status=status, reason=reason, finished=batch["finished"])
        for entry in series:
            if entry["status"] == "running":
                entry.update(status=status, reason=reason)
        save()

    try:
        if not dry_run:
            batch["status"] = "running"
            save()
        for index, (entry, output) in enumerate(zip(series, outputs), 1):
            logger.info("[%s/%s] %s -> %s", index, len(series), Path(entry["config"]).name, output.name)
            logger.info("Stages: %s", ", ".join(step["stage"] for step in entry["stages"]))
            if not dry_run:
                entry["status"] = "running"
            for step in entry["stages"][:1] if dry_run else entry["stages"]:
                current = step["stage"]
                command = [sys.executable, "-u", str(root / "run.py"), "--worker", "-c", entry["config"],
                           "--stage", current, "--report-section", section, "--series-output", str(output)]
                if dry_run:
                    command.append("--dry-run")
                else:
                    command.extend(["--log-file", str(layout.root / entry["log"])])
                if tag:
                    command.extend(["--tag", tag])
                for name in enabled_functions:
                    command.extend(["--function", name])
                if import_source and current == "edges":
                    command.extend(["--import-fc", str((root / import_source).resolve())])
                active = step
                step.update(status="running", started=timestamp())
                save()
                logger.info("%s ... %s", current, "validating" if dry_run else "running")
                result = execute(command, cwd=root)
                step["returncode"] = result.returncode
                if result.returncode:
                    fail("failed", f"{current} exited with {result.returncode}")
                    logger.error("Batch failed: %s", batch["reason"])
                    return result.returncode, None if dry_run else layout.summary_path
                step.update(status="planned" if dry_run else "complete", finished=timestamp())
                active = None
                save()
            entry["status"] = "planned" if dry_run else "complete"
            save()
        if dry_run:
            if can_compare:
                logger.info("A/N comparison target: %s (if shared products are available)", layout.compare)
            for name, _, _ in functions:
                logger.info("Batch function %s -> %s", name, layout.compare if len(outputs) > 1 else outputs[0])
            return 0, None

        if can_compare:
            from srcs.workflow.compare import REPORT_PRODUCTS, compare
            state = batch["comparison"]
            if any(all((output / "results" / name).is_file() for output in outputs) for name in REPORT_PRODUCTS):
                active = state
                state.update(status="running", started=timestamp())
                save()
                function_results = {}
                comparison = compare(*outputs, configs=[context.config for context in contexts],
                                     destination=layout.compare, function_results=function_results)
                state.update(status="complete", result=layout.relative(comparison), finished=timestamp(),
                             functions={key: layout.relative(path) for key, path in function_results.items()})
                active = None
            else:
                state.update(status="skipped", reason="No shared saved comparison products")
            save()
        for name, plugin, spec in functions:
            active = batch["functions"][name]
            active.update(status="running", started=timestamp())
            save()
            run = plugin.entry(spec["run"])
            for result in run(outputs, contexts, destination=layout.compare if len(outputs) > 1 else outputs[0]):
                active["results"].append({key: layout.relative(value) if key in {"path", "numeric"} else value
                                          for key, value in result.items()})
                save()
            active.update(status="complete", finished=timestamp())
            active = None
            save()
        batch.update(status="complete", finished=timestamp())
        save()
        logger.info("Batch complete: %s", layout.summary_path)
        return 0, layout.summary_path
    except BaseException as exc:
        fail("interrupted" if isinstance(exc, KeyboardInterrupt) else "failed", str(exc))
        logger.error("Batch %s: %s", batch["status"], exc)
        raise
    finally:
        if not dry_run:
            close_log_file()
