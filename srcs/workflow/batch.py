"""P2-style sequential series/stage coordination with durable batch summaries."""
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
from uuid import uuid4

from srcs.utils.storage import save_json


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


def run_batch(configs, *, project_root, stage=None, start=None, dry_run=False,
              import_source=None, report_section="all", report_root=None, execute=None, tag=None,
              enabled_functions=()):
    root = Path(project_root).resolve()
    from srcs.core.runtime import select_steps
    if import_source and (len(configs) != 1 or (stage or start) != "edges"):
        raise ValueError("FC import requires one configuration and an edges start")
    configs = [(root / config).resolve() for config in configs]
    for path in configs:
        if not path.is_file():
            raise FileNotFoundError(f"Configuration not found: {path}")
    from srcs.utils.config import load_context
    contexts = [load_context(path, create_dirs=False, enabled_functions=enabled_functions) for path in configs]
    for context in contexts:
        context.args.tag = tag
    workflows = [context.workflow for context in contexts]
    from srcs.workflow import REPORT_SECTIONS
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
    functions = list(batch_functions(contexts)) if full_report and not stage else []
    outputs = [Path(context.output_path) for context in contexts]
    execute = execute or subprocess.run
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%SZ")
    summary_root = root / report_root if report_root is not None else Path(contexts[0].output_root) / "_batches"
    summary = summary_root / f"{stamp}_{uuid4().hex[:8]}" / "batch_summary.json"
    batch = dict(schema="workflow_batch", status="running", project_root=str(root),
                 python=sys.executable, dry_run=dry_run, series=[])
    if functions:
        batch["functions"] = {name: dict(status="planned", inputs=[str(output) for output in outputs], results=[])
                              for name, _, _ in functions}
        if "seedfc_supplement" in batch["functions"]:
            batch["supplement"] = batch["functions"]["seedfc_supplement"]
    def save():
        save_json(summary, batch)
    save()
    try:
        for index, config in enumerate(configs, 1):
            workflow, stages = workflows[index - 1], plans[index - 1]
            if workflow.source and import_source:
                raise ValueError("Cannot import FC while an external source is enabled")
            batch["stages"] = stages
            entry = dict(config=str(config), status="running", stages=[])
            batch["series"].append(entry)
            print(f"[{index}/{len(configs)}] {config.name}", flush=True)
            # Dry runs validate only the first stage, then record the full plan.
            execution_stages = stages[:1] if dry_run else stages
            for current in execution_stages:
                command = [sys.executable, "-u", str(root / "run.py"), "--worker", "-c", str(config),
                           "--stage", current, "--report-section", report_section]
                if dry_run:
                    command.append("--dry-run")
                if tag:
                    command.extend(["--tag", tag])
                for name in enabled_functions:
                    command.extend(["--function", name])
                if import_source and current == "edges":
                    command.extend(["--import-fc", str((root / import_source).resolve())])
                step = dict(stage=current, status="running", started=datetime.now(timezone.utc).isoformat())
                entry["stages"].append(step)
                save()
                print(f"{current} ... {'validating' if dry_run else 'running'}", flush=True)
                result = execute(command, cwd=root)
                step.update(returncode=result.returncode,
                            status="planned" if dry_run and result.returncode == 0 else "complete" if result.returncode == 0 else "failed",
                            finished=datetime.now(timezone.utc).isoformat())
                if result.returncode:
                    entry.update(status="failed", reason=f"{current} exited with {result.returncode}")
                    batch["status"] = "failed"
                    save()
                    print(f"Batch failed; summary: {summary}", flush=True)
                    return result.returncode, summary
                save()
            if dry_run:
                entry["stages"].extend(dict(stage=s, status="deferred") for s in stages[1:])
            entry["status"] = "planned" if dry_run else "complete"
            save()
        if not dry_run and len(configs) == 2 and full_report:
            from srcs.workflow.compare import REPORT_PRODUCTS, compare
            identities = [context.data_identity["scope"] for context in contexts]
            if (identities[0]["series"] != identities[1]["series"]
                    and identities[0]["branch"] == identities[1]["branch"]
                    and any(all((output/"results"/name).is_file() for output in outputs) for name in REPORT_PRODUCTS)):
                function_results = {}
                comparison = compare(*outputs, configs=[context.config for context in contexts],
                                     function_results=function_results)
                batch["comparison"] = dict(status="complete", path=str(comparison), functions=function_results)
                save()
        if not dry_run:
            for name, plugin, spec in functions:
                state = batch["functions"][name]
                state["status"] = "running"
                save()
                run = plugin.entry(spec["run"])
                for result in run(outputs, contexts):
                    state["results"].append(result)
                    save()
                state["status"] = "complete"
                save()
    except BaseException as exc:
        status = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        batch.update(status=status, reason=str(exc))
        for state in batch.get("functions", {}).values():
            if state["status"] == "running":
                state.update(status=status, reason=str(exc))
        if batch["series"] and batch["series"][-1]["status"] == "running":
            active = batch["series"][-1]
            active.update(status=status, reason=str(exc))
            if active["stages"] and active["stages"][-1]["status"] == "running":
                active["stages"][-1].update(status=status, reason=str(exc))
        save()
        raise
    batch["status"] = "planned" if dry_run else "complete"
    save()
    print(f"Batch {batch['status']}; summary: {summary}", flush=True)
    return 0, summary
