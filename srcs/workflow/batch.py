"""P2-style sequential series/stage coordination with durable batch summaries."""
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
from uuid import uuid4

from srcs.utils.storage import save_json
from srcs.workflow import STAGES


def run_batch(configs, *, project_root, stage=None, start="fc", dry_run=False,
              import_source=None, report_section="all", report_root="reports/automatic", execute=None, tag=None):
    root = Path(project_root).resolve()
    from srcs.core.runtime import select_steps
    from srcs.utils.config import load_config
    from srcs.workflow.sources import active_source
    if import_source and (len(configs) != 1 or (stage or start) != "edges"):
        raise ValueError("FC import requires one configuration and an edges start")
    configs = [(root / config).resolve() for config in configs]
    for path in configs:
        if not path.is_file():
            raise FileNotFoundError(f"Configuration not found: {path}")
    execute = execute or subprocess.run
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%SZ")
    summary = root / report_root / f"{stamp}_{uuid4().hex[:8]}" / "batch_summary.json"
    batch = dict(schema="workflow_batch_v1", status="running", project_root=str(root),
                 python=sys.executable, dry_run=dry_run, series=[])
    def save():
        save_json(summary, batch)
    save()
    try:
        for index, config in enumerate(configs, 1):
            source = active_source(load_config(config),root)
            available = ("seedfc", *STAGES[STAGES.index("edges"):]) if source else STAGES
            if source and import_source:
                raise ValueError("Cannot import FC while an external source is enabled")
            stages = select_steps(available, stage=stage, start=None if stage else (start or ("seedfc" if source else "fc")))
            from srcs.workflow import validate_plan
            validate_plan(stages)
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
    except BaseException as exc:
        status = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        batch.update(status=status, reason=str(exc))
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
