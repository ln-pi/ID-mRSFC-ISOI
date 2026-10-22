"""P2-style sequential series/stage coordination with durable batch summaries."""
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
from uuid import uuid4

from srcs.utils.results import save_json
from srcs.workflow import STAGES


def run_batch(configs, *, project_root, stage=None, start="fc", dry_run=False,
              import_source=None, report_section="all", report_root="reports/automatic", execute=None):
    root = Path(project_root).resolve()
    if (stage is not None and stage not in STAGES) or start not in STAGES:
        raise ValueError("Unknown workflow stage")
    stages = [stage] if stage else list(STAGES[STAGES.index(start):])
    if import_source and (len(configs) != 1 or stages[0] != "edges"):
        raise ValueError("FC import requires one configuration and an edges start")
    configs = [(root / config).resolve() for config in configs]
    for path in configs:
        if not path.is_file():
            raise FileNotFoundError(f"Configuration not found: {path}")
    execute = execute or subprocess.run
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%SZ")
    summary = root / report_root / f"{stamp}_{uuid4().hex[:8]}" / "batch_summary.json"
    batch = dict(schema="workflow_batch_v1", status="running", project_root=str(root),
                 python=sys.executable, dry_run=dry_run, stages=stages, series=[])
    def save():
        save_json(summary, batch)
    save()
    try:
        for index, config in enumerate(configs, 1):
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
                entry["stages"].extend(dict(stage=s, status="planned") for s in stages[1:])
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
