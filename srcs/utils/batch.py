"""Per-record execution and durable batch summaries."""
from collections import Counter
from datetime import datetime
from pathlib import Path
from uuid import uuid4

from srcs.utils.io import save_json
from srcs.utils.progress import emit_progress


class BatchReport:
    def __init__(self, ctx, stage, total=None):
        self.logger = ctx.logger
        self.rows = []
        self.total = total
        emit_progress(0, total, 'records')
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        self.path = Path(ctx.output_dir("BatchReports")) / f"{stage}_{stamp}_{uuid4().hex[:8]}.json"
        self.persisted = False
        self.flush()

    def flush(self):
        counts = dict.fromkeys(("success", "missing", "failed"), 0)
        counts.update(Counter(row["status"] for row in self.rows))
        payload = {"counts": counts, "records": self.rows}
        self.persisted = False
        try:
            save_json(self.path, payload)
        except PermissionError as exc:
            # A viewer or antivirus may hold the current destination open on Windows.
            alternate = self.path.with_name(self.path.stem + "_recovered_" + uuid4().hex[:8] + ".json")
            self.logger.warning("Report replacement denied (%s); trying a new report: %s", exc, alternate)
            try:
                save_json(alternate, payload)
            except OSError:
                self.logger.exception("Report could not be saved; retaining records in memory and continuing. Target: %s", self.path)
                return
            self.path = alternate
        except OSError:
            self.logger.exception("Report could not be saved; retaining records in memory and continuing. Target: %s", self.path)
            return
        self.persisted = True

    def run(self, action, **identity):
        row = dict(identity)
        try:
            outcome = action()
        except FileNotFoundError as exc:
            row.update(status="missing", reason=str(exc), error_type=type(exc).__name__)
            self.logger.warning("Skipping missing record %s: %s", identity, exc)
        except Exception as exc:
            row.update(status="failed", reason=str(exc), error_type=type(exc).__name__)
            self.logger.exception("Skipping failed record %s", identity)
        else:
            status = outcome.get("status") if isinstance(outcome, dict) else outcome
            if isinstance(outcome, dict) and "input_path" in outcome:
                row["input_path"] = outcome["input_path"]
            row.update(status="reused" if status == "reused" else "success", reason="")
        self.rows.append(row)
        self.flush()
        emit_progress(len(self.rows), self.total, 'records')
        return row["status"] in ("success", "reused")

    def finish(self):
        self.flush()
        if not self.persisted:
            self.logger.error("Batch report remains UNSAVED: %s; records: %s", self.path, self.rows)
            return
        self.logger.info("Batch summary: %s; report: %s", dict(Counter(r["status"] for r in self.rows)), self.path)
