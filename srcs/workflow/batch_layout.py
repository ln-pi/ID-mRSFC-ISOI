"""Allocate one readable output directory per workflow invocation."""
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
import re

from srcs.core.artifacts import write_json


def component(value, fallback):
    text = re.sub(r"[^A-Za-z0-9._-]+", "-", str(value or "")).strip("._-")
    return (text or fallback)[:60].rstrip("._-")


def series_label(scope):
    return "-".join(filter(None, (scope["series"], scope["branch"])))


def dataset_label(scopes):
    scopes = list(scopes)
    datasets = list(dict.fromkeys(s.get("dataset") for s in scopes))
    series = list(dict.fromkeys(s["series"] for s in scopes))
    branches = list(dict.fromkeys(s["branch"] for s in scopes))
    groups = "AN" if set(series) == {"Aseries", "Nseries"} else "-".join(series)
    return "-".join(filter(None, [*(s for s in datasets if s), groups, *branches]))


@dataclass(frozen=True)
class BatchLayout:
    root: Path

    @staticmethod
    def candidates(output_root, data_label, task, *, now=None):
        output_root = Path(output_root).resolve()
        stamp = (now or datetime.now(timezone.utc)).strftime("%Y%m%d-%H%M%SZ")
        name = f"Batch_{stamp}_{component(data_label, 'Series')}_{component(task, 'Core')}"
        yield output_root / name
        sequence = 2
        while True:
            yield output_root / f"{name}_{sequence:02d}"
            sequence += 1

    @classmethod
    def preview(cls, output_root, data_label, task, *, now=None):
        return cls(next(path for path in cls.candidates(output_root, data_label, task, now=now)
                        if not path.exists()))

    @classmethod
    def create(cls, output_root, data_label, task, *, now=None):
        for path in cls.candidates(output_root, data_label, task, now=now):
            try:
                path.mkdir(parents=True)
            except FileExistsError:
                continue
            return cls(path)

    @property
    def summary_path(self):
        return self.root / "batch_summary.json"

    @property
    def compare(self):
        return self.root / "compare"

    def relative(self, path):
        return Path(path).resolve().relative_to(self.root.resolve()).as_posix()

    def save_summary(self, data):
        write_json(self.summary_path, {"schema": "workflow_batch_v2", **data})

    def logger_path(self, name):
        return self.root / "logs" / f"{component(name, 'batch')}.log"


def series_input(project_root, value):
    """Accept an explicit series path, a path relative to output, or a legacy basename."""
    project = Path(project_root).resolve()
    path = Path(value)
    output = project / "data" / "output"
    candidates = [(project / path).resolve(), (output / path).resolve()]
    bare = not path.is_absolute() and len(path.parts) == 1
    if bare:
        candidates.extend((batch / path.name).resolve() for batch in sorted(output.glob("Batch_*")) if batch.is_dir())
    matches = [candidate for candidate in dict.fromkeys(candidates)
               if (candidate / "metadata" / "analysis_config.json").is_file()]
    if len(matches) == 1:
        return matches[0]
    if matches:
        raise ValueError(f"Multiple batches or saved outputs contain {path.name}; specify the full series path")
    raise FileNotFoundError(f"No saved series found: {value}")


def saved_scopes(roots, configs):
    import json
    for root, config in zip(roots, configs):
        identity = root / "metadata" / "data_identity.json"
        if identity.is_file():
            yield json.loads(identity.read_text(encoding="utf-8"))["scope"]
        else:
            variables = config["variables"]
            seedfc = config.get("plugins", {}).get("source", {}).get("seedfc", {})
            branch = seedfc.get("options", {}).get("branch") if seedfc.get("enabled") else None
            yield {"dataset": variables.get("dataset"), "series": variables.get("series") or root.name,
                   "branch": branch or ("GSR" if config.get("preprocess", {}).get("gsr", True) else "noGSR")}
