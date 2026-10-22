"""Build/verify the common complete-axis edge store without selecting a cohort."""
import json
from pathlib import Path

import pandas as pd

from srcs.analysis.connectome import EdgeStore
from srcs.utils.fc_import import verify_import
from srcs.utils.output_paths import fc_reports
from srcs.utils.provenance import code_identity, sha256
from srcs.utils.results import save_json
from srcs.utils.dataset import table_path


def store_path(ctx, atlas):
    return Path(ctx.output_path) / "results/edges" / f"{atlas}.npz"


def source_identity(ctx, atlas):
    verify_import(ctx)
    reports = fc_reports(ctx.output_path)
    if not reports:
        raise FileNotFoundError("FC audit required; run fc or explicitly --import-fc SOURCE")
    report = max(reports, key=lambda p: (p.stat().st_mtime_ns, p.name))
    rows = json.loads(report.read_text(encoding="utf-8"))["records"]
    records = {tuple(row[k] for k in ("week", "mouse", "run", "atlas")): row for row in rows}
    expected = {(w, m, r, "_".join(a)) for w in ctx.weeks for m in ctx.mice for r in ctx.runs for a in ctx.atlases}
    if len(records) != len(rows) or set(records) != expected:
        raise ValueError("FC audit scope differs from configured axes")
    sources = []
    for week in ctx.weeks:
        path = table_path(ctx, week, atlas)
        sources.append(dict(week=week, path=str(path.resolve()), sha256=sha256(path) if path.exists() else None))
    provenance = Path(ctx.output_path) / "metadata/fc_import_provenance.json"
    return dict(config_tag=ctx.analysis_tag, sources=sources,
                audit=dict(path=str(report.resolve()), sha256=sha256(report)),
                preprocessing_provenance=json.loads(provenance.read_text(encoding="utf-8")) if provenance.exists()
                else {"origin": "current_fc_pipeline", "config": ctx.config}), records


def load_store(ctx, atlas):
    """Reject stale source/config/code products instead of silently reusing them."""
    path = store_path(ctx, atlas)
    manifest = json.loads(path.with_suffix(".json").read_text(encoding="utf-8"))
    if sha256(path) != manifest["sha256"]:
        raise ValueError("Edge store bytes changed; rebuild --stage edges")
    store = EdgeStore.load(path)
    current, _ = source_identity(ctx, atlas)
    if store.metadata["input_identity"] != current:
        raise ValueError("Edge store inputs changed; rebuild --stage edges")
    identity = code_identity()
    if store.metadata["code"]["source_sha256"] != identity["source_sha256"]:
        raise ValueError("Edge store code changed; rebuild --stage edges")
    if (store.mice, store.weeks, store.runs) != (ctx.mice, ctx.weeks, ctx.runs):
        raise ValueError("Edge store labels differ from configuration")
    return store


def run(ctx):
    code = code_identity()
    for atlas in ["_".join(a) for a in ctx.atlases]:
        identity, audit = source_identity(ctx, atlas)
        tables = {row["week"]: pd.read_pickle(row["path"]) if row["sha256"] else None for row in identity["sources"]}
        store = EdgeStore.from_tables(tables, ctx.mice, ctx.weeks, ctx.runs, atlas, audit)
        store.metadata.update(input_identity=identity, code=code)
        store.save(store_path(ctx, atlas))
        save_json(store_path(ctx, atlas).with_suffix(".json"),
                  dict(metadata=store.metadata, sha256=sha256(store_path(ctx, atlas))))
        ctx.logger.info("Edges %s: %s; valid records %s/%s", atlas, store.F_r.shape,
                        store.valid_record.sum(), store.valid_record.size)
