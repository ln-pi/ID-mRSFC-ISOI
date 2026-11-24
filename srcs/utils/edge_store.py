"""Build/verify the common complete-axis edge store without selecting a cohort."""
import json
from pathlib import Path

import pandas as pd

from dataclasses import dataclass
import numpy as np
from srcs.analysis.connectome import EDGE_CLIP, SCHEMA
from srcs.utils.storage import save_arrays, load_arrays
from srcs.utils.fc_import import verify_import
from srcs.utils.layout import fc_reports
from srcs.utils.provenance import code_identity, sha256
from srcs.utils.storage import save_json
from srcs.utils.dataset import table_path


@dataclass
class EdgeStore:
    metadata: dict
    F_r: np.ndarray
    F_z: np.ndarray
    valid_record: np.ndarray
    invalid_reason: np.ndarray
    edge_index: np.ndarray

    @property
    def mice(self):
        return self.metadata["mouse_ids"]

    @property
    def weeks(self):
        return self.metadata["week_ids"]

    @property
    def runs(self):
        return self.metadata["run_ids"]

    @property
    def roi_labels(self):
        return self.metadata["roi_labels"]

    @classmethod
    def from_tables(cls, tables, mice, weeks, runs, atlas, record_audit=None):
        # Fix the ROI/edge basis using the first valid labeled FC in config order.
        from srcs.analysis.connectivity import valid_fc
        labels = None
        for week in weeks:
            table = tables.get(week)
            if table is None:
                continue
            if (not isinstance(table, pd.DataFrame) or list(table.index) != list(mice)
                    or list(table.columns) != list(runs) or not table.index.is_unique
                    or not table.columns.is_unique):
                raise ValueError(f"FC table axes differ from configuration: {week}")
            for mouse in mice:
                for run in runs:
                    cell = table.at[mouse, run]
                    if (isinstance(cell, pd.DataFrame) and cell.index.is_unique
                            and len(cell) >= 3 and valid_fc(cell, list(cell.index))):
                        if labels is None:
                            labels = list(cell.index)
        if labels is None:
            raise ValueError(f"No valid labeled FC matrices for {atlas}")
        edge_index = np.asarray(np.triu_indices(len(labels), k=1))
        shape = (len(mice), len(weeks), len(runs))
        raw = np.full((*shape, edge_index.shape[1]), np.nan)
        valid = np.zeros(shape, bool)
        reasons = np.full(shape, "missing_record", dtype="U256")
        for u, week in enumerate(weeks):
            table = tables.get(week)
            if table is None:
                reasons[:, u, :] = "missing_week_table"
                continue
            for i, mouse in enumerate(mice):
                for r, run in enumerate(runs):
                    cell = table.at[mouse, run]
                    if record_audit is not None:
                        row = record_audit[(week, mouse, run, atlas)]
                        if row["status"] not in ("success", "reused"):
                            reasons[i, u, r] = f"audit_{row['status']}: {row.get('reason', '')}"
                            continue
                    if cell is None or (pd.api.types.is_scalar(cell) and pd.isna(cell)):
                        continue
                    if not isinstance(cell, pd.DataFrame) or not valid_fc(cell, labels):
                        reasons[i, u, r] = "invalid_fc_or_roi_order"
                        continue
                    raw[i, u, r] = cell.to_numpy(dtype=float)[tuple(edge_index)]
                    valid[i, u, r], reasons[i, u, r] = True, ""
        meta = dict(schema=SCHEMA, atlas=atlas, axes=["mouse", "week", "run", "edge"],
                    shape=list(raw.shape), mouse_ids=list(mice), week_ids=list(weeks),
                    run_ids=list(runs), roi_labels=labels, edge_order="numpy.triu_indices(k=1)",
                    edge_transform=dict(method="atanh", clip=[-EDGE_CLIP, EDGE_CLIP]),
                    missing_policy="retain_configured_axes_nan_no_imputation",
                    algorithm="fc-edges-1")
        return cls(meta, raw, np.arctanh(np.clip(raw, -EDGE_CLIP, EDGE_CLIP)),
                   valid, reasons, edge_index)

    def save(self, path):
        save_arrays(path, self.metadata, F_r=self.F_r, F_z=self.F_z,
                    valid_record=self.valid_record, invalid_reason=self.invalid_reason,
                    edge_index=self.edge_index)

    @classmethod
    def load(cls, path):
        meta, arrays = load_arrays(path, SCHEMA)
        obj = cls(meta, **arrays)
        shape = (len(obj.mice), len(obj.weeks), len(obj.runs), len(obj.edge_index[0]))
        if (obj.F_r.shape != shape or obj.F_z.shape != shape
                or obj.valid_record.shape != shape[:3] or obj.invalid_reason.shape != shape[:3]
                or not np.array_equal(obj.edge_index, np.asarray(np.triu_indices(len(obj.roi_labels), 1)))):
            raise ValueError("Edge store dimensions/labels are inconsistent")
        return obj

    def subset(self, mice, weeks, runs, scale="z"):
        indices = [self.mice.index(m) for m in mice], [self.weeks.index(w) for w in weeks], [self.runs.index(r) for r in runs]
        if scale not in ("r", "z"):
            raise ValueError("Edge scale must be r or z")
        values = self.F_z if scale == "z" else self.F_r
        return values[np.ix_(*indices, range(values.shape[-1]))]

    def average_halves(self, mice, weeks, runs, split_axis, groups):
        """Two mouse x edge mean-F_z arrays; never reconstruct nested matrices."""
        from srcs.analysis.selection import validate_groups
        if split_axis not in ("week", "run"):
            raise ValueError("split_axis must be week or run")
        values = self.subset(mice, weeks, runs)
        if not np.isfinite(values).all():
            raise ValueError("Combined selection requires complete FC records")
        axis = 1 if split_axis == "week" else 2
        validate_groups(groups, values.shape[axis])
        return tuple(np.take(values, group, axis=axis).mean(axis=(1, 2)) for group in groups)


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
                else {"origin": "current_fc_pipeline", "data_identity": ctx.data_identity}), records


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
    if store.metadata.get("numerical_code") != ctx.data_identity["stages"]["edges"]["code"]:
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
        store.metadata.update(input_identity=identity, code=code,
                              numerical_code=ctx.data_identity["stages"]["edges"]["code"])
        store.save(store_path(ctx, atlas))
        save_json(store_path(ctx, atlas).with_suffix(".json"),
                  dict(metadata=store.metadata, sha256=sha256(store_path(ctx, atlas))))
        ctx.logger.info("Edges %s: %s; valid records %s/%s", atlas, store.F_r.shape,
                        store.valid_record.sum(), store.valid_record.size)
