"""Canonical week x mouse x run data and axis-independent analysis views."""
import json
from itertools import permutations
from pathlib import Path

import pandas as pd


def recording_folders(root, week, mouse, run):
    return list(dict.fromkeys(Path(root).joinpath(*parts)
                             for parts in permutations((week, mouse, run))))


def resolve_recording(ctx, week, mouse, run, raw=False):
    root = Path(ctx.raw_root) / ctx.dataset if raw else Path(ctx.processed_path)
    folders = recording_folders(root, week, mouse, run)
    if raw:
        matches = [p for p in folders if (p / "info.txt").is_file() or any(p.glob("img*.bin"))]
    else:
        name = f"{mouse}_{run}_{week}_FiltDRR.npy"
        matches = [p / name for p in folders if (p / name).is_file()]
    matches = list(dict.fromkeys(p.resolve() for p in matches))
    if len(matches) > 1:
        raise ValueError(f"Ambiguous recording {(week, mouse, run)}: {matches}")
    if not matches:
        raise FileNotFoundError(f"Missing recording {(week, mouse, run)} under {root}; "
                                + ("expected info.txt/img*.bin" if raw else f"expected {mouse}_{run}_{week}_FiltDRR.npy"))
    return matches[0]


def preprocessing_folder(ctx, week, mouse, run):
    folders = [p for p in recording_folders(ctx.processed_path, week, mouse, run) if p.is_dir()]
    if len(folders) > 1:
        raise ValueError(f"Multiple processed recording folders for {(week, mouse, run)}: {folders}")
    return folders[0] if folders else Path(ctx.processed_dir(week, mouse, run))


def table_path(ctx, week, atlas):
    return Path(ctx.correlation_path) / f"{week}_correlation_matrices_{atlas}_excluded.pkl"


def load_week_tables(ctx, atlas, mice, weeks, runs):
    return {w: pd.read_pickle(table_path(ctx, w, atlas)).loc[mice, runs] for w in weeks}


def axis_frames(tables, split_axis, mice, weeks, runs):
    """Each frame has mouse rows and the opposite axis as columns."""
    if split_axis == "week":
        frames = {w: tables[w].loc[mice, runs].copy() for w in weeks}
    else:
        frames = {}
        for r in runs:
            frame = pd.DataFrame(index=mice, columns=weeks, dtype=object)
            for w in weeks:
                for m in mice:
                    frame.at[m, w] = tables[w].at[m, r]
            frames[r] = frame
    for name, frame in frames.items():
        frame.attrs["axis_id"] = name
    return frames


def select_complete_cohort(ctx):
    from srcs.comp.correlation import valid_fc
    mode = ctx.section("ID").get("completeness_mode", "strict")
    if mode not in ("strict", "cover"):
        raise ValueError("ID.completeness_mode must be strict or cover")
    weeks, mice, runs = ctx.weeks, ctx.mice, ctx.runs
    if not ctx.atlases:
        raise ValueError("At least one atlas is required")
    # Validate public options before reading data; cover always operates on weeks.
    _ = ctx.split_axis, ctx.selection_mode
    reports = list((Path(ctx.output_path) / "BatchReports").glob("02_correlation_*.json"))
    if not reports:
        raise FileNotFoundError("03 requires a current-config 02_correlation_*.json BatchReport")
    report_path = max(reports, key=lambda p: (p.stat().st_mtime_ns, p.name))
    records = {}
    for row in json.loads(report_path.read_text(encoding="utf-8"))["records"]:
        key = tuple(row[k] for k in ("week", "mouse", "run", "atlas"))
        if key in records:
            raise ValueError(f"Duplicate record in BatchReport: {key}")
        records[key] = row
    tags = ["_".join(a) for a in ctx.atlases]
    expected = {(w, m, r, a) for w in weeks for m in mice for r in runs for a in tags}
    if set(records) != expected:
        raise ValueError(f"Latest 02 BatchReport is incomplete or has different scope: {report_path}; "
                         f"missing={len(expected-set(records))}, extra={len(set(records)-expected)}")
    reasons, paths = {m: [] for m in mice}, []
    for atlas in tags:
        reference = None
        for week in weeks:
            path = table_path(ctx, week, atlas)
            frame = pd.read_pickle(path)  # Missing whole tables are errors, never cohort exclusions.
            if not isinstance(frame, pd.DataFrame) or not frame.index.is_unique or not frame.columns.is_unique:
                raise ValueError(f"Invalid aggregate table: {path}")
            paths.append(str(path))
            for mouse in mice:
                for run in runs:
                    row = records[(week, mouse, run, atlas)]
                    issues = []
                    if row["status"] not in ("success", "reused"):
                        issues.append(f"report {row['status']}: {row.get('reason', '')}")
                    cell = frame.at[mouse, run] if mouse in frame.index and run in frame.columns else None
                    if not isinstance(cell, pd.DataFrame):
                        issues.append("missing matrix in PKL")
                    elif not valid_fc(cell, list(cell.index)):
                        issues.append("invalid matrix in PKL")
                    elif reference is None:
                        reference = list(cell.index)
                    elif list(cell.index) != reference:
                        issues.append("inconsistent ROI labels/order")
                    if issues:
                        reasons[mouse].append(dict(week=week, run=run, atlas=atlas, reasons=issues))
    bad = {m: {r["week"] for r in reasons[m]} for m in mice}
    spans = [(0, len(weeks))] if mode == "strict" else [
        (start, stop) for start in range(len(weeks)) for stop in range(start + 1, len(weeks) + 1)]
    candidates = []
    for start, stop in spans:
        selected = weeks[start:stop]
        kept = [m for m in mice if not bad[m].intersection(selected)]
        candidates.append(dict(start_index=start, stop_index_exclusive=stop, weeks=selected,
                               mice=kept, score=len(kept)*len(selected), eligible=len(kept)>=2))
    eligible = [c for c in candidates if c["eligible"]]
    chosen = max(eligible, key=lambda c: (c["score"], len(c["mice"]), -c["start_index"])) if eligible else None
    kept, selected = (chosen["mice"], chosen["weeks"]) if chosen else ([], [])
    audit = dict(policy=mode, schema="week_mouse_run_v1", batch_report=str(report_path), tables=paths,
                 configured_mice=mice, retained_mice=kept, configured_weeks=weeks, retained_weeks=selected,
                 excluded_weeks=[w for w in weeks if w not in selected], runs=runs, atlases=tags,
                 excluded_mice={m: reasons[m] for m in mice if m not in kept}, invalid_records=reasons,
                 selected_score=chosen["score"] if chosen else 0, candidates=candidates,
                 tie_break="score, mouse count, earliest configured start",
                 status="complete" if chosen else "insufficient_complete_mice")
    return kept, audit


def analysis_plan(ctx, weeks):
    ids = weeks if ctx.split_axis == "run" else ctx.runs
    axis = "week" if ctx.split_axis == "run" else "run"
    if ctx.selection_mode == "single":
        return {v: dict(weeks=[v] if axis == "week" else weeks,
                        runs=[v] if axis == "run" else ctx.runs) for v in ids}
    suffix = "Week" if axis == "week" else "Run"
    return {f"{n}{suffix}": dict(weeks=ids[:n] if axis == "week" else weeks,
                                 runs=ids[:n] if axis == "run" else ctx.runs)
            for n in range(1, len(ids)+1)}


def statistics_metadata(ctx, stat):
    meta = stat.attrs.get("analysis")
    if not meta or meta.get("split_axis") != ctx.split_axis or meta.get("selection_mode") != ctx.selection_mode:
        raise ValueError("Statistics aggregate does not match the current analysis axes")
    if list(meta["conditions"]) != list(stat.columns):
        raise ValueError("Statistics condition labels do not match analysis metadata")
    return meta
