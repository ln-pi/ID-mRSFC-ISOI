"""Persist per-Run tensors and the outputs of pure pairwise analysis."""
from pathlib import Path
import pandas as pd
from srcs.analysis.connectome import pairwise_tensor, save_arrays
from srcs.analysis.pairwise import selections, week_positions, analyze_views, METRICS, COLUMNS
from srcs.utils.provenance import sha256
from srcs.utils.results import save_json, save_table
from srcs.utils.edge_store import load_store, store_path


def run(ctx):
    views, positions = selections(ctx), week_positions(ctx)
    rows, audits, products = [], [], []
    for atlas in ["_".join(a) for a in ctx.atlases]:
        store = load_store(ctx, atlas)
        for run_id in store.runs:
            meta, arrays = pairwise_tensor(store, run_id)
            meta.update(edge_store=dict(path=str(store_path(ctx, atlas)), sha256=sha256(store_path(ctx, atlas))),
                        code=store.metadata["code"])
            path = Path(ctx.output_dir("results", "pairwise", atlas)) / f"{run_id}.npz"
            save_arrays(path, meta, **arrays)
            save_json(path.with_suffix(".json"), meta)
            products.append(dict(path=str(path), sha256=sha256(path), atlas=atlas, run=run_id))
            part, selection_audit = analyze_views(store, run_id, arrays, views, positions)
            rows.extend(part)
            audits.extend(selection_audit)
    folder = Path(ctx.output_dir("results", "pairwise"))
    frame = pd.DataFrame(rows, columns=COLUMNS)
    save_table(folder / "subjects_by_pair.csv", frame, index=False)
    keys = ["atlas", "run", "view", "cohort", "mouse", "scale"]
    if len(frame):
        group = frame.groupby(keys, sort=False, dropna=False)
        summary = group[METRICS].mean().join(group.size().rename("n_valid_pairs"))
        summary = summary.join(group.n_candidates.min().rename("min_candidates"))
        summary = summary.join(group.n_candidates.max().rename("max_candidates")).reset_index()
        expected = {name: len(pairs) for name, pairs in views.items()}
        summary["n_requested_pairs"] = summary.view.map(expected)
    else:
        summary = pd.DataFrame(columns=[*keys, *METRICS, "n_valid_pairs", "min_candidates", "max_candidates", "n_requested_pairs"])
    save_table(folder / "subjects_by_window.csv", summary, index=False)
    save_json(folder / "selection.json", dict(schema="pairwise_selection_v1", selections=audits,
        products=products, phases=ctx.section("pairwise").get("phases", {}), week_positions=positions,
        phase_mean="mean_of_per_pair_per_mouse_metrics; not_metrics_of_mean_similarity",
        interpretation=ctx.section("pairwise").get("window_interpretation", "observation_windows"),
        scale="pattern_Pearson_r_of_F_z", inference="descriptive_only_mouse_is_unit",
        missing="no_imputation; available pair and complete selected window reported separately"))
    ctx.logger.info("Pairwise: %s per-mouse comparisons; %s tensors", len(frame), len(products))
