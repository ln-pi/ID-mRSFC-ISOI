"""Persist statistical tables; rendering is performed only by report."""
from pathlib import Path

import pandas as pd

from srcs.analysis.inference import observations, repeated_measures, posthoc
from srcs.utils.combined_results import export_subject_metrics
from srcs.utils.dataset import statistics_metadata
from srcs.utils.provenance import sha256
from srcs.utils.results import save_json, save_table


def run(ctx):
    stat = pd.read_pickle(ctx.statistics_path)
    metadata = statistics_metadata(ctx, stat)
    axis = "week" if metadata["split_axis"] == "run" else "run"
    within = ["condition", axis] if len(stat.columns) > 1 else ["condition"]
    folder = Path(ctx.output_dir("results", "statistics"))
    tables = folder / "tables"
    export_subject_metrics(ctx, stat, metadata)
    save_json(Path(ctx.output_dir("Metadata")) / "statistics_analysis.json", metadata)
    for atlas in stat.index:
        frame = observations(stat, atlas, metadata)
        save_table(tables / f"observations_{atlas}.csv", frame, index=False)
        anova, corrected = repeated_measures(frame, within)
        save_table(tables / f"anova_{atlas}.csv", anova)
        save_table(tables / f"rm_anova_{atlas}.csv", corrected, index=False)
    excluded = set(ctx.section("analysis").get("statistics_exclude", []))
    for condition in stat.columns:
        save_table(tables / f"posthoc_{condition}.csv", posthoc(stat, condition, excluded), index=False)
    save_json(folder / "methods.json", dict(schema="combined_statistics_v1", analysis=metadata,
        input=dict(path=str(ctx.statistics_path), sha256=sha256(ctx.statistics_path)),
        within=within, inference="existing_within_mouse_combined_only", posthoc_correction="FDR_BH",
        undefined_posthoc="record_NaN_without_poisoning_finite_p_value_correction"))
    ctx.logger.info("Statistics: %s atlases; %s conditions", len(stat.index), len(stat.columns))
