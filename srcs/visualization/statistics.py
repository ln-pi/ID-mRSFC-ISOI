"""Render saved combined subject metrics without recalculating inference."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from srcs.utils.dataset import statistics_metadata
from srcs.visualization.combined import plot


def run(ctx):
    stat = pd.read_pickle(ctx.statistics_path)
    metadata = statistics_metadata(ctx, stat)
    path = Path(ctx.output_path) / "results/statistics/individual_identifiability.csv"
    frame = pd.read_csv(path)
    enabled = ctx.render_options()["save"] or ctx.render_options()["show"]
    axis = "week" if metadata["split_axis"] == "run" else "run"
    for atlas in stat.index:
        if not enabled:
            continue
        style = ctx.viz_style("stat_mean_error")
        style.update(title=atlas, xlabel=("Cumulative " if metadata["selection_mode"] == "cumulative" else "") + axis.capitalize(),
                     xticklabels=list(stat.columns))
        plot(ctx, ([cell[0] for cell in stat.loc[atlas]], [cell[1] for cell in stat.loc[atlas]]), "mean_error",
             f"Line_Group-Indiv_{atlas}_gradient", Path(ctx.output_dir("Statistics")), style)
        for metric in ("self_match_advantage", "self_match_percentile"):
            selected = frame[frame.atlas == atlas]
            if selected.empty:
                continue  # An explicitly excluded atlas has no subject metrics.
            table = selected.pivot(index="mouse", columns="analysis_condition", values=metric)
            table = table.reindex(index=metadata["mice"], columns=[str(c) for c in stat.columns])
            fig, ax = plt.subplots(figsize=(max(6, len(stat.columns)*.65), 4.5))
            x = np.arange(len(table.columns))
            for mouse, values in table.iterrows():
                ax.plot(x, values, marker="o", ms=3, lw=1, alpha=.6, label=mouse)
            ax.set_xticks(x, table.columns, rotation=45 if len(x) > 5 else 0)
            ax.set(xlabel="Analysis condition", ylabel=metric.replace("_", " "), title=atlas)
            ax.axhline(0 if metric == "self_match_advantage" else .5, color="gray", ls=":", lw=1)
            if metric == "self_match_percentile":
                ax.set_ylim(-.025, 1.025)
            ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
            fig.tight_layout()
            if ctx.render_options()["save"]:
                folder = Path(ctx.output_dir("Statistics", "IndividualIdentifiability"))
                fig.savefig(folder / f"{metric}_{atlas}.png", dpi=ctx.figure_dpi())
            if ctx.render_options()["show"]:
                plt.show()
            plt.close(fig)
    return [str(ctx.statistics_path), str(path)]
