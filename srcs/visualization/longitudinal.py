"""Saved trajectories, model estimates and explicit diagnostic views."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from srcs.analysis.longitudinal import SCHEMA
from srcs.utils.analysis_products import verified_index
from srcs.utils.provenance import sha256
from srcs.utils.results import save_json


def run(ctx):
    source = Path(ctx.output_path)/"results/longitudinal/index.json"
    index = verified_index(source)
    if index["schema"] != SCHEMA:
        raise ValueError("Unsupported longitudinal schema")
    folder = Path(ctx.output_dir("report", "longitudinal"))
    dispersion = pd.read_csv(source.parent/"dispersion.csv")
    subjects = pd.read_csv(source.parent/"subjects.csv")
    models = pd.read_csv(source.parent/"random_slopes.csv")
    slopes = pd.read_csv(source.parent/"mouse_slopes.csv")
    diagnostics = pd.read_csv(source.parent/"model_diagnostics.csv")
    written = []
    def save(fig, name):
        for suffix in ("png", "pdf"):
            fig.savefig(folder/f"{name}.{suffix}", dpi=ctx.figure_dpi(), bbox_inches="tight")
            written.append(f"{name}.{suffix}")
        plt.close(fig)
    for (atlas, view), data in dispersion.groupby(["atlas", "view"], sort=False):
        fig, axes = plt.subplots(1, 4, figsize=(20, 4.8), constrained_layout=True)
        weeks = list(dict.fromkeys(data.week))
        for cohort, part in data.groupby("cohort", sort=False):
            part = part.set_index("week").reindex(weeks)
            axes[0].plot(part.time, part.dispersion, "-o", ms=3, label=cohort)
            axes[0].fill_between(part.time.to_numpy(), part.ci_low.to_numpy(), part.ci_high.to_numpy(), alpha=.15)
            axes[1].plot(part.time, part.n_mice, "-o", ms=3, label=cohort)
            axes[3].plot(part.time, part.group_mean_fc, "-o", ms=3, label=cohort)
        sub = subjects[(subjects.atlas == atlas) & (subjects.view == view) & (subjects.cohort == "fixed_window")]
        for mouse, part in sub.groupby("mouse", sort=False):
            if not np.isfinite(part.deviation_rms).any():
                continue
            axes[2].plot(part.time, part.deviation_rms, "-o", ms=2, label=mouse)
        for ax, title in zip(axes, ("Dispersion + animal-block 95% intervals", "Valid mice", "Individual LOO RMS | fixed window", "Group mean FC")):
            ax.set_title(title, fontsize=10)
            ax.set_xticks(data.time.unique(), weeks, rotation=90)
            ax.legend(fontsize=6)
        axes[0].set_ylabel("Fisher-z squared")
        axes[2].set_ylabel("RMS Fisher-z FC")
        axes[3].set_ylabel("Mean Fisher-z FC across edges")
        fig.suptitle(f"{atlas} | {view} | trajectories; fixed and available cohorts kept separate")
        save(fig, f"trajectories_{atlas}_{view}")
    for atlas, data in models.groupby("atlas", sort=False):
        fig, axes = plt.subplots(1, 3, figsize=(17, 6), constrained_layout=True)
        labels = [f"{row.target} | {row.status}" + (" | Hessian warning" if getattr(row, "hessian_warning", False) is True else "")
                  for row in data.itertuples()]
        estimates = data.fixed_slope if "fixed_slope" in data else np.full(len(data), np.nan)
        axes[0].barh(range(len(data)), estimates)
        axes[0].set_yticks(range(len(data)), labels, fontsize=8)
        axes[0].set(xlabel="Fixed FC slope per numeric time unit", title="Network REML model status")
        axes[0].axvline(0, color="gray", lw=.7)
        if "target" in slopes:
            part = slopes[(slopes.atlas == atlas) & (slopes.target == "all_edges")]
            axes[1].barh(range(len(part)), part.total_slope)
            axes[1].set_yticks(range(len(part)), part.mouse)
            axes[1].set(xlabel="Fixed + conditional mouse slope", title="All-edge FC mean | mouse slopes")
        if "target" in diagnostics:
            part = diagnostics[(diagnostics.atlas == atlas) & (diagnostics.target == "all_edges")]
            for run_id, group in part.groupby("run"):
                axes[2].scatter(group.time, group.residual, s=12, alpha=.5, label=run_id)
            axes[2].axhline(0, color="gray", lw=.7)
            axes[2].set(xlabel="Numeric time", ylabel="Conditional residual", title="Check temporal structure / heteroscedasticity")
            axes[2].legend(fontsize=8)
        fig.suptitle(f"{atlas} | linear-time mouse random intercept/slope + mouse-week effect | no causal interpretation")
        save(fig, f"random_slopes_{atlas}")
    save_json(folder/"index.json", dict(schema=SCHEMA, files=written, source=str(source), source_sha256=sha256(source)))
    return [str(source)]+[item["path"] for item in index["files"]+index["products"]]
