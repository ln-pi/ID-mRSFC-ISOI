"""Independently redraw saved pairwise numerical products on shared r scales."""
from math import ceil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from srcs.utils.storage import load_arrays
from srcs.analysis.connectome import SIM_SCHEMA
from srcs.core.artifacts import artifact_path
from srcs.workflow.products import verified_index
from srcs.utils.provenance import sha256
from srcs.utils.storage import save_json
from srcs.utils.visualization import displayed_weeks
from srcs.visualization.scales import resolve


def display_options(ctx):
    """Compatibility for callers inspecting the historical pairwise plot settings."""
    viz = ctx.viz()
    legacy = ctx.config.get("visualization", {}).get("pairwise", {})
    limits = (legacy.get("vmin", -1.), legacy.get("vmax", 1.))
    scale = resolve(limits, viz, "pairwise", fixed=limits)
    if not -1 <= scale.vmin < scale.vmax <= 1:
        raise ValueError("visualization.pairwise requires -1 <= vmin < vmax <= 1")
    return scale.vmin, scale.vmax, scale.cmap


def range_audit(values, mask, low, high):
    finite = np.asarray(values)[np.asarray(mask, bool) & np.isfinite(values)]
    under, over = int(np.sum(finite < low)), int(np.sum(finite > high))
    return dict(finite_values=len(finite), negative_values=int(np.sum(finite < 0)),
                minimum=float(finite.min()) if finite.size else None,
                maximum=float(finite.max()) if finite.size else None,
                below_scale=under, above_scale=over,
                extend="both" if under and over else "min" if under else "max" if over else "neither")


def run(ctx):
    source = Path(ctx.output_dir("results","pairwise",create=False))/"selection.json"
    selection = verified_index(source)
    products = {(p["atlas"], p["run"]): p for p in selection["products"]}
    folder = Path(ctx.output_dir("report", "pairwise"))
    written = []
    viz = ctx.viz()
    reference = ctx.section("descriptive").get("reference_week")
    legacy = ctx.config.get("visualization", {}).get("pairwise", {})
    fixed = display_options(ctx)[:2] if {"vmin", "vmax"} & set(legacy) and "pairwise" not in ctx.config.get("visualization", {}).get("metrics", {}) else None
    ranges = []
    prepared = {}
    populations = {}
    for (atlas, run_id), product in products.items():
        path = Path(product["path"])
        meta, arrays = load_arrays(path, SIM_SCHEMA)
        saved_weeks = meta["week_ids"]
        weeks = displayed_weeks(saved_weeks, viz, reference)
        positions = [saved_weeks.index(week) for week in weeks]
        values = np.where(arrays["statistical_mask"], arrays["S_r"], np.nan)
        within = values[np.arange(len(meta["mouse_ids"])), np.arange(len(meta["mouse_ids"]))][:, positions][:, :, positions]
        views = []
        for view in selection["selections"]:
            if view["atlas"] != atlas or view["run"] != run_id or view["view"] == "all_weeks":
                continue
            pairs = [pair for pair in view["pairs"] if pair["week_u"] in weeks and pair["week_v"] in weeks]
            matrices = [values[:, :, saved_weeks.index(pair["week_u"]), saved_weeks.index(pair["week_v"])] for pair in pairs]
            if pairs:
                views.append((view["view"], pairs, matrices))
        populations.setdefault(atlas, []).extend([within.ravel(), *(m.ravel() for _, _, group in views for m in group)])
        prepared[atlas, run_id] = meta["mouse_ids"], weeks, within, views
    scales = {atlas: resolve(np.concatenate(values), viz, "pairwise", bounds=(-1, 1), fixed=fixed)
              for atlas, values in populations.items()}
    def save(fig, name):
        fig.savefig(folder / f"{name}.png", dpi=ctx.figure_dpi(), bbox_inches="tight")
        fig.savefig(folder / f"{name}.pdf", bbox_inches="tight")
        plt.close(fig)
        written.extend([f"{name}.png", f"{name}.pdf"])
    for (atlas, run_id), (mice, weeks, within, views) in prepared.items():
        scale = scales[atlas]
        if scale is None:
            continue
        low, high = scale.vmin, scale.vmax
        audit = scale.record(atlas=atlas, run=run_id, weeks=weeks, below_scale=scale.under, above_scale=scale.over)
        ranges.append(audit)
        if audit["below_scale"] or audit["above_scale"]:
            ctx.logger.warning("Pairwise %s/%s: %s below and %s above display range [%s, %s]; "
                               "out-of-range values retain distinct colors and colorbar extensions",
                               atlas, run_id, audit["below_scale"], audit["above_scale"], low, high)
        cols = min(4, len(mice))
        fig, axes = plt.subplots(ceil(len(mice)/cols), cols, figsize=(4*cols, 3.6*ceil(len(mice)/cols)), squeeze=False,
                                 constrained_layout=True)
        for i, ax in enumerate(axes.flat):
            if i >= len(mice):
                ax.set_visible(False)
                continue
            image = ax.imshow(within[i], **scale.image, interpolation="none")
            ax.set_xticks(range(len(weeks)), weeks, rotation=90, fontsize=7)
            ax.set_yticks(range(len(weeks)), weeks, fontsize=7)
            ax.set_title(mice[i])
            ax.set(xlabel="week v", ylabel="week u")
        fig.suptitle(f"{atlas} | {run_id} | within-mouse Pearson r; same-record diagonal masked")
        fig.colorbar(image, ax=axes.ravel().tolist(), label="Pattern Pearson r", shrink=.7, extend=scale.extend)
        save(fig, f"stability_{atlas}_{run_id}")
        for view, pairs, matrices in views:
            cols = min(3, len(pairs))
            fig, axes = plt.subplots(ceil(len(pairs)/cols), cols, figsize=(4*cols, 3.8*ceil(len(pairs)/cols)),
                                     squeeze=False, constrained_layout=True)
            for k, ax in enumerate(axes.flat):
                if k >= len(pairs):
                    ax.set_visible(False)
                    continue
                pair = pairs[k]
                image = ax.imshow(matrices[k], **scale.image, interpolation="none")
                ax.set_xticks(range(len(mice)), mice, rotation=90, fontsize=7)
                ax.set_yticks(range(len(mice)), mice, fontsize=7)
                ax.set(xlabel=f"week v: {pair['week_v']}", ylabel=f"week u: {pair['week_u']}")
                ax.set_title(f"{pair['week_u']} to {pair['week_v']} | gap {pair['week_gap']} weeks")
            fig.suptitle(f"{atlas} | {run_id} | {view}")
            fig.colorbar(image, ax=axes.ravel().tolist(), label="Pattern Pearson r", shrink=.7, extend=scale.extend)
            save(fig, f"{view}_{atlas}_{run_id}")
    from srcs.visualization.phase_pairwise import render
    phase_files = render([(ctx.variables.get("series", "series"), p) for p in selection.get("phase_products", [])],
                         folder, ctx.figure_dpi(), viz=viz, reference_week=reference)
    written.extend(Path(p).name for p in phase_files)
    save_json(folder / "index.json", dict(files=written, source=artifact_path(source, folder), source_sha256=sha256(source),
        scale="pattern_Pearson_r_of_F_z", cmap=viz["colormaps"][viz["metrics"]["pairwise"]["color_mode"]],
        out_of_range="distinct under/over colors; values unchanged", data_ranges=ranges,
        nan="masked", same_record="excluded",
        layout="standard_directed_panels; rows=mouse_i_at_u, columns=mouse_j_at_v"))
    save_json(folder / "color_scales.json", ranges)
    ctx.logger.info("Report: %s saved figures", len(written))
