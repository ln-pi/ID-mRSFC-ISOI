"""Independently redraw saved pairwise numerical products on shared r scales."""
import json
from math import ceil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from srcs.analysis.connectome import load_arrays, SIM_SCHEMA
from srcs.utils.provenance import sha256
from srcs.utils.results import save_json


def display_options(ctx):
    cfg = ctx.section("visualization").get("pairwise", {})
    low, high = cfg.get("vmin", -1.), cfg.get("vmax", 1.)
    if (isinstance(low, bool) or isinstance(high, bool)
            or not isinstance(low, (int, float)) or not isinstance(high, (int, float))
            or not np.isfinite([low, high]).all() or not -1 <= low < high <= 1):
        raise ValueError("visualization.pairwise requires -1 <= vmin < vmax <= 1")
    cmap = plt.get_cmap(cfg.get("cmap", "inferno")).with_extremes(
        under=cfg.get("under_color", "#d73027"), over=cfg.get("over_color", "#762a83"))
    return float(low), float(high), cmap


def range_audit(values, mask, low, high):
    finite = np.asarray(values)[np.asarray(mask, bool) & np.isfinite(values)]
    under, over = int(np.sum(finite < low)), int(np.sum(finite > high))
    return dict(finite_values=len(finite), negative_values=int(np.sum(finite < 0)),
                minimum=float(finite.min()) if finite.size else None,
                maximum=float(finite.max()) if finite.size else None,
                below_scale=under, above_scale=over,
                extend="both" if under and over else "min" if under else "max" if over else "neither")


def run(ctx):
    source = Path(ctx.output_path) / "results/pairwise/selection.json"
    selection = json.loads(source.read_text(encoding="utf-8"))
    products = {(p["atlas"], p["run"]): p for p in selection["products"]}
    folder = Path(ctx.output_dir("report", "pairwise"))
    written = []
    low, high, cmap = display_options(ctx)
    ranges = []
    def save(fig, name):
        fig.savefig(folder / f"{name}.png", dpi=ctx.figure_dpi(), bbox_inches="tight")
        fig.savefig(folder / f"{name}.pdf", bbox_inches="tight")
        plt.close(fig)
        written.extend([f"{name}.png", f"{name}.pdf"])
    for (atlas, run_id), product in products.items():
        path = Path(product["path"])
        if sha256(path) != product["sha256"]:
            raise ValueError(f"Saved pairwise tensor changed: {path}; rerun pairwise")
        meta, arrays = load_arrays(path, SIM_SCHEMA)
        audit = dict(atlas=atlas, run=run_id, **range_audit(arrays["S_r"], arrays["statistical_mask"], low, high))
        ranges.append(audit)
        if audit["below_scale"] or audit["above_scale"]:
            ctx.logger.warning("Pairwise %s/%s: %s below and %s above display range [%s, %s]; "
                               "out-of-range values retain distinct colors and colorbar extensions",
                               atlas, run_id, audit["below_scale"], audit["above_scale"], low, high)
        mice, weeks = meta["mouse_ids"], meta["week_ids"]
        cols = min(4, len(mice))
        fig, axes = plt.subplots(ceil(len(mice)/cols), cols, figsize=(4*cols, 3.6*ceil(len(mice)/cols)), squeeze=False,
                                 constrained_layout=True)
        for i, ax in enumerate(axes.flat):
            if i >= len(mice):
                ax.set_visible(False)
                continue
            matrix = arrays["S_r"][i, i].copy()
            matrix[~arrays["statistical_mask"][i, i]] = np.nan
            image = ax.imshow(matrix, vmin=low, vmax=high, cmap=cmap, interpolation="none")
            ax.set_xticks(range(len(weeks)), weeks, rotation=90, fontsize=7)
            ax.set_yticks(range(len(weeks)), weeks, fontsize=7)
            ax.set_title(mice[i])
            ax.set(xlabel="week v", ylabel="week u")
        fig.suptitle(f"{atlas} | {run_id} | within-mouse Pearson r; same-record diagonal masked")
        fig.colorbar(image, ax=axes.ravel().tolist(), label="Pattern Pearson r", shrink=.7, extend=audit["extend"])
        save(fig, f"stability_{atlas}_{run_id}")
        for view in selection["selections"]:
            if view["atlas"] != atlas or view["run"] != run_id or view["view"] == "all_weeks" or not view["pairs"]:
                continue
            pairs = view["pairs"]
            cols = min(3, len(pairs))
            fig, axes = plt.subplots(ceil(len(pairs)/cols), cols, figsize=(4*cols, 3.8*ceil(len(pairs)/cols)),
                                     squeeze=False, constrained_layout=True)
            for k, ax in enumerate(axes.flat):
                if k >= len(pairs):
                    ax.set_visible(False)
                    continue
                pair = pairs[k]
                u, v = weeks.index(pair["week_u"]), weeks.index(pair["week_v"])
                image = ax.imshow(arrays["S_r"][:, :, u, v], vmin=low, vmax=high, cmap=cmap, interpolation="none")
                ax.set_xticks(range(len(mice)), mice, rotation=90, fontsize=7)
                ax.set_yticks(range(len(mice)), mice, fontsize=7)
                ax.set(xlabel=f"week v: {pair['week_v']}", ylabel=f"week u: {pair['week_u']}")
                ax.set_title(f"{pair['week_u']} to {pair['week_v']} | gap {pair['week_gap']} weeks")
            fig.suptitle(f"{atlas} | {run_id} | {view['view']}")
            fig.colorbar(image, ax=axes.ravel().tolist(), label="Pattern Pearson r", shrink=.7, extend=audit["extend"])
            save(fig, f"{view['view']}_{atlas}_{run_id}")
    save_json(folder / "index.json", dict(files=written, source=str(source), source_sha256=sha256(source),
        scale="pattern_Pearson_r_of_F_z", limits=[low, high], cmap=cmap.name,
        out_of_range="distinct under/over colors; values unchanged", data_ranges=ranges,
        nan="masked", same_record="excluded",
        layout="standard_directed_panels; rows=mouse_i_at_u, columns=mouse_j_at_v"))
    ctx.logger.info("Report: %s saved figures", len(written))
