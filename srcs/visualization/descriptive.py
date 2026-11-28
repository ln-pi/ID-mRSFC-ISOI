"""Amplitude and topology figures from verified saved products only."""
from math import ceil
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from srcs.utils.storage import load_arrays
from srcs.analysis.magnitude import SCHEMA
from srcs.analysis.topology import edge_matrix
from srcs.core.artifacts import artifact_entries, artifact_path
from srcs.workflow.products import verified_index
from srcs.utils.provenance import sha256
from srcs.utils.storage import save_json
from srcs.utils.visualization import displayed_weeks
from srcs.visualization.scales import resolve


def _heat(ax, matrix, labels, *, title, cmap, vmin, vmax):
    im = ax.imshow(matrix, cmap=cmap, vmin=vmin, vmax=vmax, interpolation="none", aspect="auto")
    step = max(1, ceil(len(labels)/12))
    ticks = list(range(0, len(labels), step))
    ax.set_xticks(ticks, [labels[i] for i in ticks], rotation=90, fontsize=6)
    ax.set_yticks(ticks, [labels[i] for i in ticks], fontsize=6)
    ax.set_title(title, fontsize=9)
    return im


def run(ctx):
    source = Path(ctx.output_path) / "results/descriptive/index.json"
    index = verified_index(source)
    if index["schema"] != SCHEMA:
        raise ValueError("Unsupported descriptive report schema")
    paths = [source.resolve(), *(item["path"] for key in ("products", "files")
                                 for item in artifact_entries(index[key], source.parent))]
    folder = Path(ctx.output_dir("report", "descriptive"))
    regions = pd.read_csv(source.parent / "regional.csv")
    subjects = pd.read_csv(source.parent / "subjects.csv")
    dispersion = pd.read_csv(source.parent / "dispersion.csv")
    viz = ctx.viz()
    written, examples, scales = [], [], []

    def save(fig, name):
        for suffix in ("png", "pdf"):
            target = folder / f"{name}.{suffix}"
            fig.savefig(target, dpi=ctx.figure_dpi(), bbox_inches="tight")
            written.append(target.name)
        plt.close(fig)

    def atlas_panels(vectors, titles, labels, edges, name, title, signed=False, scale=None):
        cols = min(4, len(vectors))
        rows = ceil(len(vectors)/cols)
        fig, axes = plt.subplots(rows, cols, figsize=(4*cols, 3.6*rows), squeeze=False,
                                 constrained_layout=True)
        metric = "signed_offset" if signed else "regional_variance"
        scale = scale or resolve(vectors, viz, metric, anchor_zero=not signed)
        if scale is None:
            plt.close(fig)
            return
        scales.append(scale.record(figure=name, metric=metric))
        im = None
        for k, ax in enumerate(axes.flat):
            if k >= len(vectors):
                ax.set_visible(False)
                continue
            matrix = edge_matrix(vectors[k], edges, len(labels))
            im = _heat(ax, matrix, labels, title=titles[k], **scale.image)
        fig.suptitle(title)
        fig.colorbar(im, ax=list(axes.flat), shrink=.7, extend=scale.extend,
                     label="Delta Fisher-z FC (positive=stronger)" if signed else "Edge variance (Fisher-z squared)")
        save(fig, name)

    for product in index["products"]:
        if product["cohort"] != "available_week":
            continue
        meta, arrays = load_arrays(product["path"], SCHEMA)
        atlas = meta["atlas"]
        mice, saved_weeks, views = meta["mouse_ids"], meta["week_ids"], meta["view_ids"]
        weeks = displayed_weeks(saved_weeks, viz, ctx.section("descriptive").get("reference_week"))
        positions = [saved_weeks.index(week) for week in weeks]
        region_atlas = regions[(regions.atlas == atlas) & (regions.cohort == "available_week") &
                               regions.week.isin(weeks)]
        region_scales = {
            level: resolve(part.mean_value, viz, "regional_variance", anchor_zero=True)
            for level, part in region_atlas[region_atlas.quantity == "edge_variance"].groupby("level")
        }
        edge_scale = resolve(arrays["edge_variance"][positions], viz, "regional_variance", anchor_zero=True)
        for r, view in enumerate(views):
            prefix = f"{atlas}_{view}"
            base = (dispersion.atlas == atlas) & (dispersion.view == view)
            fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), constrained_layout=True)
            for cohort, style in (("available_week", "-o"), ("fixed_window", "--s")):
                if viz.get("display_cohorts") and cohort not in viz["display_cohorts"]:
                    continue
                data = dispersion[base & (dispersion.cohort == cohort)].set_index("week").reindex(weeks)
                axes[0].plot(range(len(weeks)), data.dispersion, style, label=cohort)
                axes[1].plot(range(len(weeks)), data.n_mice, style, label=cohort)
            ss = subjects[(subjects.atlas == atlas) & (subjects.view == view) & (subjects.cohort == "available_week")]
            for mouse in mice:
                data = ss[ss.mouse == mouse].set_index("week").reindex(weeks)
                axes[2].plot(range(len(weeks)), data.rms_from_group_loo, "-o", ms=3, label=mouse)
            for ax, title in zip(axes, ("Group dispersion (denominator N)", "Valid mice", "LOO RMS per mouse | available")):
                ax.set_title(title)
                ax.set_xticks(range(len(weeks)), weeks, rotation=90)
                ax.legend(fontsize=7)
            axes[0].set_ylabel("Fisher-z squared")
            axes[2].set_ylabel("RMS Fisher-z FC")
            fig.suptitle(f"{atlas} | {view} | descriptive trajectories")
            save(fig, f"trajectories_{prefix}")

            good = np.flatnonzero(arrays["n_valid_mice"][positions, r] >= 2)
            if not len(good):
                continue
            u = positions[int(good[0])]
            example_week = saved_weeks[u]
            fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), constrained_layout=True)
            for key, ax in zip(("euclidean", "rms", "correlation"), axes):
                matrix = arrays[f"distance_{key}"][u, r]
                metric = "correlation_distance" if key == "correlation" else "distance"
                scale = resolve(matrix, viz, metric, anchor_zero=True)
                if scale is None:
                    ax.set_axis_off()
                    continue
                scales.append(scale.record(metric=metric, atlas=atlas, view=view, week=example_week))
                im = _heat(ax, matrix, mice, title=key, **scale.image)
                fig.colorbar(im, ax=ax, shrink=.8, extend=scale.extend, label="1 - pattern r" if key == "correlation" else "Fisher-z FC")
            fig.suptitle(f"{atlas} | {view} | {example_week} | available mice | distances")
            save(fig, f"distances_{prefix}")

            atlas_panels(arrays["edge_variance"][positions, r], weeks, meta["roi_labels"], arrays["edge_index"],
                         f"edge_variance_{prefix}", f"{atlas} | {view} | sample edge variance (ddof=1)", scale=edge_scale)
            # All mice at the example week are shown; the most deviant mouse is not excluded.
            delta = arrays["residual_group_loo"][:, u, r]
            scores = np.sqrt(np.mean(delta ** 2, axis=1))
            titles = [f"{mouse} | LOO RMS {score:.3g}" for mouse, score in zip(mice, scores)]
            atlas_panels(delta, titles, meta["roi_labels"], arrays["edge_index"],
                         f"delta_fc_{prefix}", f"{atlas} | {view} | {example_week} | same-week LOO Delta FC", signed=True)
            top = int(np.nanargmax(scores))
            examples.append(dict(atlas=atlas, view=view, week=example_week, mouse=mice[top],
                                 selection="first_week_with_at_least_two_mice_then_largest_LOO_RMS"))

            selected = region_atlas[region_atlas.view == view]
            fig, axes = plt.subplots(1, 3, figsize=(20, 7), constrained_layout=True)
            for ax, level in zip(axes[:2], ("roi", "network_pair")):
                data = selected[(selected.quantity == "edge_variance") & (selected.level == level)]
                if data.empty:
                    ax.text(.5, .5, "No explicit network mapping available", ha="center")
                    ax.set_axis_off()
                    continue
                matrix = data.pivot(index="region", columns="week", values="mean_value").reindex(columns=weeks)
                order = list(dict.fromkeys(data.region))
                matrix = matrix.reindex(index=order)
                scale = region_scales.get(level)
                if scale is None:
                    ax.set_axis_off()
                    continue
                scales.append(scale.record(metric="regional_variance", atlas=atlas, view=view, level=level, weeks=weeks))
                im = ax.imshow(matrix, aspect="auto", **scale.image)
                ax.set_xticks(range(len(weeks)), weeks, rotation=90, fontsize=8)
                labels = [f"{meta.get('roi_network_mapping', {}).get(roi, '')} | {roi}" for roi in matrix.index] if level == "roi" else matrix.index
                ax.set_yticks(range(len(matrix)), labels, fontsize=6)
                if level == "roi":
                    assigned = [meta.get("roi_network_mapping", {}).get(roi) for roi in matrix.index]
                    for position in range(1, len(assigned)):
                        if assigned[position] != assigned[position-1]:
                            ax.axhline(position-.5, color="white", linewidth=1)
                ax.set_title(f"{level}: mean edge variance")
                fig.colorbar(im, ax=ax, shrink=.7, extend=scale.extend)
            direction = selected[(selected.quantity == "residual_group_loo") &
                                 (selected.level == "network_pair") & (selected.week == example_week) &
                                 (selected.mouse == mice[top])]
            if len(direction):
                y = np.arange(len(direction))
                axes[2].barh(y-.2, direction.positive_rms, height=.4, color="firebrick", label="stronger (+)")
                axes[2].barh(y+.2, -direction.negative_rms, height=.4, color="steelblue", label="weaker (-)")
                axes[2].set_yticks(y, direction.region, fontsize=7)
                axes[2].axvline(0, color="gray", lw=.7)
                axes[2].legend(fontsize=8)
                axes[2].set_xlabel("Directional RMS Fisher-z FC (negative side shown left)")
                axes[2].set_title(f"Example: {mice[top]} | {example_week} | largest LOO RMS")
            else:
                axes[2].set_axis_off()
            fig.suptitle(f"{atlas} | {view} | {meta['network_kind']} | available mice")
            save(fig, f"regions_direction_{prefix}")
    save_json(folder / "index.json", dict(schema=SCHEMA, files=written, examples=examples,
        source=artifact_path(source, folder), source_sha256=sha256(source),
        note="Example is descriptive, not an outlier test; arrays/tables contain all weeks, mice and cohorts"))
    save_json(folder / "color_scales.json", scales)
    ctx.logger.info("Descriptive report: %s figures", len(written))
    return [str(p) for p in paths]
