"""Shared-coordinate comparisons of saved series products."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from srcs.core.artifacts import sha256, write_json
from srcs.visualization.phase_pairwise import render as render_phases
from srcs.analysis.phases import weekly_metric_means
from srcs.utils.visualization import displayed_weeks, filter_cohorts, settings
from srcs.utils.storage import save_table
from srcs.visualization.scales import resolve


TRAJECTORY_METRICS = {
    "dispersion": "Dispersion", "n_mice": "Valid mice", "group_mean_fc": "Mean FC",
    "deviation_rms": "Mean individual LOO RMS",
}


def rebuild_phase_frames(frames, options):
    """Summarize visible weeks from saved metrics, including a partially shown phase."""
    weekly = frames["longitudinal_dispersion"]
    if weekly.empty:
        return
    configured = options.get("phases", {})
    if not configured:
        return
    old = frames.get("longitudinal_phase_dispersion")
    names = options.get("summary_phases") or (list(dict.fromkeys(old.phase)) if old is not None else list(configured))
    shown = set(weekly.week)
    phases = {name: [week for week in configured[name] if week in shown] for name in names if name in configured}
    phases = {name: weeks for name, weeks in phases.items() if weeks}
    keys = ["atlas", "cohort", "view"]
    keys += [key for key in ("side", "series") if key in weekly]
    frames["longitudinal_phase_dispersion"] = weekly_metric_means(
        weekly, phases, keys, ["dispersion", "n_mice", "group_mean_fc"])
    if "longitudinal_subjects" in frames:
        frames["longitudinal_phase_subjects"] = weekly_metric_means(
            frames["longitudinal_subjects"], phases, [*keys, "mouse"], ["deviation_rms"])


def render_trajectories(frames, colors, viz, save):
    plots = {}
    keys = ["atlas", "cohort", "side", "series", "view"]
    for prefix, axis in (("", "week"), ("phase_", "phase")):
        data = frames.get(f"longitudinal_{prefix}dispersion")
        subjects = frames.get(f"longitudinal_{prefix}subjects")
        if data is None or subjects is None or data.empty or not set(keys) <= set(data):
            continue
        animal_means = subjects.groupby([*keys, axis], sort=False).deviation_rms.mean().reset_index()
        plots[prefix] = data.merge(animal_means, on=[*keys, axis], how="left", validate="one_to_one")
    if not plots:
        return
    for (atlas, cohort), population in pd.concat(plots.values()).groupby(["atlas", "cohort"], sort=False):
        limits = {}
        for metric in TRAJECTORY_METRICS:
            finite = population.loc[np.isfinite(population[metric]), metric]
            if len(finite):
                low, high = finite.min(), finite.max()
                pad = (high - low) * .05 or max(abs(low) * .05, .01)
                limits[metric] = low - pad, high + pad
        markers = viz["view_markers"]
        for prefix, table in plots.items():
            axis = "phase" if prefix else "week"
            part = table[(table.atlas == atlas) & (table.cohort == cohort)]
            if part.empty:
                continue
            labels = list(dict.fromkeys(part[axis]))
            views = list(dict.fromkeys(part.view))
            for mode in ("mean_fc", "per_run"):
                subset = part[(part.view == "mean_fc") == (mode == "mean_fc")]
                if subset.empty:
                    continue
                fig, axes = plt.subplots(1, len(TRAJECTORY_METRICS), figsize=viz["trajectory_figsize"],
                                         squeeze=False, constrained_layout=True)
                for (side, series, view), group in subset.groupby(["side", "series", "view"], sort=False):
                    group = group.set_index(axis).reindex(labels)
                    marker = markers[views.index(view) % len(markers)]
                    for ax, metric in zip(axes[0], TRAJECTORY_METRICS):
                        ax.plot(labels, group[metric], color=colors[side], marker=marker, label=f"{series} | {view}")
                for ax, (metric, title) in zip(axes[0], TRAJECTORY_METRICS.items()):
                    ax.set_title(title)
                    if metric in limits:
                        ax.set_ylim(*limits[metric])
                    ax.tick_params(axis="x", rotation=45)
                axes[0, 0].legend(fontsize="small")
                title = "Phase means of weekly metrics" if prefix else "Weekly metrics"
                fig.suptitle(f"{atlas} | {cohort} | {mode} | {title}")
                save(fig, f"{prefix}longitudinal_{atlas}_{cohort}_{mode}")


def prepare_and_render(left, right, folder, first_label, second_label, *, viz=None, config=None, products=None):
    config = config or {}
    viz = settings(config) if viz is None else viz
    reference = config.get("descriptive", {}).get("reference_week")
    folder = Path(folder)
    destination = folder.parents[1]/"report"/"series_overlays"
    destination.mkdir(parents=True, exist_ok=True)
    sources, phase_products, frames, figures, tables = {}, [], {}, [], []
    for side, root, label in (("first", left, first_label), ("second", right, second_label)):
        for stage, filename in (("pairwise", "selection.json"), ("longitudinal", "index.json"), ("descriptive", "index.json")):
            if products is not None and f"{stage}/{filename}" not in products:
                continue
            path = root/"results"/stage/filename
            if not path.is_file():
                continue
            index = json.loads(path.read_text(encoding="utf-8"))
            sources[str(path)] = sha256(path)
            if stage == "pairwise":
                phase_products.extend((label, item) for item in index.get("phase_products", []))
                continue
            for item in index.get("files", []) + index.get("products", []):
                if sha256(item["path"]) != item["sha256"]:
                    raise ValueError(f"Changed saved product: {item['path']}")
            for name in (("dispersion", "subjects", "phase_dispersion", "phase_subjects") if stage == "longitudinal" else ("regional",)):
                csv = path.parent/(name+".csv")
                if csv.is_file():
                    data = pd.read_csv(csv)
                    if "atlas" in data and len(data):
                        frames.setdefault(stage+"_"+name, []).append(data.assign(series=label, side=side))
                        sources[str(csv)] = sha256(csv)
    # A combined comparison requires saved phase products from both series.
    if len({label for label, _ in phase_products}) == 2:
        figures.extend(render_phases(phase_products, destination/"pairwise", viz=viz, reference_week=reference))
    combined = {}
    for name, pieces in frames.items():
        table = pd.concat(pieces, ignore_index=True)
        combined[name] = table
        path = folder / (name + ".csv")
        save_table(path, table, index=False)
        tables.append(path)
    colors = dict(zip(("first", "second"), viz["series_colors"]))
    scales = []
    plot_frames = {name: filter_cohorts(frame, viz) for name, frame in combined.items()}
    for name, frame in plot_frames.items():
        if "week" in frame:
            weeks = displayed_weeks(frame.week.unique(), viz, reference)
            plot_frames[name] = frame[frame.week.isin(weeks)]
    if viz["time_start"] == "baseline" and "longitudinal_dispersion" in plot_frames:
        rebuild_phase_frames(plot_frames, config.get("pairwise", {}))
    def save(fig, name):
        for suffix in ("png", "pdf"):
            path = destination/(name+"."+suffix)
            fig.savefig(path, dpi=viz["defaults"]["dpi"], bbox_inches="tight")
            figures.append(str(path))
        plt.close(fig)
    render_trajectories(plot_frames, colors, viz, save)
    regional = plot_frames.get("descriptive_regional")
    if regional is not None and {"quantity", "level", "view", "cohort"} <= set(regional):
        for (atlas, cohort), part in regional.groupby(["atlas", "cohort"], sort=False):
            for level in ("roi", "network_pair"):
                variance = part[(part.quantity == "edge_variance") & (part.level == level)]
                if variance.empty:
                    continue
                weeks = list(dict.fromkeys(variance.week))
                regions = list(dict.fromkeys(variance.region))
                panels = list(variance.groupby(["side", "series", "view"], sort=False))
                scale = resolve(variance.mean_value, viz, "regional_variance", anchor_zero=True)
                if scale is None:
                    continue
                scales.append(scale.record(metric="regional_variance", atlas=atlas, cohort=cohort, level=level, weeks=weeks))
                fig, axes = plt.subplots(len(panels), 1, figsize=(10, max(4, 2.8*len(panels))), squeeze=False, constrained_layout=True)
                for row, ((side, series, view), group) in enumerate(panels):
                    values = group.pivot(index="region", columns="week", values="mean_value").reindex(index=regions, columns=weeks)
                    image = axes[row, 0].imshow(values, aspect="auto", **scale.image)
                    axes[row, 0].set_yticks(range(len(regions)), regions, fontsize=6)
                    axes[row, 0].set_xticks(range(len(weeks)), weeks, rotation=90, fontsize=7)
                    axes[row, 0].set_title(f"{series} | {view}")
                fig.colorbar(image, ax=list(axes.flat), label="Mean edge variance", shrink=.65, extend=scale.extend)
                save(fig, "regions_variance_"+level+"_"+atlas+"_"+cohort)
        direction_files, direction_scales = render_directions(regional, destination, viz=viz)
        figures.extend(direction_files)
    else:
        direction_scales = {}
    write_json(destination/"color_scales.json", scales + list(direction_scales.values()))
    manifest = destination / "manifest.json"
    write_json(manifest, dict(inputs=sources, direction_color_scales=direction_scales,
               figures={p: sha256(p) for p in figures},
               labels={"first": first_label, "second": second_label},
               aggregation="LOO references remain within series; group color and Run marker shared"))
    return [*tables, manifest]


def direction_values(regional):
    """Match the scale population exactly to the cells displayed in the figure."""
    data = regional[(regional.quantity == "residual_group_loo") & (regional.level == "network_pair")]
    keys = ["atlas", "cohort", "side", "series", "view", "region", "week"]
    values = data.groupby(keys, sort=False, as_index=False)[["positive_rms", "negative_rms"]].mean()
    limits = {}
    for atlas, part in values.groupby("atlas", sort=False):
        finite = part[["positive_rms", "negative_rms"]].to_numpy().ravel()
        finite = finite[np.isfinite(finite)]
        limits[atlas] = max(float(finite.max()) if len(finite) else 0., np.finfo(float).eps)
    return values, limits


def render_directions(regional, destination, dpi=None, *, viz=None):
    viz = settings() if viz is None else viz
    dpi = viz["defaults"]["dpi"] if dpi is None else dpi
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    values, _ = direction_values(regional)
    resolved = {atlas: resolve(part[["positive_rms", "negative_rms"]], viz, "regions_direction", anchor_zero=True)
                for atlas, part in values.groupby("atlas", sort=False)}
    files = []
    for (atlas, cohort), data in values.groupby(["atlas", "cohort"], sort=False):
        scale = resolved[atlas]
        if scale is None:
            continue
        weeks = list(dict.fromkeys(data.week))
        regions = list(dict.fromkeys(data.region))
        panels = list(data.groupby(["side", "series", "view"], sort=False))
        fig, axes = plt.subplots(len(panels), 2, figsize=(12, max(4, 2.8*len(panels))),
                                 squeeze=False, constrained_layout=True)
        for row, ((side, series, view), group) in enumerate(panels):
            for column, metric in enumerate(("positive_rms", "negative_rms")):
                matrix = group.pivot(index="region", columns="week", values=metric).reindex(index=regions, columns=weeks)
                image = axes[row, column].imshow(matrix, **scale.image, aspect="auto")
                axes[row, column].set_yticks(range(len(regions)), regions, fontsize=6)
                axes[row, column].set_xticks(range(len(weeks)), weeks, rotation=90, fontsize=7)
                axes[row, column].set_title(f"{series} | {view} | {metric}")
        fig.colorbar(image, ax=list(axes.flat), label="Mean individual directional RMS", shrink=.65, extend=scale.extend)
        fig.suptitle(f"{atlas} | {cohort} | shared linear scale across series, Runs and cohorts")
        for extension in ("png", "pdf"):
            path = destination/f"regions_direction_{atlas}_{cohort}.{extension}"
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            files.append(str(path))
        plt.close(fig)
    return files, {atlas: scale.record(atlas=atlas, metric="regions_direction", vmin=scale.vmin, vmax=scale.vmax,
                     population="displayed animal means", weeks=list(dict.fromkeys(values[values.atlas == atlas].week)),
                     shared_across="series, Runs, signs and cohorts", normalization="linear")
                   for atlas, scale in resolved.items() if scale is not None}
