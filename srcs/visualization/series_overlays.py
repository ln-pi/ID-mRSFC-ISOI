"""Shared-coordinate comparisons of saved series products."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from srcs.core.artifacts import sha256, write_json
from srcs.visualization.phase_pairwise import render as render_phases


def prepare_and_render(left, right, folder, first_label, second_label):
    folder = Path(folder)
    destination = folder.parents[1]/"report"/"series_overlays"
    destination.mkdir(parents=True, exist_ok=True)
    sources, phase_products, frames, figures = {}, [], {}, []
    for side, root, label in (("first", left, first_label), ("second", right, second_label)):
        for stage, filename in (("pairwise", "selection.json"), ("longitudinal", "index.json"), ("descriptive", "index.json")):
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
        figures.extend(render_phases(phase_products, destination/"pairwise"))
    combined = {}
    for name, pieces in frames.items():
        table = pd.concat(pieces, ignore_index=True)
        combined[name] = table
        table.to_csv(folder/(name+".csv"), index=False)
    colors = {"first": "#356caa", "second": "#bd502f"}
    def save(fig, name):
        for suffix in ("png", "pdf"):
            path = destination/(name+"."+suffix)
            fig.savefig(path, dpi=150, bbox_inches="tight")
            figures.append(str(path))
        plt.close(fig)
    for prefix in ("", "phase_"):
        data = combined.get("longitudinal_"+prefix+"dispersion")
        subjects = combined.get("longitudinal_"+prefix+"subjects")
        if data is None or subjects is None or not {"view", "cohort"} <= set(data):
            continue
        axis = "phase" if prefix else "week"
        labels = list(dict.fromkeys(data[axis]))
        for (atlas, cohort), part in data.groupby(["atlas", "cohort"], sort=False):
            fig, axes = plt.subplots(1, 4, figsize=(19, 4.5), constrained_layout=True)
            sub = subjects[(subjects.atlas == atlas) & (subjects.cohort == cohort)]
            markers = {view: marker for view, marker in zip(data.view.unique(), ("o", "^", "s", "D", "x"))}
            for (side, series, view), group in part.groupby(["side", "series", "view"], sort=False):
                group = group.set_index(axis).reindex(labels)
                style = dict(color=colors[side], marker=markers[view], label=f"{series} | {view}")
                for ax, metric in zip(axes, ("dispersion", "n_mice", "group_mean_fc")):
                    ax.plot(labels, group[metric], **style)
                animal = sub[(sub.side == side) & (sub.view == view)]
                # Animal means follow per-animal phase means, preserving equal animal weights.
                mean = animal.groupby(axis).deviation_rms.mean().reindex(labels)
                axes[3].plot(labels, mean, **style)
            for ax, title in zip(axes, ("Dispersion", "Valid mice", "Mean FC", "Mean individual LOO RMS")):
                ax.set_title(title)
                ax.tick_params(axis="x", rotation=45)
            axes[0].legend(fontsize=6)
            fig.suptitle(f"{atlas} | {cohort} | {'Mean weekly metrics by phase' if prefix else 'Weekly metrics'}")
            save(fig, prefix+"longitudinal_"+atlas+"_"+cohort)
    regional = combined.get("descriptive_regional")
    if regional is not None and {"quantity", "level", "view", "cohort"} <= set(regional):
        for (atlas, cohort), part in regional.groupby(["atlas", "cohort"], sort=False):
            for level in ("roi", "network_pair"):
                variance = part[(part.quantity == "edge_variance") & (part.level == level)]
                if variance.empty:
                    continue
                weeks = list(dict.fromkeys(variance.week))
                regions = list(dict.fromkeys(variance.region))
                panels = list(variance.groupby(["side", "series", "view"], sort=False))
                limit = max(float(variance.mean_value.max()), np.finfo(float).eps)
                fig, axes = plt.subplots(len(panels), 1, figsize=(10, max(4, 2.8*len(panels))), squeeze=False, constrained_layout=True)
                for row, ((side, series, view), group) in enumerate(panels):
                    values = group.pivot(index="region", columns="week", values="mean_value").reindex(index=regions, columns=weeks)
                    image = axes[row, 0].imshow(values, vmin=0, vmax=limit, aspect="auto",
                        cmap=plt.get_cmap("inferno").with_extremes(bad="lightgray"))
                    axes[row, 0].set_yticks(range(len(regions)), regions, fontsize=6)
                    axes[row, 0].set_xticks(range(len(weeks)), weeks, rotation=90, fontsize=7)
                    axes[row, 0].set_title(f"{series} | {view}")
                fig.colorbar(image, ax=list(axes.flat), label="Mean edge variance", shrink=.65)
                save(fig, "regions_variance_"+level+"_"+atlas+"_"+cohort)
        direction_files, direction_scales = render_directions(regional, destination)
        figures.extend(direction_files)
    else:
        direction_scales = {}
    write_json(destination/"manifest.json", dict(inputs=sources, direction_color_scales=direction_scales,
               figures={p: sha256(p) for p in figures},
               labels={"first": first_label, "second": second_label},
               aggregation="LOO references remain within series; group color and Run marker shared"))


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


def render_directions(regional, destination, dpi=150):
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=True)
    values, limits = direction_values(regional)
    files = []
    for (atlas, cohort), data in values.groupby(["atlas", "cohort"], sort=False):
        weeks = list(dict.fromkeys(data.week))
        regions = list(dict.fromkeys(data.region))
        panels = list(data.groupby(["side", "series", "view"], sort=False))
        fig, axes = plt.subplots(len(panels), 2, figsize=(12, max(4, 2.8*len(panels))),
                                 squeeze=False, constrained_layout=True)
        for row, ((side, series, view), group) in enumerate(panels):
            for column, metric in enumerate(("positive_rms", "negative_rms")):
                matrix = group.pivot(index="region", columns="week", values=metric).reindex(index=regions, columns=weeks)
                image = axes[row, column].imshow(matrix, vmin=0, vmax=limits[atlas],
                    cmap=plt.get_cmap("inferno").with_extremes(bad="lightgray"), aspect="auto")
                axes[row, column].set_yticks(range(len(regions)), regions, fontsize=6)
                axes[row, column].set_xticks(range(len(weeks)), weeks, rotation=90, fontsize=7)
                axes[row, column].set_title(f"{series} | {view} | {metric}")
        fig.colorbar(image, ax=list(axes.flat), label="Mean individual directional RMS", shrink=.65)
        fig.suptitle(f"{atlas} | {cohort} | shared linear scale across series, Runs and cohorts")
        for extension in ("png", "pdf"):
            path = destination/f"regions_direction_{atlas}_{cohort}.{extension}"
            fig.savefig(path, dpi=dpi, bbox_inches="tight")
            files.append(str(path))
        plt.close(fig)
    return files, {atlas: dict(vmin=0, vmax=limit, population="displayed animal means",
                              shared_across="series, Runs, signs and cohorts", normalization="linear")
                   for atlas, limit in limits.items()}
