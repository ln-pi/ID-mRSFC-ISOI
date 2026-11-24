"""Draw saved stage similarities, with one common scale across series and Runs."""
from pathlib import Path
import json
from itertools import product
import matplotlib.pyplot as plt
import numpy as np
from srcs.core.artifacts import sha256
from srcs.utils.storage import load_arrays, save_json
from srcs.utils.visualization import displayed_weeks, settings
from srcs.visualization.scales import resolve


def render(products, folder, dpi=None, *, viz=None, reference_week=None):
    viz = settings() if viz is None else viz
    dpi = viz["defaults"]["dpi"] if dpi is None else dpi
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    loaded = []
    for series, item in products:
        if sha256(item["path"]) != item["sha256"]:
            raise ValueError("Changed stage-similarity product")
        path = Path(item["path"])
        meta, arrays = load_arrays(path, "phase_similarity")
        if viz["time_start"] == "baseline":
            selection = json.loads((path.parents[1] / "selection.json").read_text(encoding="utf-8"))
            weeks = sorted(selection["week_positions"], key=selection["week_positions"].__getitem__)
            allowed = set(displayed_weeks(weeks, viz, reference_week))
            phases = {case["phase_u"]: case["weeks_u"] for case in meta["cases"]}
            visible = {name: [week for week in phase if week in allowed] for name, phase in phases.items()}
            visible = {name: phase for name, phase in visible.items() if phase}
            if not visible:
                continue
            if phases != visible:
                from srcs.analysis.phases import similarities
                from srcs.utils.edge_store import EdgeStore
                weekly = next(entry for entry in selection["products"] if entry["atlas"] == meta["atlas"])
                if sha256(weekly["path"]) != weekly["sha256"]:
                    raise ValueError("Changed weekly pairwise product for phase redraw")
                with np.load(weekly["path"], allow_pickle=False) as archive:
                    edge = json.loads(str(archive["metadata"]))["edge_store"]
                if sha256(edge["path"]) != edge["sha256"]:
                    raise ValueError("Changed edge store for phase redraw")
                store = EdgeStore.load(edge["path"])
                cases, values, counts = similarities(store, dict(phases=visible, summary_phases=list(visible)))
                meta = {**meta, "cases": cases}
                arrays = dict(S_r=values, n_pairs=counts)
        loaded.append((series, meta, arrays))
    files, scales = [], []
    atlases = list(dict.fromkeys(meta["atlas"] for _, meta, _ in loaded))
    for atlas in atlases:
        selected = [(series, meta, arrays) for series, meta, arrays in loaded if meta["atlas"] == atlas]
        cases = selected[0][1]["cases"]
        if any(meta["cases"] != cases for _, meta, _ in selected):
            raise ValueError("Combined stage panels require identical phase weeks and Run definitions")
        scale = resolve(np.concatenate([arrays["S_r"].ravel() for _, _, arrays in selected]),
                        viz, "pairwise", bounds=(-1, 1))
        if scale is None:
            continue
        phase_weeks = {case["phase_u"]: case["weeks_u"] for case in cases}
        scales.append(scale.record(metric="phase_pairwise", atlas=atlas, phases=phase_weeks))
        settings = list(dict.fromkeys((c["cohort"], c["method"], c["query_run"], c["reference_run"]) for c in cases))
        phases = list(dict.fromkeys(c["phase_u"] for c in cases))
        for cohort, method, query, reference in settings:
            fig, axes = plt.subplots(len(selected)*len(phases), len(phases), squeeze=False,
                figsize=(3.4*len(phases), 3*len(phases)*len(selected)), constrained_layout=True)
            for group, (series, meta, arrays) in enumerate(selected):
                lookup = {(c["cohort"], c["method"], c["query_run"], c["reference_run"], c["phase_u"], c["phase_v"]): i
                          for i, c in enumerate(meta["cases"])}
                mice = meta["mouse_ids"]
                for u, v in product(range(len(phases)), repeat=2):
                    ax = axes[group*len(phases)+u, v]
                    key = (cohort, method, query, reference, phases[u], phases[v])
                    if key not in lookup:
                        ax.set_axis_off()
                        continue
                    matrix = arrays["S_r"][lookup[key]]
                    image = ax.imshow(matrix, **scale.image)
                    ax.set_xticks(range(len(mice)), mice, rotation=90, fontsize=5)
                    ax.set_yticks(range(len(mice)), mice, fontsize=5)
                    ax.set_title(f"{series} | {phases[u]} to {phases[v]}", fontsize=8)
                    ax.set(xlabel=reference+" reference", ylabel=query+" query")
            fig.colorbar(image, ax=list(axes.flat), label="Pattern Pearson r", shrink=.65, extend=scale.extend)
            fig.suptitle(f"{atlas} | {cohort} | {method} | {query} to {reference}")
            for extension in ("png", "pdf"):
                path = folder/f"phase_{atlas}_{cohort}_{method}_{query}_{reference}.{extension}"
                fig.savefig(path, dpi=dpi, bbox_inches="tight")
                files.append(str(path))
            plt.close(fig)
    save_json(folder / "phase_color_scales.json", scales)
    return files
