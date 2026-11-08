"""Draw saved stage similarities, with one common scale across series and Runs."""
from pathlib import Path
from itertools import product
import matplotlib.pyplot as plt
import numpy as np
from srcs.core.artifacts import sha256
from srcs.utils.storage import load_arrays


def render(products, folder, dpi=150):
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=True)
    loaded = []
    for series, item in products:
        if sha256(item["path"]) != item["sha256"]:
            raise ValueError("Changed stage-similarity product")
        meta, arrays = load_arrays(item["path"], "phase_similarity")
        loaded.append((series, meta, arrays))
    files = []
    atlases = list(dict.fromkeys(meta["atlas"] for _, meta, _ in loaded))
    for atlas in atlases:
        selected = [(series, meta, arrays) for series, meta, arrays in loaded if meta["atlas"] == atlas]
        cases = selected[0][1]["cases"]
        signatures = [meta["cases"] for _, meta, _ in selected]
        if any(value != signatures[0] for value in signatures[1:]):
            raise ValueError("Combined stage panels require identical phase weeks and Run definitions")
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
                    image = ax.imshow(matrix, vmin=-1, vmax=1,
                        cmap=plt.get_cmap("inferno").with_extremes(bad="lightgray"))
                    ax.set_xticks(range(len(mice)), mice, rotation=90, fontsize=5)
                    ax.set_yticks(range(len(mice)), mice, fontsize=5)
                    ax.set_title(f"{series} | {phases[u]} to {phases[v]}", fontsize=8)
                    ax.set(xlabel=reference+" reference", ylabel=query+" query")
            fig.colorbar(image, ax=list(axes.flat), label="Pattern Pearson r", shrink=.65)
            fig.suptitle(f"{atlas} | {cohort} | {method} | {query} to {reference}")
            for extension in ("png", "pdf"):
                path = folder/f"phase_{atlas}_{cohort}_{method}_{query}_{reference}.{extension}"
                fig.savefig(path, dpi=dpi, bbox_inches="tight")
                files.append(str(path))
            plt.close(fig)
    return files
