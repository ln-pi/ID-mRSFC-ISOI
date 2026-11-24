"""Compare saved FC reports using shared axes and compatible acquisition metadata."""
import json
from pathlib import Path

import numpy as np

from srcs.core.artifacts import publish, sha256, write_json
from srcs.core.identity import digest
from srcs.core.runtime import exclusive_run
from srcs.utils.config import PROJECT_ROOT
from srcs.utils.edge_store import EdgeStore
from srcs.utils.visualization import saved_config, settings, validate
from srcs.workflow.plugins import enabled_plugins
from srcs.workflow.products import verified_index

REPORT_PRODUCTS = ("pairwise/selection.json", "descriptive/index.json", "longitudinal/index.json")
SHARED_DISPLAY = ("defaults", "time_start", "colormaps", "color_scale", "metrics", "missing_color",
                  "under_color", "over_color", "display_cohorts", "series_colors",
                  "view_markers", "trajectory_figsize")


def load_series(root):
    identity = root / "metadata/data_identity.json"
    scope = json.loads(identity.read_text(encoding="utf-8"))["scope"]
    hashes = {str(identity): sha256(identity)}
    products = set()
    positions = None
    for relative in REPORT_PRODUCTS:
        path = root / "results" / relative
        if not path.is_file():
            continue
        index = verified_index(path)
        for product in index.get("phase_products", ()):
            if sha256(product["path"]) != product["sha256"]:
                raise ValueError(f"Changed phase product: {product['path']}")
        if relative.startswith("pairwise/"):
            positions = index["week_positions"]
        products.add(relative)
        hashes[str(path)] = sha256(path)
    stores = {}
    for partition in scope["partitions"]:
        atlas = partition["name"]
        path = root / "results/edges" / f"{atlas}.npz"
        sidecar = path.with_suffix(".json")
        checksum = sha256(path)
        if checksum != json.loads(sidecar.read_text(encoding="utf-8"))["sha256"]:
            raise ValueError(f"Changed edge store: {path}")
        store = EdgeStore.load(path)
        if (store.mice, store.weeks, store.runs) != (scope["mice"], scope["weeks"], scope["runs"]):
            raise ValueError(f"Edge store acquisition axes differ from saved identity: {path}")
        stores[atlas] = store
        hashes[str(path)] = checksum
    return scope, hashes, products, positions, stores


def comparison_mode(left, right):
    a, _, _, positions_a, stores_a = left
    b, _, _, positions_b, stores_b = right
    if any(a[key] != b[key] for key in ("partitions", "weeks", "runs")):
        raise ValueError("Comparison requires identical partitions, weeks and Runs")
    if positions_a is not None and positions_b is not None and positions_a != positions_b:
        raise ValueError("Comparison requires identical numerical week positions")
    paired = a["series"] == b["series"]
    if paired == (a["branch"] == b["branch"]):
        raise ValueError("Compare two series at one branch or two branches within one series")
    for atlas, x in stores_a.items():
        y = stores_b[atlas]
        if x.roi_labels != y.roi_labels or not np.array_equal(x.edge_index, y.edge_index):
            raise ValueError("Comparison requires identical ordered ROIs and edges")
        if paired and (set(x.mice) != set(y.mice) or not np.array_equal(
                x.valid_record, y.valid_record[[y.mice.index(mouse) for mouse in x.mice]])):
            raise ValueError("Paired branches differ in animals or valid acquisition membership")
    return "paired" if paired else "independent"


def compare(first, second, output_root=None, *, configs=None, function_results=None):
    first, second = (Path(path).resolve() for path in (first, second))
    configs = configs or [saved_config(root) for root in (first, second)]
    display_config = configs[0]
    viz, other_viz = [settings(config) for config in configs]
    validate(viz)
    validate(other_viz)
    if any(viz.get(key) != other_viz.get(key) for key in SHARED_DISPLAY):
        raise ValueError("Comparison requires matching shared visualization settings")
    reference, other_reference = [c.get("descriptive", {}).get("reference_week") for c in configs]
    phases, other_phases = [c.get("pairwise", {}).get("phases", {}) for c in configs]
    if viz["time_start"] == "baseline" and (reference != other_reference or phases != other_phases):
        raise ValueError("Comparison requires matching plot baselines and phase definitions")
    inputs = [load_series(root) for root in (first, second)]
    mode = comparison_mode(*inputs)
    common = inputs[0][2] & inputs[1][2]
    if not common:
        raise ValueError("No shared pairwise, descriptive or longitudinal products to compare")
    a, b = [item[0] for item in inputs]
    labels = [item["series"] if mode == "independent" else f"{item['series']}/{item['branch']}" for item in (a,b)]
    base = Path(output_root).resolve() if output_root else first.parent
    renderer = [sha256(PROJECT_ROOT / "srcs/visualization" / module)
                for module in ("series_overlays.py", "phase_pairwise.py", "scales.py")]
    identity = dict(schema="series_comparison", mode=mode, inputs=[item[1] for item in inputs],
                    visualization=viz, pairwise=display_config.get("pairwise", {}), reference=reference,
                    renderer=renderer, comparison_code=sha256(__file__))
    parent = base / "_compare" / f"{first.name}__vs__{second.name}_{digest(identity)[:8]}"
    folder = parent / "results/compare"
    with exclusive_run(parent / "compare.lock"):
        from srcs.visualization.series_overlays import prepare_and_render
        folder.mkdir(parents=True, exist_ok=True)
        paths = prepare_and_render(first, second, folder, *labels, viz=viz, config=display_config,
                                   products=common)
        marker = folder / "identity.json"
        write_json(marker, identity)
        publish(folder/"index.json", "comparison", "series_comparison", [marker, *paths])
    results = function_results if function_results is not None else {}
    functions = [{p.spec["id"]: p for p in enabled_plugins(config, "function", PROJECT_ROOT)} for config in configs]
    for name, plugin in functions[0].items():
        spec = plugin.spec.get("comparison")
        if not spec or name not in functions[1]:
            continue
        missing = [root / relative for root in (first, second) for relative in spec["requires"]
                   if not (root / relative).is_file()]
        if missing:
            raise FileNotFoundError(f"{name} comparison requires saved products: {missing}")
        results[name] = str(plugin.entry(spec["run"])(first, second, output_root=base, configs=configs))
    return folder
