"""Compare saved FC reports using shared axes and compatible acquisition metadata."""
import json
from pathlib import Path

import numpy as np

from srcs.core.artifacts import artifact_path, publish, sha256, write_json
from srcs.core.identity import digest
from srcs.core.runtime import exclusive_run
from srcs.utils.config import PROJECT_ROOT, load_config
from srcs.utils.edge_store import EdgeStore
from srcs.utils.logger import close_log_file, setup_logger
from srcs.utils.visualization import merge, saved_config, settings, validate
from srcs.workflow.batch_layout import BatchLayout, dataset_label, saved_scopes, series_input
from srcs.workflow.plugins import enabled_plugins, enable_functions
from srcs.workflow.products import verified_index

REPORT_PRODUCTS = ("pairwise/selection.json", "descriptive/index.json", "longitudinal/index.json")
SHARED_DISPLAY = ("defaults", "time_start", "colormaps", "color_scale", "metrics", "missing_color",
                  "under_color", "over_color", "display_cohorts", "series_colors",
                  "view_markers", "trajectory_figsize")


def load_series(root):
    identity = root / "metadata/data_identity.json"
    scope = json.loads(identity.read_text(encoding="utf-8"))["scope"]
    hashes = {identity.relative_to(root).as_posix(): sha256(identity)}
    products = set()
    positions = None
    for relative in REPORT_PRODUCTS:
        path = root / "results" / relative
        if not path.is_file():
            continue
        index = verified_index(path)
        if relative.startswith("pairwise/"):
            positions = index["week_positions"]
        products.add(relative)
        hashes[path.relative_to(root).as_posix()] = sha256(path)
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


def _comparison_scope_mode(first, second):
    """Classify the requested pair before checking saved products and acquisition axes."""
    paired = first["series"] == second["series"]
    if paired == (first["branch"] == second["branch"]):
        raise ValueError("Compare two series at one branch or two branches within one series")
    return "paired" if paired else "independent"


def comparison_mode(left, right):
    a, _, _, positions_a, stores_a = left
    b, _, _, positions_b, stores_b = right
    if any(a[key] != b[key] for key in ("partitions", "weeks", "runs")):
        raise ValueError("Comparison requires identical partitions, weeks and Runs")
    if positions_a is not None and positions_b is not None and positions_a != positions_b:
        raise ValueError("Comparison requires identical numerical week positions")
    mode = _comparison_scope_mode(a, b)
    for atlas, x in stores_a.items():
        y = stores_b[atlas]
        if x.roi_labels != y.roi_labels or not np.array_equal(x.edge_index, y.edge_index):
            raise ValueError("Comparison requires identical ordered ROIs and edges")
        if mode == "paired" and (set(x.mice) != set(y.mice) or not np.array_equal(
                x.valid_record, y.valid_record[[y.mice.index(mouse) for mouse in x.mice]])):
            raise ValueError("Paired branches differ in animals or valid acquisition membership")
    return mode


def compare(first, second, *, destination, configs=None, function_results=None):
    """Write core and optional plugin comparisons into one batch comparison root."""
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
    renderer = [sha256(PROJECT_ROOT / "srcs/visualization" / module)
                for module in ("series_overlays.py", "phase_pairwise.py", "scales.py")]
    identity = dict(schema="series_comparison", mode=mode, inputs=[item[1] for item in inputs],
                    visualization=viz, pairwise=display_config.get("pairwise", {}), reference=reference,
                    renderer=renderer, comparison_code=sha256(__file__))
    parent = Path(destination).resolve()
    folder = parent / "results/core"
    report = parent / "report/core"
    batch_root = parent.parent
    with exclusive_run(parent / "metadata/core.lock"):
        from srcs.visualization.series_overlays import prepare_and_render
        folder.mkdir(parents=True, exist_ok=True)
        paths = prepare_and_render(first, second, folder, *labels, viz=viz, config=display_config,
                                   products=common, report_destination=report, batch_root=batch_root)
        marker = folder / "identity.json"
        refs = [artifact_path(root, folder) if root.is_relative_to(batch_root) else str(root)
                for root in (first, second)]
        write_json(marker, {**identity, "comparison_key": digest(identity), "input_roots": refs})
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
        results[name] = str(plugin.entry(spec["run"])(first, second, configs=configs, destination=parent))
    return folder


def compare_saved(first, second, *, project_root, enabled_functions=(), visualization_config=None,
                  output_root=None, dry_run=False):
    """Create a standalone comparison batch from two saved series."""
    from srcs.workflow.batch import timestamp

    project = Path(project_root).resolve()
    roots = [series_input(project, name) for name in (first, second)]
    configs = [enable_functions(saved_config(root), enabled_functions) for root in roots]
    scopes = list(saved_scopes(roots, configs))
    _comparison_scope_mode(*scopes)
    if visualization_config:
        override = load_config(project / visualization_config)
        configs = [merge(config, {"visualization": override.get("visualization", override)}) for config in configs]
    parent = project / (output_root or configs[0]["paths"]["output"])
    allocate = BatchLayout.preview if dry_run else BatchLayout.create
    layout = allocate(parent, dataset_label(scopes), "Compare")
    print(f"Batch output: {layout.root}", flush=True)
    if dry_run:
        for config, root in zip(configs, roots):
            validate(settings(config))
            print(f"Saved series: {root}")
        return 0

    summary = dict(kind="compare", status="running", inputs=[str(root) for root in roots],
                   comparison=dict(status="running", path="compare"), started=timestamp())
    layout.save_summary(summary)
    logger = setup_logger(str(layout.logger_path("batch")))
    try:
        function_results = {}
        result = compare(*roots, configs=configs, function_results=function_results, destination=layout.compare)
        summary["comparison"].update(status="complete", result=layout.relative(result),
                                     functions={name: layout.relative(path) for name, path in function_results.items()})
        summary["status"] = "complete"
        logger.info("Comparison saved under %s", layout.compare)
        return 0
    except BaseException as exc:
        status = "interrupted" if isinstance(exc, KeyboardInterrupt) else "failed"
        summary.update(status=status, reason=str(exc))
        summary["comparison"].update(status=status, reason=str(exc))
        logger.error("Comparison %s: %s", status, exc)
        raise
    finally:
        summary["finished"] = timestamp()
        layout.save_summary(summary)
        close_log_file()
