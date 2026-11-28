"""Validate numerical inputs and track products at the worker boundary."""
import json
from time import monotonic
from pathlib import Path
from importlib import import_module
from tqdm import tqdm

from srcs.core.artifacts import artifact_path, resolve_artifact, sha256, write_json, publish
from srcs.core.contracts import StageResult, ProductRef, StageContext
from dataclasses import asdict
from srcs.core.identity import code_closure, digest
from srcs.core.runtime import cache_key, reusable
from srcs.utils.config import PROJECT_ROOT
from srcs.utils.dataset import resolve_recording
from srcs.workflow.plugins import function_options


def record_path(ctx, stage):
    return Path(ctx.output_path) / "metadata/stages" / f"{stage}.json"


def record(ctx, stage):
    path = record_path(ctx, stage)
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def _product_refs(state, anchor):
    return tuple(ProductRef(**dict(item, path=resolve_artifact(item["path"], anchor)))
                 for item in state.get("products", []))


def _output_ref(ctx, path, anchor):
    # Processed signals can be written outside the series by preprocessing.
    path = Path(path).resolve()
    return artifact_path(path, anchor) if path.is_relative_to(Path(ctx.output_path).resolve()) else str(path)


def recording_inputs(ctx, raw=False):
    paths = set(Path(ctx.mask_root).glob("*.png"))
    if raw:
        paths = {path for path in paths if path.name.casefold()=="chamber_mask.png"}
        paths.update(PROJECT_ROOT/"srcs/atlas"/name for name in
                     ("Bregma.pkl","Raw_Area_ID_bin4.npy","Brain_Area_Masks_bin4.pkl"))
    for week in ctx.weeks:
        for mouse in ctx.mice:
            for run in ctx.runs:
                try:
                    source = resolve_recording(ctx, week, mouse, run, raw=raw)
                except FileNotFoundError:
                    continue
                if raw:
                    paths.update(source.glob("img*.bin"))
                    paths.update(source.glob("info.txt"))
                else:
                    paths.add(source)
    paths = sorted(paths)
    total_bytes = sum(path.stat().st_size for path in paths)
    ctx.logger.info("Verifying input contents: %d files, %.2f GiB; FC computation starts after verification",
                    len(paths), total_bytes / 1024**3)
    started = monotonic()
    inputs = {}
    with tqdm(total=total_bytes, desc="Verifying inputs", unit="B", unit_scale=True) as progress:
        for path in paths:
            progress.set_postfix_str(path.name, refresh=True)
            inputs[str(path)] = sha256(path)
            progress.update(path.stat().st_size)
    ctx.logger.info("Input verification completed in %.1f seconds", monotonic() - started)
    return inputs


def signature(ctx, stage):
    inputs = {}
    if stage in ("preprocess", "fc"):
        inputs = recording_inputs(ctx, raw=stage == "preprocess")
    elif stage == "edges":
        from srcs.utils.edge_store import source_identity
        inputs = {"_".join(a): digest(source_identity(ctx, "_".join(a))[0]) for a in ctx.atlases}
    definition = ctx.workflow.stages[stage]
    if definition.source:
        source = ctx.workflow.source
        inputs = source.entry(":input_files")(ctx, source.options)
    parents = definition.parents
    if stage == "fc" and record_path(ctx,"preprocess").exists():
        parents = ("preprocess",)
    for parent in parents:
        previous = record(ctx, parent)
        if not previous and parent == "fc" and (Path(ctx.output_path) / "metadata/fc_import_provenance.json").exists():
            continue
        if not previous or not reusable(previous, signature(ctx, parent)["cache"]["key"],
                                        record_path(ctx, parent).parent):
            raise ValueError(f"Inputs or products changed: rerun {parent}")
        inputs[parent] = digest(previous["cache"]["outputs"])
    option_keys = import_module(definition.module).OPTION_KEYS
    options = {}
    for key in option_keys:
        owner = ctx.workflow.option_sections.get(key)
        options[key] = function_options(ctx.config, owner, key) if owner else ctx.section(key)
    if definition.source:
        options["plugins"] = {"source": {source.spec["id"]: {"options": source.options}}}
    code = digest(code_closure(PROJECT_ROOT, definition.module))
    return dict(schema="stage_record", stage=stage, code=code, options=options,
                option_keys=option_keys, options_digest=digest(options), inputs=inputs,
                cache={"key": cache_key(code, options, inputs)})


def products(ctx, stage):
    root = Path(ctx.output_path)
    if stage in ("preprocess", "fc") or ctx.workflow.stages[stage].source:
        prefix = "01_preprocess" if stage == "preprocess" else "02_correlation"
        audit = max((root / "metadata/batch_reports").glob(prefix + "_*.json"),
                    key=lambda p: (p.stat().st_mtime_ns, p.name))
        files = [audit]
        if stage != "preprocess":
            files.extend(Path(ctx.correlation_path).glob("*_excluded.pkl"))
        else:
            files.extend(Path(p) for p in recording_inputs(ctx) if p.endswith("_FiltDRR.npy"))
    else:
        files = [p for p in (root / "results" / stage).rglob("*") if p.is_file()]
    if not files:
        raise ValueError(f"{stage} produced no files")
    return {str(p.resolve()): sha256(p) for p in files}


def run_cached(ctx, stage, action):
    state = signature(ctx, stage)
    path = record_path(ctx, stage)
    previous = record(ctx, stage)
    if reusable(previous, state["cache"]["key"], path.parent):
        ctx.logger.info("Reusing verified products: %s", stage)
        return StageResult(_product_refs(previous, path.parent))
    state["status"] = "running"
    write_json(path, state)
    try:
        refs = {}
        for parent in state["inputs"]:
            if parent in ctx.workflow.stages:
                for ref in _product_refs(record(ctx, parent), record_path(ctx, parent).parent):
                    refs[ref.name] = ref
        ctx.stage = StageContext(refs,state["options"],state["option_keys"],
                                 Path(ctx.output_path),ctx.analysis_tag)
        try:
            action(ctx)
        finally:
            ctx.stage = None
        results = products(ctx, stage)
        state["cache"]["outputs"] = {_output_ref(ctx, file, path.parent): checksum
                                     for file, checksum in results.items()}
        refs = []
        outputs = list(results)
        for name, patterns in ctx.workflow.stages[stage].products.items():
            selected = [p for p in outputs if any(Path(p).match(glob) for glob in patterns)]
            if "@aggregate" in patterns:
                selected = [ctx.statistics_path]
            if not selected:
                raise ValueError(f"No files for declared product {name}")
            index = path.parent/"products"/(name+".json")
            ref = publish(index, name, name, selected, artifact_root=ctx.output_path)
            state["cache"]["outputs"][_output_ref(ctx, index, path.parent)] = ref.sha256
            refs.append(ref)
        state["products"] = [dict(asdict(ref), path=_output_ref(ctx, ref.path, path.parent)) for ref in refs]
        state["status"] = "complete"
    except BaseException as exc:
        state.update(status="failed", reason=str(exc))
        raise
    finally:
        write_json(path, state)
    return StageResult(tuple(refs))
