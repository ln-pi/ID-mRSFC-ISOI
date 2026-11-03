"""Validate numerical inputs and track products at the worker boundary."""
import json
from pathlib import Path
from importlib import import_module

from srcs.core.artifacts import sha256, write_json, publish
from srcs.core.contracts import StageResult, ProductRef, StageContext
from dataclasses import asdict
from srcs.core.identity import code_closure, digest
from srcs.core.runtime import cache_key, reusable
from srcs.utils.config import PROJECT_ROOT
from srcs.utils.dataset import resolve_recording
from srcs.workflow import MODULES, PRODUCTS


PARENTS = {"edges": ("fc",), "geometry": ("edges",), "inference": ("geometry",),
           "pairwise": ("edges",), "descriptive": ("edges",),
           "stability": ("edges", "descriptive"), "longitudinal": ("edges", "descriptive"),
           "combined": ("edges",), "statistics": ("combined",)}


def record_path(ctx, stage):
    return Path(ctx.output_path) / "metadata/stages" / f"{stage}.json"


def record(ctx, stage):
    path = record_path(ctx, stage)
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


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
    return {str(path): sha256(path) for path in sorted(paths)}


def signature(ctx, stage):
    inputs = {}
    if stage in ("preprocess", "fc"):
        inputs = recording_inputs(ctx, raw=stage == "preprocess")
    elif stage == "edges":
        from srcs.utils.edge_store import source_identity
        inputs = {"_".join(a): digest(source_identity(ctx, "_".join(a))[0]) for a in ctx.atlases}
    if stage == "seedfc":
        from plugins.seedfc import input_files
        inputs = input_files(ctx,ctx.config["plugins"]["seedfc"]["options"])
    parents = PARENTS.get(stage, ())
    if stage == "fc" and record_path(ctx,"preprocess").exists():
        parents = ("preprocess",)
    if stage == "edges" and "seedfc" in ctx.data_identity["stages"]:
        parents = ("seedfc",)
    for parent in parents:
        previous = record(ctx, parent)
        if not previous and parent == "fc" and (Path(ctx.output_path) / "metadata/fc_import_provenance.json").exists():
            continue
        if not previous or not reusable(previous, signature(ctx, parent)["cache"]["key"]):
            raise ValueError(f"Inputs or products changed: rerun {parent}")
        inputs[parent] = digest(previous["cache"]["outputs"])
    option_keys = import_module(MODULES[stage]).OPTION_KEYS
    options = {key: ctx.section(key) for key in option_keys}
    code = digest(code_closure(PROJECT_ROOT, MODULES[stage]))
    return dict(schema="stage_record_v1", stage=stage, code=code, options=options,
                option_keys=option_keys, options_digest=digest(options), inputs=inputs,
                cache={"key": cache_key(code, options, inputs)})


def products(ctx, stage):
    root = Path(ctx.output_path)
    if stage in ("preprocess", "fc", "seedfc"):
        prefix = "01_preprocess" if stage == "preprocess" else "02_correlation"
        audit = max((root / "metadata/batch_reports").glob(prefix + "_*.json"),
                    key=lambda p: (p.stat().st_mtime_ns, p.name))
        files = [audit]
        if stage in ("fc","seedfc"):
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
    if reusable(record(ctx, stage), state["cache"]["key"]):
        ctx.logger.info("Reusing verified products: %s", stage)
        return StageResult(tuple(ProductRef(**item) for item in record(ctx,stage).get("products",[])))
    state["status"] = "running"
    write_json(path, state)
    try:
        refs = {}
        for parent in state["inputs"]:
            for item in record(ctx,parent).get("products",[]) if parent in MODULES else ():
                refs[item["name"]]=ProductRef(**dict(item,path=Path(item["path"])))
        ctx.stage = StageContext(refs,state["options"],state["option_keys"],
                                 Path(ctx.output_path),ctx.analysis_tag)
        try:
            action(ctx)
        finally:
            ctx.stage = None
        state["cache"]["outputs"] = products(ctx, stage)
        refs = []
        outputs = list(state["cache"]["outputs"])
        for name, patterns in PRODUCTS["fc" if stage=="seedfc" else stage].items():
            selected = [p for p in outputs if any(Path(p).match(glob) for glob in patterns)]
            if "@aggregate" in patterns:
                selected = [ctx.statistics_path]
            if not selected:
                raise ValueError(f"No files for declared product {name}")
            index = path.parent/"products"/(name+".json")
            ref = publish(index,name,name,selected)
            state["cache"]["outputs"][str(index)] = ref.sha256
            refs.append(ref)
        state["products"] = [dict(asdict(ref),path=str(ref.path)) for ref in refs]
        state["status"] = "complete"
    except BaseException as exc:
        state.update(status="failed", reason=str(exc))
        raise
    finally:
        write_json(path, state)
    return StageResult(tuple(refs))
