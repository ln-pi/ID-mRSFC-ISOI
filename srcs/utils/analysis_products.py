"""Verify saved dependencies, including their link to the current edge store."""
import json
from pathlib import Path

from srcs.analysis.connectome import load_arrays
from srcs.analysis.magnitude import SCHEMA
from srcs.utils.provenance import sha256
from srcs.utils.results import save_json, save_table


def verified_index(path):
    path = Path(path)
    index = json.loads(path.read_text(encoding="utf-8"))
    for item in index.get("products", []) + index.get("files", []):
        if sha256(item["path"]) != item["sha256"]:
            raise ValueError(f"Saved product changed: {item['path']}; rerun its numerical stage")
    return index


def descriptive_products(ctx, atlas):
    path = Path(ctx.output_path)/"results/descriptive/v2/index.json"
    index = verified_index(path)
    if index["schema"] != SCHEMA:
        raise ValueError("Run descriptive v2 first")
    products = {}
    for item in index["products"]:
        if item["atlas"] != atlas:
            continue
        meta, arrays = load_arrays(item["path"], SCHEMA)
        edge = Path(ctx.output_path)/"results/edges"/f"{atlas}.npz"
        if sha256(edge) != meta["edge_store"]["sha256"]:
            raise ValueError("Descriptive input edge store changed; rerun descriptive")
        products[item["cohort"]] = (meta, arrays)
    if set(products) != {"available_week", "fixed_window"}:
        raise ValueError(f"Incomplete descriptive products for {atlas}")
    return products, dict(path=str(path), sha256=sha256(path))


def bootstrap_options(ctx, section):
    cfg = ctx.section(section)
    count, seed = cfg.get("n_bootstrap", 200), cfg.get("seed", 42)
    if isinstance(count, bool) or not isinstance(count, int) or count < 0 or (0 < count < 20):
        raise ValueError(f"{section}.n_bootstrap must be 0 or an integer >=20")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError(f"{section}.seed must be a nonnegative integer")
    return count, seed


def publish(folder, schema, products, tables, methods):
    import pandas as pd
    for name, rows in tables.items():
        frame = rows if isinstance(rows, pd.DataFrame) else pd.DataFrame(rows)
        if not len(frame.columns):
            frame = pd.DataFrame(columns=["atlas", "status"])
        save_table(folder/f"{name}.csv", frame, index=False)
    save_json(folder/"methods.json", methods)
    files = [folder/f"{name}.csv" for name in tables]+[folder/"methods.json"]
    save_json(folder/"index.json", dict(schema=schema, products=products,
        files=[dict(path=str(p), sha256=sha256(p)) for p in files]))
