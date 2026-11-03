"""Summarize directional matching at animal level with leave-one-animal sensitivity."""
from pathlib import Path
import numpy as np
import pandas as pd

from srcs.core.artifacts import publish
from srcs.utils.results import save_table


def run(ctx):
    folder = Path(ctx.output_dir("results","inference","v1"))
    source = Path(ctx.output_path)/"results/geometry/v1/subjects.csv"
    subjects = pd.read_csv(source)
    endpoint = ctx.section("inference").get("primary_endpoint","q_mean")
    if endpoint not in ("q_mean","margin_mean"):
        raise ValueError("Unknown primary endpoint")
    metric = endpoint.removesuffix("_mean")
    rows = []
    keys = ["atlas","week","cohort","run_a","run_b","metric"]
    for labels,frame in subjects.groupby(keys,sort=False):
        def animal_mean(group):
            return group.groupby("mouse")[metric].mean().mean()
        full = frame[frame.rule == "full"]
        estimates = [animal_mean(group) for _,group in frame[frame.rule=="leave_one_out"].groupby("removed")]
        finite = np.asarray([v for v in estimates if np.isfinite(v)])
        low, high = np.quantile(finite,[0.05,0.95]) if len(finite) else (np.nan,np.nan)
        rows.append(dict(zip(keys,labels), endpoint=endpoint, estimate=animal_mean(full),
                         low=low, high=high, interval="not_confidence_interval", n=full.mouse.nunique()))
    path = folder/"summary.csv"
    save_table(path,pd.DataFrame(rows),index=False)
    publish(folder/"index.json","inference_summary_v1","inference_summary_v1",[path])
