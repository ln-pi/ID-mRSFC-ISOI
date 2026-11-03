"""Summarize directional matching at animal level with leave-one-animal sensitivity."""
from pathlib import Path
import numpy as np
import pandas as pd

from srcs.core.artifacts import publish
from srcs.utils.storage import save_table
from srcs.analysis.longitudinal import linear_slope
from srcs.analysis.pairwise import week_positions


from srcs.core.runtime import numerical_stage

OPTION_KEYS = ("inference", "pairwise",)


@numerical_stage("inference", OPTION_KEYS)
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
    full = subjects[subjects.rule=="full"]
    positions=week_positions(ctx.weeks, ctx.section("pairwise").get("week_positions"))
    options=ctx.section("inference")
    rng=np.random.default_rng(options.get("seed",42))
    trends=[]
    for labels,frame in full.groupby(["atlas","cohort","run_a","run_b","metric"],sort=False):
        by_mouse=frame.pivot_table(index="mouse",columns="week",values=metric).reindex(columns=ctx.weeks)
        if by_mouse.empty:
            continue
        values=by_mouse.to_numpy()
        time=np.array([positions[w] for w in ctx.weeks])
        def estimate(sample):
            averages=pd.DataFrame(values[sample]).mean(axis=0).to_numpy()
            return linear_slope(time,averages)
        n_bootstrap=options.get("n_bootstrap",200)
        if not isinstance(n_bootstrap,int) or n_bootstrap<0:
            raise ValueError("inference.n_bootstrap must be a nonnegative integer")
        slopes=[estimate(rng.integers(len(values),size=len(values))) for _ in range(n_bootstrap)]
        finite=np.asarray([v for v in slopes if np.isfinite(v)])
        low,high=np.quantile(finite,[.05,.95]) if len(finite) else (np.nan,np.nan)
        trends.append(dict(zip(["atlas","cohort","run_a","run_b","metric"],labels),
                           endpoint=endpoint,slope=estimate(np.arange(len(values))),
                           low=low,high=high,n=len(values),interval="not_confidence_interval"))
    trend_path=folder/"trends.csv"
    save_table(trend_path,pd.DataFrame(trends),index=False)
    publish(folder/"index.json","inference_summary_v1","inference_summary_v1",[path,trend_path])
