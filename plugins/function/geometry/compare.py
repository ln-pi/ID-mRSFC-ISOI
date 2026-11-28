"""Compare completed series in their common labeled edge space."""
import ast
import json
from inspect import getsource
from pathlib import Path

import numpy as np
import pandas as pd

from srcs.core.artifacts import artifact_entries, publish, sha256, write_json
from srcs.core.identity import code_closure, digest
from srcs.core.staging import staging_directory
from srcs.utils.config import PROJECT_ROOT
from srcs.utils.storage import save_table
from .discriminability import score
from srcs.analysis.pairwise import week_positions


def saved_positions(root, weeks):
    marker = root / "metadata/stages/inference.json"
    if marker.is_file():
        options = json.loads(marker.read_text(encoding="utf-8"))["options"].get("pairwise", {})
    else:
        metadata = root / "metadata/analysis_config.json"
        options = json.loads(metadata.read_text(encoding="utf-8"))["config"].get("pairwise", {})
    return week_positions(weeks, options.get("week_positions"))


def load_products(root):
    root = Path(root)
    identity = json.loads((root/"metadata/data_identity.json").read_text(encoding="utf-8"))
    hashes = {}
    for name in ("geometry", "inference"):
        path = root/"results"/name/"index.json"
        index = json.loads(path.read_text(encoding="utf-8"))
        files = [*artifact_entries(index["files"], path.parent),
                 *artifact_entries(index.get("products", []), path.parent)]
        checksums = {}
        for item in files:
            product = item["path"]
            if sha256(product) != item["sha256"]:
                raise ValueError(f"Changed product: {product}")
            role = product.relative_to(path.parent).as_posix() if product.is_relative_to(path.parent) else product.name
            checksums[role] = item["sha256"]
        hashes[name] = digest(checksums)
    folder=root/"results/geometry/"
    basis=json.loads((folder/"basis.json").read_text(encoding="utf-8"))
    points=pd.read_csv(folder/"points.csv").drop_duplicates(["atlas","cohort","week","run","mouse"])
    subjects=pd.read_csv(folder/"subjects.csv")
    scores=subjects[(subjects.rule=="full")&(subjects.metric=="rms")].groupby(
        ["atlas","cohort","week","query","mouse"],as_index=False)[["q","margin"]].mean()
    points=points.merge(scores.rename(columns={"query":"run"}),how="left",
                        on=["atlas","cohort","week","run","mouse"],validate="one_to_one")
    points["point"]=points.point.map(lambda s: np.array(ast.literal_eval(s)))
    return identity, hashes, basis, points


def compare(first, second, *, destination, configs, seed=None, repeats=None):
    from srcs.utils.visualization import settings
    from srcs.workflow.plugins import function_options
    display_config = configs[0]
    viz = settings(display_config)
    reference = display_config.get("descriptive", {}).get("reference_week")
    options = function_options(display_config, "geometry", "compare")
    if options != function_options(configs[1], "geometry", "compare"):
        raise ValueError("Geometry comparison requires matching comparison options")
    seed = options.get("seed", 42) if seed is None else seed
    repeats = options.get("repeats", 200) if repeats is None else repeats
    left,right = map(lambda p:Path(p).resolve(),(first,second))
    a,ha,ba,x=load_products(left)
    b,hb,bb,y=load_products(right)
    if ba!=bb or a["scope"]["weeks"]!=b["scope"]["weeks"] or a["scope"]["runs"]!=b["scope"]["runs"]:
        raise ValueError("Comparison requires identical partitions, ordered ROI/edges, weeks and runs")
    positions = saved_positions(left, a["scope"]["weeks"])
    if positions != saved_positions(right, b["scope"]["weeks"]):
        raise ValueError("Comparison requires identical numerical week positions")
    paired=a["scope"]["series"]==b["scope"]["series"]
    if paired == (a["scope"]["branch"]==b["scope"]["branch"]):
        raise ValueError("Compare two series at one branch or two branches within one series")
    keys=["atlas","cohort","week","run","mouse"]
    if paired and (set(a["scope"]["mice"])!=set(b["scope"]["mice"]) or
                   set(map(tuple,x[keys].values))!=set(map(tuple,y[keys].values))):
        raise ValueError("Paired branches differ in mouse or valid acquisition membership")
    from srcs.core.runtime import exclusive_run

    numerical_code = [getsource(func) for func in (load_products, saved_positions, calculate_comparison)]
    numerical_code.append(code_closure(PROJECT_ROOT, f"{__package__}.discriminability"))
    payload = dict(schema="compare_identity", a=digest(a), b=digest(b), inputs=[ha, hb],
                   mode="paired" if paired else "independent",
                   mice={"first": a["scope"]["mice"], "second": b["scope"]["mice"]},
                   code=digest(numerical_code),
                   positions=positions, seed=seed, repeats=repeats, candidate_limit=options.get("candidate_limit", 9))
    target = Path(destination).resolve()
    numeric = target / "results/geometry"
    with exclusive_run(target / "metadata/geometry.lock"):
        save_batch_numeric(numeric, x, y, left, right, a, b, paired, positions, payload)
        render_comparison(numeric, a, b, paired, viz, reference, report_destination=target / "report/geometry")
    return numeric


def save_batch_numeric(folder, x, y, left, right, a, b, paired, positions, payload):
    """Publish one verified numerical result, independent of report settings."""
    from srcs.workflow.products import verified_index

    marker = folder / "identity.json"
    try:
        published_identity = json.loads(marker.read_text(encoding="utf-8"))
        if digest(published_identity) == digest(payload):
            verified_index(folder / "index.json")
            return
    except (OSError, ValueError, KeyError, TypeError):
        pass
    folder.parent.mkdir(parents=True, exist_ok=True)
    with staging_directory(folder.parent, "geometry") as staging:
        staged = staging / "results"
        products = calculate_comparison(staged, x, y, left, right, a, b, paired, positions,
                                        payload["seed"], payload["repeats"], payload["candidate_limit"], payload)
        publish(staged / "index.json", "geometry_comparison", "geometry_comparison",
                [*products, staged / "identity.json"])
        old = Path(staging) / "previous"
        if folder.exists():
            folder.rename(old)
        try:
            staged.rename(folder)
        except OSError:
            if old.exists():
                old.rename(folder)
            raise


def render_comparison(folder, a, b, paired, viz, reference, *, report_destination):
    first_label, second_label = a["scope"]["series"], b["scope"]["series"]
    if first_label == second_label:
        first_label += "/" + a["scope"]["branch"]
        second_label += "/" + b["scope"]["branch"]
    from .comparison_figures import render
    render(folder, first_label, second_label, paired, viz=viz, reference_week=reference,
           destination=report_destination)


def calculate_comparison(folder, x, y, left, right, a, b, paired, positions, seed, repeats, candidate_limit, payload):
    keys = ["atlas", "cohort", "week", "run", "mouse"]
    groups=["atlas","cohort","week","run"]
    common_n=pd.concat([x.groupby(groups).mouse.nunique(),y.groupby(groups).mouse.nunique()],
                       axis=1,join="inner").min(axis=1)
    cross_matrices=[]
    for side,root,frame in (("first",left,x),("second",right,y)):
        matrix=pd.read_csv(root/"results/geometry/distances.csv")
        matrix=matrix[(matrix.metric=="rms")&(matrix.query!=matrix.reference)]
        cross_matrices.append(matrix.assign(side=side))
        own=matrix[matrix.mouse==matrix.other].groupby(["atlas","cohort","week","query","mouse"]).distance.mean()
        frame["self_distance"]=own.reindex(pd.MultiIndex.from_frame(frame[keys])).to_numpy()
        scores=[]
        for key,part in matrix.groupby(["atlas","cohort","week","query","reference","mouse"]):
            if key[:4] not in common_n:
                continue
            own=np.flatnonzero(part.other.to_numpy()==key[-1])
            if len(own)!=1:
                raise ValueError(f"Missing or duplicated self reference: {key}")
            others=np.flatnonzero(part.other.to_numpy()!=key[-1])
            count=min(candidate_limit,int(common_n[key[:4]])-1,len(others))
            subsampler=np.random.default_rng(int(digest([seed,*key])[:16],16))
            selected=subsampler.choice(others,count,replace=False)
            result=score(part.distance.to_numpy(),own[0],selected)
            scores.append(dict(zip(["atlas","cohort","week","run","reference","mouse"],key),
                               q=result["q"],margin=result["margin"]))
        if scores:
            updated=pd.DataFrame(scores).groupby(keys).agg({"q":"mean","margin":"mean"})
            aligned=updated.reindex(pd.MultiIndex.from_frame(frame[keys]))
            frame[["q","margin"]]=aligned[["q","margin"]].to_numpy()
    rng=np.random.default_rng(seed)
    ids_a=np.array(sorted(a["scope"]["mice"]))
    ids_b=np.array(sorted(b["scope"]["mice"]))
    draws_a=rng.choice(ids_a,(repeats,len(ids_a)))
    draws_b=draws_a if paired else rng.choice(ids_b,(repeats,len(ids_b)))
    rows={k:[] for k in ("centroid_gap","dispersion_gap","q_gap","group_time")}
    def measures(points,scores,indices):
        if not len(indices):
            return None
        sample=points[indices]
        return (sample.mean(axis=0), np.var(sample,axis=0,ddof=1).mean() if len(sample)>1 else np.nan,
                *pd.DataFrame(scores[indices]).mean().to_numpy())
    def gap(ma,mb):
        if ma is None or mb is None:
            return [np.nan]*4
        return [np.sqrt(np.mean((ma[0]-mb[0])**2)),mb[1]-ma[1],mb[2]-ma[2],mb[3]-ma[3]]
    def intervals(values):
        finite=values[np.isfinite(values)]
        return np.quantile(finite,[.05,.95]).tolist() if len(finite) else [np.nan,np.nan]
    group_b={key:part for key,part in y.groupby(groups)}
    trajectories={}
    for key,p in x.groupby(groups):
        if key not in group_b:
            continue
        sides=[]
        for part,draws in ((p,draws_a),(group_b[key],draws_b)):
            pts=np.stack(part.point)
            scores=part[["q","margin"]].to_numpy()
            lookup={m:i for i,m in enumerate(part.mouse)}
            sides.append([measures(pts,scores,np.arange(len(pts)))]+[
                measures(pts,scores,[lookup[m] for m in sample if m in lookup]) for sample in draws])
        values=np.array([gap(ma,mb) for ma,mb in zip(*sides)])
        for j,name in enumerate(("centroid_gap","dispersion_gap","q_gap")):
            low,high=intervals(values[1:,j])
            rows[name].append(dict(zip(groups,key),estimate=values[0,j],low=low,high=high,
                                   interval="sensitivity_range"))
            trajectories.setdefault((key[0],key[1],key[3],name),[]).append((key[2],values[:,j]))
        margin_low,margin_high=intervals(values[1:,3])
        rows["q_gap"][-1].update(margin_gap=values[0,3],margin_low=margin_low,margin_high=margin_high,
                                 n_foreign=min(candidate_limit,int(common_n[key])-1))
        trajectories.setdefault((key[0],key[1],key[3],"margin_gap"),[]).append((key[2],values[:,3]))
    for key,items in trajectories.items():
        labels,values=zip(*items)
        t=np.array([positions[w] for w in labels])
        slopes=[]
        for trajectory in np.stack(values,axis=1):
            keep=np.isfinite(trajectory)
            slopes.append(np.polyfit(t[keep],trajectory[keep],1)[0] if keep.sum()>=2 else np.nan)
        low,high=intervals(np.array(slopes[1:]))
        rows["group_time"].append(dict(atlas=key[0],cohort=key[1],run=key[2],endpoint=key[3],
                                      slope=slopes[0],low=low,high=high,interval="sensitivity_range"))
    folder.mkdir(parents=True,exist_ok=True)
    paths=[]
    from sklearn.decomposition import PCA
    display=pd.concat([x.assign(side="first"),y.assign(side="second")],ignore_index=True)
    for _,part in display.groupby(["atlas","cohort"]):
        if len(part)<2:
            continue
        model=PCA(n_components=2,svd_solver="full")
        display.loc[part.index,["pc1","pc2"]]=model.fit_transform(np.stack(part.point))
        display.loc[part.index,"variance_pc1"]=model.explained_variance_ratio_[0]
        display.loc[part.index,"variance_pc2"]=model.explained_variance_ratio_[1]
    projection=folder/"projection.csv"
    matrices=folder/"distances.csv"
    save_table(matrices,pd.concat(cross_matrices),index=False)
    paths.append(matrices)
    views=[]
    for key,part in display.groupby(["side",*groups]):
        pts=np.stack(part.point)
        views.append(dict(zip(["side",*groups],key),n=len(part),
                          dispersion=np.var(pts,axis=0,ddof=1).mean() if len(part)>1 else np.nan,
                          self_distance=part.self_distance.mean(),q=part.q.mean(),margin=part.margin.mean()))
    summary=folder/"trajectories.csv"
    save_table(summary,pd.DataFrame(views),index=False)
    paths.append(summary)
    save_table(projection,display.drop(columns="point"),index=False)
    paths.append(projection)
    for name,table in rows.items():
        path=folder/(name+".csv")
        save_table(path,pd.DataFrame(table),index=False)
        paths.append(path)
    write_json(folder/"identity.json",payload)
    return paths
