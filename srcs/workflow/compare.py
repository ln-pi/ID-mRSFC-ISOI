"""Compare completed series in their common labeled edge space."""
import ast
import json
from pathlib import Path

import numpy as np
import pandas as pd

from srcs.core.artifacts import publish, sha256, write_json
from srcs.core.identity import code_closure, digest
from srcs.utils.config import PROJECT_ROOT
from srcs.utils.results import save_table


def load_products(root):
    root = Path(root)
    identity = json.loads((root/"metadata/data_identity.json").read_text(encoding="utf-8"))
    hashes = {}
    for name in ("geometry/v1","inference/v1","descriptive/v2"):
        path = root/"results"/name/"index.json"
        index = json.loads(path.read_text(encoding="utf-8"))
        files = index["files"]
        if isinstance(files,list):
            files = {v["path"]:v["sha256"] for v in files+index.get("products",[])}
        for item, expected in files.items():
            if sha256(item) != expected:
                raise ValueError(f"Changed product: {item}")
        hashes[name] = sha256(path)
    folder=root/"results/geometry/v1"
    basis=json.loads((folder/"basis.json").read_text(encoding="utf-8"))
    points=pd.read_csv(folder/"points.csv").drop_duplicates(["atlas","cohort","week","run","mouse"])
    subjects=pd.read_csv(folder/"subjects.csv")
    scores=subjects[(subjects.rule=="full")&(subjects.metric=="rms")].groupby(
        ["atlas","cohort","week","query","mouse"],as_index=False)[["q","margin"]].mean()
    points=points.merge(scores.rename(columns={"query":"run"}),how="left",
                        on=["atlas","cohort","week","run","mouse"],validate="one_to_one")
    points["point"]=points.point.map(lambda s: np.array(ast.literal_eval(s)))
    return identity, hashes, basis, points


def compare(first,second,output_root=None,seed=42,repeats=200):
    left,right = map(lambda p:Path(p).resolve(),(first,second))
    a,ha,ba,x=load_products(left)
    b,hb,bb,y=load_products(right)
    if ba!=bb or a["scope"]["weeks"]!=b["scope"]["weeks"] or a["scope"]["runs"]!=b["scope"]["runs"]:
        raise ValueError("Comparison requires identical partitions, ordered ROI/edges, weeks and runs")
    paired=a["scope"]["series"]==b["scope"]["series"]
    if paired == (a["scope"]["branch"]==b["scope"]["branch"]):
        raise ValueError("Compare two series at one branch or two branches within one series")
    keys=["atlas","cohort","week","run","mouse"]
    if paired and (set(a["scope"]["mice"])!=set(b["scope"]["mice"]) or
                   set(map(tuple,x[keys].values))!=set(map(tuple,y[keys].values))):
        raise ValueError("Paired branches differ in mouse or valid acquisition membership")
    payload=dict(schema="compare_identity_v1",a=digest(a),b=digest(b),inputs=[ha,hb],mode="paired" if paired else "independent",
                 code=digest(code_closure(PROJECT_ROOT,"srcs.workflow.compare")),seed=seed,repeats=repeats)
    folder=(Path(output_root) if output_root else left.parent)/"_compare"/(left.name+"__vs__"+right.name+"_"+digest(payload)[:8])/"results/compare/v1"
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
        return sample.mean(axis=0), np.var(sample,axis=0,ddof=1).mean() if len(sample)>1 else np.nan, scores[indices].mean()
    def gap(ma,mb):
        if ma is None or mb is None:
            return [np.nan]*3
        return [np.sqrt(np.mean((ma[0]-mb[0])**2)),mb[1]-ma[1],mb[2]-ma[2]]
    def intervals(values):
        finite=values[np.isfinite(values)]
        return np.quantile(finite,[.05,.95]).tolist() if len(finite) else [np.nan,np.nan]
    groups=["atlas","cohort","week","run"]
    group_b={key:part for key,part in y.groupby(groups)}
    trajectories={}
    for key,p in x.groupby(groups):
        if key not in group_b:
            continue
        sides=[]
        for part,draws in ((p,draws_a),(group_b[key],draws_b)):
            pts=np.stack(part.point)
            scores=part.q.to_numpy()
            lookup={m:i for i,m in enumerate(part.mouse)}
            sides.append([measures(pts,scores,np.arange(len(pts)))]+[
                measures(pts,scores,[lookup[m] for m in sample if m in lookup]) for sample in draws])
        values=np.array([gap(ma,mb) for ma,mb in zip(*sides)])
        for j,name in enumerate(("centroid_gap","dispersion_gap","q_gap")):
            low,high=intervals(values[1:,j])
            rows[name].append(dict(zip(groups,key),estimate=values[0,j],low=low,high=high,
                                   interval="sensitivity_range"))
            trajectories.setdefault((key[0],key[1],key[3],name),[]).append((key[2],values[:,j]))
    for key,items in trajectories.items():
        labels,values=zip(*items)
        cfg=json.loads((left/"metadata/analysis_config.json").read_text(encoding="utf-8"))["config"]
        positions=cfg.get("pairwise",{}).get("week_positions",{})
        t=np.array([positions[w] if w in positions else float(w.rstrip("w")) for w in labels])
        slopes=[]
        for trajectory in np.stack(values,axis=1):
            keep=np.isfinite(trajectory)
            slopes.append(np.polyfit(t[keep],trajectory[keep],1)[0] if keep.sum()>=2 else np.nan)
        low,high=intervals(np.array(slopes[1:]))
        rows["group_time"].append(dict(atlas=key[0],cohort=key[1],run=key[2],endpoint=key[3],
                                      slope=slopes[0],low=low,high=high,interval="sensitivity_range"))
    folder.mkdir(parents=True,exist_ok=True)
    paths=[]
    for name,table in rows.items():
        path=folder/(name+".csv")
        save_table(path,pd.DataFrame(table),index=False)
        paths.append(path)
    write_json(folder/"identity.json",payload)
    publish(folder/"index.json","comparison_v1","comparison_v1",paths+[folder/"identity.json"])
    return folder
