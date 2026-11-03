"""Publish per-week geometry, directed matching and reproducible differences."""
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from srcs.analysis.geometry.crossval import crossvalidated
from srcs.analysis.geometry.distances import distances, dispersion
from srcs.analysis.geometry.discriminability import score
from srcs.core.artifacts import publish, write_json
from srcs.core.identity import digest
from srcs.utils.edge_store import load_store
from srcs.utils.storage import save_table


from srcs.core.runtime import numerical_stage

OPTION_KEYS = ("geometry",)


@numerical_stage("geometry", OPTION_KEYS)
def run(ctx):
    folder = Path(ctx.output_dir("results", "geometry", "v1"))
    cfg = ctx.section("geometry")
    k = cfg.get("matched_candidates", 9)
    if not isinstance(k, int) or k < 1:
        raise ValueError("matched_candidates must be a positive integer")
    matrices, subjects, summaries, differences, points, displacements = [], [], [], [], [], []
    bases = {}
    for atlas in map("_".join, ctx.atlases):
        store = load_store(ctx, atlas)
        bases[atlas] = dict(roi_labels=list(store.roi_labels), edge_index=store.edge_index.tolist())
        for a,b in combinations(range(len(ctx.runs)), 2):
            valid = store.valid_record[:,:,a] & store.valid_record[:,:,b]
            for cohort in ("available_week", "fixed_window"):
                for w,week in enumerate(ctx.weeks):
                    keep = valid[:,w] if cohort == "available_week" else valid.all(axis=1)
                    ids = np.array(ctx.mice)[keep]
                    fields = dict(atlas=atlas, week=week, cohort=cohort, n=len(ids),
                                  run_a=ctx.runs[a], run_b=ctx.runs[b], edges=store.F_z.shape[-1])
                    x,y = (store.F_z[keep,w,r,:] for r in (a,b))
                    for r,values in ((a,x),(b,y)):
                        centroid_gap = np.sqrt(np.mean((x.mean(axis=0)-y.mean(axis=0))**2)) if len(ids) else np.nan
                        summaries.append(dict(fields, run=ctx.runs[r], dispersion=dispersion(values),
                                              run_centroid_gap=centroid_gap))
                        same_rms,same_shape=distances(values,values)
                        for mouse, point in zip(ids, values):
                            points.append(dict(fields, run=ctx.runs[r], mouse=mouse, point=point.tolist()))
                        for i,mouse in enumerate(ids):
                            for j,other in enumerate(ids):
                                for metric,matrix in (("rms",same_rms),("correlation",same_shape)):
                                    matrices.append(dict(fields,query=ctx.runs[r],reference=ctx.runs[r],
                                                         mouse=mouse,other=other,metric=metric,distance=matrix[i,j]))
                        for previous in range(w):
                            both=keep & valid[:,previous]
                            earlier=store.F_z[both,previous,r]
                            later=store.F_z[both,w,r]
                            rms,shape=distances(later,earlier)
                            for i,mouse in enumerate(np.array(ctx.mice)[both]):
                                displacements.append(dict(fields,run=ctx.runs[r],mouse=mouse,
                                    previous_week=ctx.weeks[previous],rms=rms[i,i],correlation=shape[i,i]))
                    for mode in ("identity", "shrinkage", "animal_holdout"):
                        cv, counts = crossvalidated(x, y, mode)
                        for i,j in combinations(range(len(ids)),2):
                            differences.append(dict(fields, first=ids[i], second=ids[j], mode=mode,
                                                    estimate=cv[i,j], n_noise=counts[i,j]))
                    for query,reference,rm in ((b,a,distances(y,x)),(a,b,distances(x,y))):
                        for metric,matrix in zip(("rms","correlation"),rm):
                            for i,mouse in enumerate(ids):
                                base = dict(fields, query=ctx.runs[query], reference=ctx.runs[reference],
                                            mouse=mouse, metric=metric)
                                for j,other in enumerate(ids):
                                    matrices.append(dict(base, other=other, distance=matrix[i,j]))
                                others = np.flatnonzero((np.arange(len(ids)) != i) & np.isfinite(matrix[i]))
                                seed = int(digest([cfg.get("seed",42),week,query,reference,str(mouse)])[:16],16)
                                rng = np.random.default_rng(seed)
                                matched = rng.choice(others,min(k,len(others)),replace=False)
                                pools = [("full","",others), ("matched","",matched)]
                                pools.extend(("leave_one_out",ids[j],others[others!=j]) for j in others)
                                for rule,removed,pool in pools:
                                    result = score(matrix[i],i,pool)
                                    result["nearest"] = "|".join(ids[result.pop("nearest_indices")])
                                    subjects.append(dict(base, rule=rule, removed=removed,
                                                         candidate_ids="|".join(ids[pool]), **result))
    basis_path = folder/"basis.json"
    write_json(basis_path,bases)
    files = [basis_path]
    for name,rows in (("distances",matrices),("subjects",subjects),("dispersion",summaries),
                      ("crossval",differences),("points",points),("displacement",displacements)):
        path = folder / (name + ".csv")
        columns = {"subjects": ["atlas","week","cohort","run_a","run_b","metric","rule","mouse","removed","q","margin"],
                   "points":["atlas","cohort","run_a","run_b","week","run","mouse","point"],
                   "distances":["atlas","cohort","run_a","run_b","week","metric","query","reference","mouse","other","distance"]}
        frame = pd.DataFrame(rows) if rows else pd.DataFrame(columns=columns.get(name,["atlas","status"]))
        if name=="points":
            for _,group in frame.groupby(["atlas","cohort","run_a","run_b"]):
                if len(group)<2:
                    continue
                projection=PCA(n_components=2,svd_solver="full")
                coords=projection.fit_transform(np.stack(group.point))
                frame.loc[group.index,["pc1","pc2"]]=coords
                frame.loc[group.index,"variance_pc1"]=projection.explained_variance_ratio_[0]
                frame.loc[group.index,"variance_pc2"]=projection.explained_variance_ratio_[1]
        save_table(path, frame, index=False)
        files.append(path)
    publish(folder/"index.json","geometry_metrics_v1","geometry_metrics_v1",files)
