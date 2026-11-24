"""Render saved geometry without refitting FC or changing cohort membership."""
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from srcs.core.artifacts import sha256, write_json
from srcs.utils.visualization import displayed_weeks, filter_cohorts
from srcs.visualization.scales import resolve


def render(ctx):
    source = Path(ctx.output_path)/"results/geometry/"
    index = json.loads((source/"index.json").read_text(encoding="utf-8"))
    for path,expected in index["files"].items():
        if sha256(path) != expected:
            raise ValueError(f"Geometry product changed: {path}")
    target = Path(ctx.output_dir("report","geometry_figures"))
    tables = {p.stem: pd.read_csv(p) for p in map(Path,index["files"]) if p.suffix==".csv"}
    viz = ctx.viz()
    weeks_to_plot = displayed_weeks(ctx.weeks, viz, ctx.section("descriptive").get("reference_week"))
    tables = {key: filter_cohorts(frame[frame.week.isin(weeks_to_plot)] if "week" in frame else frame, viz)
              for key, frame in tables.items()}
    figures, scales = [], []
    def save(fig,name):
        path = target/(name+".png")
        fig.savefig(path,dpi=ctx.render_dpi(),bbox_inches="tight")
        plt.close(fig)
        figures.append(path)
    for (atlas,cohort,run_a,run_b),points in tables["points"].groupby(["atlas","cohort","run_a","run_b"],sort=False):
        prefix = f"{atlas}_{cohort}_{run_a}_{run_b}"
        def select(key):
            f = tables[key]
            return f[(f.atlas==atlas)&(f.cohort==cohort)&(f.run_a==run_a)&(f.run_b==run_b)]
        if len(points)<2:
            continue
        points = points.copy()
        weeks = list(points.week.unique())
        fig,axs = plt.subplots(1,len(weeks),figsize=(3.5*len(weeks),3.6),sharex=True,sharey=True,squeeze=False)
        for ax,(week,part) in zip(axs[0],points.groupby("week",sort=False)):
            for _,mouse in part.groupby("mouse"):
                ax.plot(mouse.pc1,mouse.pc2,color=".7",linewidth=.8)
            markers = viz["view_markers"]
            for i, run in enumerate((run_a, run_b)):
                marker = markers[i % len(markers)]
                group = part[part.run==run]
                ax.scatter(group.pc1,group.pc2,marker=marker,label=run)
            ax.set_title(f"{week} (n={part.mouse.nunique()})")
            ax.set_xlabel("PC1")
        axs[0,0].set_ylabel("PC2")
        axs[0,0].legend()
        ratios=points.iloc[0][["variance_pc1","variance_pc2"]]
        fig.suptitle("Pooled PCA | variance "+", ".join(f"{v:.1%}" for v in ratios))
        save(fig,prefix+"_projection")
        distances = select("distances")
        dm = distances[(distances.metric=="rms") & (distances.query==run_b) & (distances.reference==run_a)]
        scale = resolve(dm.distance, viz, "distance", anchor_zero=True)
        if scale:
            scales.append(scale.record(atlas=atlas, cohort=cohort, query=run_b, reference=run_a,
                                       metric="distance", weeks=weeks))
        for week,part in dm.groupby("week",sort=False) if scale else ():
            matrix = part.pivot(index="mouse",columns="other",values="distance").reindex(index=ctx.mice,columns=ctx.mice)
            fig,ax = plt.subplots(figsize=(5,4))
            im=ax.imshow(matrix, **scale.image)
            ax.set_xticks(range(len(matrix.columns)),matrix.columns,rotation=90)
            ax.set_yticks(range(len(matrix.index)),matrix.index)
            valid=matrix.notna().any(axis=1).to_numpy()
            ax.scatter(np.nanargmin(matrix.values[valid],axis=1),np.flatnonzero(valid),marker="x",color="cyan")
            ax.set(xlabel=run_a+" reference",ylabel=run_b+" query",title=week+" | RMS")
            fig.colorbar(im,ax=ax,extend=scale.extend)
            save(fig,prefix+"_"+week+"_matrix")
        sub = select("subjects")
        sub = sub[(sub.metric=="rms")&(sub.rule=="full")]
        fig,axs = plt.subplots(1,2,figsize=(8,4))
        for row in sub.itertuples():
            axs[0].plot([0,1],[row.self_distance,row.nearest_distance],alpha=.35)
            axs[1].scatter(row.week,row.margin,marker=markers[int(row.query != run_a) % len(markers)],alpha=.5)
        axs[0].set_xticks([0,1],["Self","Nearest other"])
        axs[0].set_ylabel("RMS")
        axs[1].axhline(0,color=".5")
        axs[1].set_ylabel("Separation margin")
        save(fig,prefix+"_matching")
        fig,axs=plt.subplots(2,2,figsize=(10,7))
        disp=select("dispersion").groupby("week",sort=False).dispersion.mean()
        axs.flat[0].plot(disp.index,disp.values)
        axs.flat[0].set_title("Group dispersion")
        for ax,key,title in zip(list(axs.flat)[1:],["self_distance","q","margin"],
                                ["Within-mouse RMS","Discriminability","Separation margin"]):
            for _,group in sub.groupby("mouse"):
                line=group.groupby("week",sort=False)[key].mean().reindex(weeks)
                ax.plot(weeks,line.values,color=".7",linewidth=.6)
            ax.plot(weeks,sub.groupby("week")[key].mean().reindex(weeks),color="black")
            ax.set_title(title)
        fig.tight_layout()
        save(fig,prefix+"_trajectories")
        other=select("subjects")
        other=other[(other.metric=="rms") & (other.rule=="leave_one_out")]
        if not other.empty:
            fig,ax=plt.subplots(figsize=(7,4))
            for removed,part in other.groupby("removed"):
                series=part.groupby("week").q.mean().reindex(weeks)
                ax.plot(weeks,series,label=str(removed),alpha=.65)
            ax.set(title="Leave-one-mouse sensitivity",ylabel="Discriminability",ylim=(0,1))
            save(fig,prefix+"_sensitivity")
        signed=select("crossval")
        signed=signed[signed["mode"]=="identity"] if "mode" in signed else signed.iloc[:0]
        if len(signed) and signed.estimate.notna().any():
            scale = resolve(signed.estimate, viz, "signed_offset")
            scales.append(scale.record(atlas=atlas, cohort=cohort, metric="crossvalidated_difference", weeks=weeks))
            fig,axs=plt.subplots(1,len(weeks),figsize=(3.5*len(weeks),3.6),squeeze=False)
            for ax,week in zip(axs[0],weeks):
                sub=signed[signed.week==week]
                matrix=pd.DataFrame(np.nan,index=ctx.mice,columns=ctx.mice)
                for row in sub.itertuples():
                    matrix.loc[row.first,row.second]=matrix.loc[row.second,row.first]=row.estimate
                im=ax.imshow(matrix, **scale.image)
                ax.set_title(week)
            fig.colorbar(im,ax=list(axs[0]),label="Crossvalidated squared difference", extend=scale.extend)
            save(fig,prefix+"_crossval")
    write_json(target/"color_scales.json", scales)
    return figures



def report_section(ctx, source):
    index = json.loads(source.read_text(encoding="utf-8"))
    return render(ctx), [str(source), *index["files"]]
