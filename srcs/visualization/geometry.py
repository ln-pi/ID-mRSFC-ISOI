"""Render saved geometry without refitting FC or changing cohort membership."""
import base64
import html
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd

from srcs.core.artifacts import sha256, write_json


def render(ctx):
    source = Path(ctx.output_path)/"results/geometry/v1"
    index = json.loads((source/"index.json").read_text(encoding="utf-8"))
    for path,expected in index["files"].items():
        if sha256(path) != expected:
            raise ValueError(f"Geometry product changed: {path}")
    target = Path(ctx.output_dir("report","geometry_figures"))
    tables = {p.stem: pd.read_csv(p) for p in map(Path,index["files"]) if p.suffix==".csv"}
    figures = []
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
            for run,marker in ((run_a,"o"),(run_b,"^")):
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
        vmax = dm.distance.max()
        for week,part in dm.groupby("week",sort=False):
            matrix = part.pivot(index="mouse",columns="other",values="distance").reindex(index=ctx.mice,columns=ctx.mice)
            fig,ax = plt.subplots(figsize=(5,4))
            im=ax.imshow(matrix,cmap=plt.get_cmap("inferno").with_extremes(bad="lightgray"),vmin=0,vmax=vmax)
            ax.set_xticks(range(len(matrix.columns)),matrix.columns,rotation=90)
            ax.set_yticks(range(len(matrix.index)),matrix.index)
            valid=matrix.notna().any(axis=1).to_numpy()
            ax.scatter(np.nanargmin(matrix.values[valid],axis=1),np.flatnonzero(valid),marker="x",color="cyan")
            ax.set(xlabel=run_a+" reference",ylabel=run_b+" query",title=week+" | RMS")
            fig.colorbar(im,ax=ax)
            save(fig,prefix+"_"+week+"_matrix")
        sub = select("subjects")
        sub = sub[(sub.metric=="rms")&(sub.rule=="full")]
        fig,axs = plt.subplots(1,2,figsize=(8,4))
        for row in sub.itertuples():
            axs[0].plot([0,1],[row.self_distance,row.nearest_distance],alpha=.35)
            axs[1].scatter(row.week,row.margin,marker="o" if row.query==run_a else "^",alpha=.5)
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
            from matplotlib.colors import TwoSlopeNorm
            limit=max(np.finfo(float).eps,float(signed.estimate.abs().max()))
            fig,axs=plt.subplots(1,len(weeks),figsize=(3.5*len(weeks),3.6),squeeze=False)
            for ax,week in zip(axs[0],weeks):
                sub=signed[signed.week==week]
                matrix=pd.DataFrame(np.nan,index=ctx.mice,columns=ctx.mice)
                for row in sub.itertuples():
                    matrix.loc[row.first,row.second]=matrix.loc[row.second,row.first]=row.estimate
                im=ax.imshow(matrix,cmap=plt.get_cmap("coolwarm").with_extremes(bad="lightgray"),
                             norm=TwoSlopeNorm(0,-limit,limit))
                ax.set_title(week)
            fig.colorbar(im,ax=list(axs[0]),label="Crossvalidated squared difference")
            save(fig,prefix+"_crossval")
    return figures


def publish_report(ctx, figures, inputs, skipped):
    folder = Path(ctx.output_dir("report"))
    entries, body = [], []
    with PdfPages(folder/"report.pdf") as pdf:
        for path in figures:
            path = Path(path)
            entries.append(dict(path=str(path),sha256=sha256(path)))
            encoded = base64.b64encode(path.read_bytes()).decode("ascii")
            title = html.escape(path.stem)
            body.append(f'<figure><img alt="{title}" src="data:image/png;base64,{encoded}"><figcaption>{title}</figcaption></figure>')
            fig,ax=plt.subplots(figsize=(11.7,8.3))
            ax.imshow(plt.imread(path))
            ax.axis("off")
            pdf.savefig(fig,bbox_inches="tight")
            plt.close(fig)
        if not figures:
            fig,ax=plt.subplots()
            ax.text(.1,.5,"No figures available")
            ax.axis("off")
            pdf.savefig(fig)
            plt.close(fig)
    (folder/"report.html").write_text('<!doctype html><meta charset="utf-8"><title>FC geometry</title>'
        '<style>body{font:16px system-ui;max-width:1100px;margin:auto}img{max-width:100%}</style>'
        '<h1>FC 高维几何分析</h1>'+''.join(body),encoding="utf-8")
    write_json(folder/"manifest.json",dict(figures=entries,inputs=inputs,skipped=skipped))
