"""Draw series comparisons from published, aligned animal observations."""
from pathlib import Path
import json
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from srcs.core.artifacts import artifact_path, sha256, write_json
from srcs.utils.visualization import displayed_weeks, filter_cohorts, settings
from srcs.visualization.scales import resolve


def render(folder, first, second, paired, *, viz=None, reference_week=None, destination=None):
    viz = settings() if viz is None else viz
    colors = dict(zip(("first", "second"), viz["series_colors"]))
    markers = viz["view_markers"]
    root = Path(destination) if destination is not None else Path(folder).parents[1] / "report"
    root.mkdir(parents=True,exist_ok=True)
    source=Path(folder)/"projection.csv"
    display=pd.read_csv(source)
    summaries=pd.read_csv(Path(folder)/"trajectories.csv")
    distances=pd.read_csv(Path(folder)/"distances.csv")
    mice=json.loads((Path(folder)/"identity.json").read_text(encoding="utf-8"))["mice"]
    weeks_to_plot = displayed_weeks(display.week.unique(), viz, reference_week)
    display = filter_cohorts(display[display.week.isin(weeks_to_plot)], viz)
    distances = filter_cohorts(distances[distances.week.isin(weeks_to_plot)], viz)
    summaries = filter_cohorts(summaries[summaries.week.isin(weeks_to_plot)], viz)
    paths, scales = [], []
    for (atlas,cohort),group in display.groupby(["atlas","cohort"]):
        weeks=list(group.week.unique())
        fig,axes=plt.subplots(1,len(weeks),squeeze=False,sharex=True,sharey=True,
                              figsize=(max(4,3*len(weeks)),4))
        for ax,week in zip(axes[0],weeks):
            part=group[group.week==week]
            for (side,mouse),samples in part.groupby(["side","mouse"]):
                color=colors[side]
                ax.plot(samples.pc1,samples.pc2,color=color,alpha=.35)
                for i,(_,run) in enumerate(samples.groupby("run",sort=False)):
                    ax.scatter(run.pc1,run.pc2,color=color,marker=markers[i % len(markers)],s=16)
            ax.set(title=week,xlabel="PC1")
        axes[0,0].set_ylabel("PC2")
        for title,color in ((first,colors["first"]),(second,colors["second"])):
            axes[0,0].scatter([],[],label=title,color=color)
        for i, run in enumerate(group.run.unique()):
            axes[0,0].scatter([], [], marker=markers[i % len(markers)], color="gray", label=run)
        axes[0,0].legend(fontsize=6)
        var=group.iloc[0]
        fig.suptitle(f"Common PCA ({var.variance_pc1:.1%}, {var.variance_pc2:.1%} explained)")
        path=root/f"{atlas}_{cohort}_projection.png"
        fig.savefig(path,dpi=viz["defaults"]["dpi"],bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
        summary=summaries[(summaries.atlas==atlas)&(summaries.cohort==cohort)]
        fig,axes=plt.subplots(2,2,figsize=(10,7))
        for ax,name,title in zip(axes.flat,["dispersion","self_distance","q","margin"],
                                 ["Dispersion","Within-mouse RMS","Discriminability","Separation margin"]):
            for side,label,color in (("first",first,colors["first"]),("second",second,colors["second"])):
                series=summary[summary.side==side].groupby("week",sort=False)[name].mean()
                ax.plot(weeks,series.reindex(weeks),color=color,label=label)
            ax.set_title(title)
        axes[0,0].legend(fontsize=6)
        fig.tight_layout()
        path=root/f"{atlas}_{cohort}_trajectories.png"
        fig.savefig(path,dpi=viz["defaults"]["dpi"],bbox_inches="tight")
        plt.close(fig)
        paths.append(path)
        dm=distances[(distances.atlas==atlas)&(distances.cohort==cohort)]
        scale = resolve(dm.distance, viz, "distance", anchor_zero=True)
        if scale:
            scales.append(scale.record(metric="distance", atlas=atlas, cohort=cohort, weeks=weeks))
        for (week,query,reference),part in dm.groupby(["week","query","reference"]) if scale else ():
            fig,axes=plt.subplots(1,2,figsize=(9,4))
            for ax,side in zip(axes,("first","second")):
                sample=part[part.side==side]
                if sample.empty:
                    ax.set_axis_off()
                    continue
                matrix=sample.pivot(index="mouse",columns="other",values="distance").reindex(
                    index=mice[side],columns=mice[side])
                im=ax.imshow(matrix, **scale.image)
                ax.set_xticks(range(len(matrix.columns)),matrix.columns,rotation=90,fontsize=7)
                ax.set_yticks(range(len(matrix)),matrix.index,fontsize=7)
                valid=matrix.notna().any(axis=1).to_numpy()
                ax.scatter(np.nanargmin(matrix.values[valid],axis=1),np.flatnonzero(valid),marker="x",color="cyan")
                ax.set(title=(first if side == "first" else second)+" | n="+str(valid.sum()),xlabel=reference,ylabel=query)
            fig.colorbar(im,ax=list(axes),label="RMS",extend=scale.extend)
            path=root/f"{atlas}_{cohort}_{week}_{query}_{reference}_rms.png"
            fig.savefig(path,dpi=viz["defaults"]["dpi"],bbox_inches="tight")
            plt.close(fig)
            paths.append(path)
        if paired:
            fig,ax=plt.subplots(figsize=(8,4))
            pairs=group.pivot(index=["week","run","mouse"],columns="side",values="q")
            for row in pairs.itertuples():
                ax.plot([0,1],[row.first,row.second],alpha=.3)
            ax.set_xticks([0,1],["First branch","Second branch"])
            ax.set(ylabel="Discriminability",ylim=(0,1),title="Paired branches | "+atlas)
            path=root/f"{atlas}_{cohort}_paired.png"
            fig.savefig(path,dpi=viz["defaults"]["dpi"],bbox_inches="tight")
            plt.close(fig)
            paths.append(path)
    scale_file = root / "color_scales.json"
    write_json(scale_file, scales)
    inputs = [source, Path(folder) / "trajectories.csv", Path(folder) / "distances.csv", Path(folder) / "identity.json"]
    write_json(root / "manifest.json", dict(first=first, second=second,
               source={artifact_path(p, root): sha256(p) for p in inputs},
               figures={artifact_path(p, root): sha256(p) for p in paths},
               color_scales={artifact_path(scale_file, root): sha256(scale_file)}))
