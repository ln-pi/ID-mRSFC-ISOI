"""Draw series comparisons from published, aligned animal observations."""
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

from srcs.core.artifacts import sha256, write_json


def render(folder, first, second, paired):
    root=Path(folder).parents[2]/"report"
    root.mkdir(parents=True,exist_ok=True)
    source=Path(folder)/"projection.csv"
    display=pd.read_csv(source)
    paths=[]
    for (atlas,cohort),group in display.groupby(["atlas","cohort"]):
        weeks=list(group.week.unique())
        fig,axes=plt.subplots(1,len(weeks),squeeze=False,sharex=True,sharey=True,
                              figsize=(max(4,3*len(weeks)),4))
        for ax,week in zip(axes[0],weeks):
            part=group[group.week==week]
            for (side,mouse),samples in part.groupby(["side","mouse"]):
                color="#356caa" if side=="first" else "#bd502f"
                ax.plot(samples.pc1,samples.pc2,color=color,alpha=.35)
                for i,(_,run) in enumerate(samples.groupby("run",sort=False)):
                    ax.scatter(run.pc1,run.pc2,color=color,marker=("o","^","s")[i%3],s=16)
            ax.set(title=week,xlabel="PC1")
        axes[0,0].set_ylabel("PC2")
        for title,color in ((first,"#356caa"),(second,"#bd502f")):
            axes[0,0].scatter([],[],label=title,color=color)
        axes[0,0].legend(fontsize=6)
        var=group.iloc[0]
        fig.suptitle(f"Common PCA ({var.variance_pc1:.1%}, {var.variance_pc2:.1%} explained)")
        path=root/f"{atlas}_{cohort}_projection.png"
        fig.savefig(path,bbox_inches="tight")
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
            fig.savefig(path,bbox_inches="tight")
            plt.close(fig)
            paths.append(path)
    write_json(root/"manifest.json",dict(first=first,second=second,source={str(source):sha256(source)},
               figures={str(p):sha256(p) for p in paths}))
