"""Descriptive trajectory views without hidden statistical calculations."""
import matplotlib.pyplot as plt


def plot_trajectories(frame,path=None):
    atlases=list(frame.atlas.unique()) if 'atlas' in frame and len(frame) else ['No valid data']
    fig,axes=plt.subplots(len(atlases),1,figsize=(10,3.5*len(atlases)),squeeze=False)
    for atlas,ax in zip(atlases,axes[:,0]):
        if len(frame):
            for mouse,group in frame[frame.atlas==atlas].groupby('mouse',sort=False):
                ax.plot(group.week_index,group.distance_from_start_z,marker='o',label=mouse)
            ticks=frame[['week_index','week']].drop_duplicates().sort_values('week_index')
            ax.set_xticks(ticks.week_index,ticks.week)
            ax.legend(ncol=5,fontsize=8)
        ax.set_title(atlas);ax.set_ylabel('RMS change of FC z edges');ax.set_xlabel('Observation week')
    fig.tight_layout()
    if path is not None:
        fig.savefig(path,dpi=180);plt.close(fig)
    return fig
