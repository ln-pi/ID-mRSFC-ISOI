import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
from pathlib import Path
import pandas as pd
from scipy.io import savemat
from tqdm import tqdm

from srcs.utils.logger import setup_logger
from srcs.analysis.statistics import individual_identifiability
from srcs.utils.io import save_table, save_json


logger = setup_logger()








def one_key_save(inpath, expath):
    """Alternatively, if HD drop problem happens, use this function."""
    from srcs.preprocessing.readers import One_Key_Graph_Reader, One_Key_AI_Reader
    raw_r_dics = One_Key_Graph_Reader(inpath, nan=False)
    raw_ai_signals = One_Key_AI_Reader(inpath)

    # NOTE: returns are packed dicts
    logger.info(f"R value dics are in format {type(raw_r_dics)}, keys are: {raw_r_dics.keys()}")
    logger.info(f"AI value dics are in format {type(raw_ai_signals)}, in shape {raw_ai_signals.shape}")





def export_individual_identifiability(ctx, stat, analysis_meta):
    """Use full similarity matrices: the aggregate alone cannot recover ranks."""
    rows = []
    excluded = set(ctx.section('analysis').get('statistics_exclude', []))
    atlases = ['_'.join(a) for a in ctx.atlases if '_'.join(a) not in excluded]
    mice = list(analysis_meta['mice'])
    if len(mice) != len(set(mice)):
        raise ValueError('Duplicate mouse IDs in analysis metadata.')
    for atlas in atlases:
        for condition in stat.columns:
            path = Path(ctx.output_dir('Legacy', create=False))/'identification'/str(condition)/'Numeric'/f'identification_{atlas}.pkl'
            with path.open('rb') as stream:
                payload = pickle.load(stream)
            if payload.get('analysis') != analysis_meta or payload.get('mice') != mice:
                raise ValueError(f'Similarity metadata/cohort differs from aggregate: {path}')
            metrics = individual_identifiability(payload['similarity'])
            if len(metrics['individual_similarity']) != len(mice):
                raise ValueError(f'Similarity size differs from mouse labels: {path}')
            cell = stat.at[atlas, condition]
            for actual, expected in [(metrics['individual_similarity'], cell[1]),
                                     (metrics['other_similarity'], cell[0])]:
                if np.shape(actual) != np.shape(expected) or not np.allclose(actual, expected, rtol=1e-9, atol=1e-10):
                    raise ValueError(f'Similarity values differ from aggregate: {path}')
            for i, mouse in enumerate(mice):
                rows.append(dict(mouse=mouse, atlas=atlas, analysis_condition=str(condition),
                                 n_mice=len(mice), **{k:float(v[i]) for k,v in metrics.items()}))
    columns = ['mouse','atlas','analysis_condition','n_mice','individual_similarity','other_similarity',
               'self_match_advantage','self_match_percentile','row_self_match_percentile','col_self_match_percentile']
    frame = pd.DataFrame(rows, columns=columns)
    save_table(Path(ctx.output_dir('Statistics','tables','legacy'))/'individual_identifiability.csv',frame,index=False)
    save_json(Path(ctx.output_dir('Statistics','tables','legacy'))/'individual_identifiability_metadata.json',dict(
        analysis=analysis_meta, status='complete' if rows else 'empty_selection',
        inference='descriptive_only', advantage='self minus mean of row/column nonself similarities',
        percentile='mean row/column fraction of nonself candidates beaten; exact tie=0.5',
        percentile_range=[0,1], similarity_scale='saved_03_similarity',
        caveats=['Shared nonself references induce dependence.',
                 'Candidate count/composition affects interpretation; not a disease severity score.']))
    options = ctx.render_options()
    if options['save'] or options['show']:
        for atlas in atlases:
            for metric in ['self_match_advantage','self_match_percentile']:
                fig, ax = plt.subplots(figsize=(max(6,len(stat.columns)*.65),4.5))
                table=frame[frame.atlas==atlas].pivot(index='mouse',columns='analysis_condition',values=metric)
                table=table.reindex(index=mice,columns=[str(c) for c in stat.columns])
                x=np.arange(len(table.columns))
                for mouse, values in table.iterrows():
                    ax.plot(x,values,marker='o',ms=3,lw=1,alpha=.6,label=mouse)
                ax.set_xticks(x,table.columns,rotation=45 if len(x)>5 else 0)
                ax.set(xlabel='Analysis condition',ylabel=metric.replace('_',' '),title=atlas)
                ax.axhline(0 if metric=='self_match_advantage' else .5,color='gray',ls=':',lw=1)
                if metric=='self_match_percentile': ax.set_ylim(-.025,1.025)
                ax.legend(bbox_to_anchor=(1.02,1),loc='upper left',fontsize=7)
                fig.tight_layout()
                if options['save']:
                    fig.savefig(Path(ctx.output_dir('Statistics','figures','legacy','individual_identifiability'))/f'{metric}_{atlas}.png',dpi=ctx.figure_dpi())
                if options['show']: plt.show()
                plt.close(fig)
    return frame

from srcs.utils.io import mkdir, save_as_npy, save_as_mat, Save_Variable, save_tensor
