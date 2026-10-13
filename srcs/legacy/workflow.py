"""Optional historical half-group analysis, isolated from primary results."""
import numpy as np
import pandas as pd
from pathlib import Path
from srcs.utils.pipeline import load_records
from srcs.utils.io import save_json,save_pickle
from srcs.utils.dataset import axis_frames,analysis_plan
from srcs.analysis.longitudinal import weekly_representatives,select_cohort
from srcs.analysis.identification import identify
from srcs.legacy.averaging import build_half_dataframe
from srcs.legacy.identification import Edgewise_Contribution
from srcs.utils.progress import progress_items


def run_legacy(ctx):
    records,audit,_=load_records(ctx)
    weekly,_=weekly_representatives(ctx,records);cohort=select_cohort(ctx,weekly)
    mice,weeks=cohort['mice'],cohort['weeks']
    folder=Path(ctx.output_dir('Legacy'))
    results={}
    if len(mice)<2:
        save_json(folder/'status.json',dict(status='unavailable',reason='insufficient complete cohort',cohort=cohort));return
    for atlas in ['_'.join(a) for a in ctx.atlases]:
        tables={}
        for week in weeks:
            frame=pd.DataFrame(index=mice,columns=ctx.runs,dtype=object)
            for m in mice:
                for r in ctx.runs:frame.at[m,r]=records[atlas,week,m,r]
            tables[week]=frame
        for label,condition in progress_items(analysis_plan(ctx,weeks).items(), phase=f'{atlas}: historical comparisons'):
            selected_weeks,selected_runs=condition['weeks'],condition['runs']
            frames=axis_frames(tables,ctx.split_axis,mice,selected_weeks,selected_runs)
            try:
                groups=ctx.group_indices(selected_weeks,selected_runs)
                repeat=ctx.section('ID').get('repeat_id')
                ids=list(frames)
                if isinstance(repeat,str):repeat=ids.index(repeat)+1
                halves=build_half_dataframe(list(frames.values()),group_mode=ctx.group_mode,
                    repeat_session=repeat,seed=ctx.session_seed,mouse_ids=mice)
                iu=np.triu_indices(len(halves.iloc[0,0]),1)
                vectors=[np.stack([np.asarray(halves.at[part,m])[iu] for m in mice]) for part in ['1st','2nd']]
                result=identify(*vectors,mice,n_iter=ctx.section('ID').get('n_iter',1000) if ctx.section('ID').get('shuffle',True) else 0,seed=ctx.session_seed)
                trial_ids=list(next(iter(frames.values())).columns)
                if len(trial_ids)<2:edge=dict(status='unavailable',reason='original Edge requires at least two trials within each grouped axis item')
                else:
                    engine=Edgewise_Contribution('',mice,list(frames),seed=ctx.section('Edge').get('seed',42),atlas_exclude=atlas.split('_'),frames=frames,trial_ids=trial_ids)
                    phi,dp=engine.run(repeat=ctx.section('Edge').get('repeat',50))
                    edge=dict(status='complete',phi=phi,dp=dp)
                results[atlas+'/'+label]=dict(status='complete',**result,edge=edge,
                    groups=halves.attrs['session_groups'],shared_axis_ids=[ids[i] for i in set(groups[0])&set(groups[1])],
                    warning='Historical odd-half shared samples are exploratory; not primary independent validation')
            except ValueError as exc:results[atlas+'/'+label]=dict(status='unavailable',reason=str(exc))
    save_pickle(folder/'statistics.pkl',results);save_json(folder/'results.json',results);save_json(folder/'cohort.json',cohort)
