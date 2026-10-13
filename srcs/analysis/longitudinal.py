"""Weekly representatives, explicit cross-week cohorts and animal trajectories."""
from itertools import combinations
import numpy as np
import pandas as pd
from srcs.analysis.identification import edge_vector, pattern_similarity, identify, edge_contribution
from srcs.utils.progress import progress_items


def weekly_representatives(ctx, records):
    weekly, rows = {}, []
    for atlas in ['_'.join(a) for a in ctx.atlases]:
        for week in progress_items(ctx.weeks, phase=f'{atlas}: weekly representatives'):
            for mouse in ctx.mice:
                runs = [r for r in ctx.runs if (atlas,week,mouse,r) in records]
                if runs:
                    frames = [records[atlas,week,mouse,r] for r in runs]
                    vectors = [edge_vector(f) for f in frames]
                    z = np.mean(vectors, axis=0)
                    matrix = np.eye(len(frames[0]))
                    iu = np.triu_indices(len(matrix),1)
                    matrix[iu]=np.tanh(z);matrix[(iu[1],iu[0])]=matrix[iu]
                    weekly[atlas,week,mouse] = dict(z_edges=z,
                        matrix=pd.DataFrame(matrix,index=frames[0].index,columns=frames[0].columns),
                        runs=runs, complete=len(runs)==len(ctx.runs))
                rows.append(dict(atlas=atlas,week=week,mouse=mouse,n_runs=len(runs),
                    runs=';'.join(runs),missing_runs=';'.join(r for r in ctx.runs if r not in runs),
                    status='complete' if len(runs)==len(ctx.runs) else 'single_or_partial' if runs else 'missing'))
    return weekly, pd.DataFrame(rows)


def select_cohort(ctx, weekly):
    atlases=['_'.join(a) for a in ctx.atlases]
    mode=ctx.section('ID').get('completeness_mode','strict')
    if mode not in ('strict','cover'): raise ValueError('completeness_mode must be strict or cover')
    spans=[ctx.weeks] if mode=='strict' else [ctx.weeks[i:j] for i in range(len(ctx.weeks)) for j in range(i+2,len(ctx.weeks)+1)]
    candidates=[]
    for weeks in spans:
        mice=[m for m in ctx.mice if all(weekly.get((a,w,m),{}).get('complete',False) for a in atlases for w in weeks)]
        candidates.append(dict(weeks=weeks,mice=mice,score=len(weeks)*len(mice)))
    eligible=[c for c in candidates if len(c['mice'])>=2 and len(c['weeks'])>=2]
    best=max(eligible,key=lambda c:(c['score'],len(c['mice']),-ctx.weeks.index(c['weeks'][0]))) if eligible else dict(mice=[],weeks=[],score=0)
    return dict(policy=mode,**best,excluded_mice=[m for m in ctx.mice if m not in best['mice']],
                excluded_weeks=[w for w in ctx.weeks if w not in best['weeks']],candidates=candidates,
                tie_break='area, mouse count, earliest configured week',
                status='complete' if eligible else 'insufficient_complete_cohort')


def comparison_plan(ctx, weeks):
    opts=ctx.section('longitudinal')
    modes=opts.get('comparisons',['observation_start','adjacent','all_pairs'])
    allowed={'observation_start','adjacent','all_pairs','windows'}
    if not isinstance(modes,list) or set(modes)-allowed: raise ValueError('Unknown longitudinal comparison')
    out=[]
    for mode in modes:
        if mode=='observation_start': pairs=[(weeks[0],w) for w in weeks[1:]] if weeks else []
        elif mode=='adjacent': pairs=list(zip(weeks,weeks[1:]))
        elif mode=='all_pairs': pairs=list(combinations(weeks,2))
        else:
            early,late=opts.get('early_weeks',[]),opts.get('late_weeks',[])
            if not early or not late or set(early)&set(late) or any(w not in ctx.weeks for w in early+late):
                raise ValueError('Early/late windows require nonempty, disjoint configured weeks')
            if any(w not in weeks for w in early+late):
                out.append(dict(name='windows',first=early,second=late,status='unavailable',reason='window outside retained cohort'))
                continue
            out.append(dict(name='windows',first=early,second=late));continue
        out.extend(dict(name=f'{mode}:{a}->{b}',first=[a],second=[b]) for a,b in pairs)
    return out


def analyze_between_weeks(ctx, weekly):
    cohort=select_cohort(ctx,weekly);results={};rows=[]
    for atlas in ['_'.join(a) for a in ctx.atlases]:
        for plan in progress_items(comparison_plan(ctx,cohort['weeks']), phase=f'{atlas}: comparisons'):
            key=atlas+'/'+plan['name']
            if plan.get('status')=='unavailable':results[key]=plan;continue
            mice=cohort['mice']
            vectors=[np.array([np.mean([weekly[atlas,w,m]['z_edges'] for w in group],axis=0) for m in mice]) for group in (plan['first'],plan['second'])]
            try:
                result=identify(*vectors,mice,n_iter=ctx.section('ID').get('n_iter',1000) if ctx.section('ID').get('shuffle',True) else 0,seed=ctx.session_seed)
            except ValueError as exc:
                results[key]=dict(status='unavailable',reason=str(exc),**plan);continue
            edge=edge_contribution(np.tanh(vectors[0]),np.tanh(vectors[1]))
            results[key]=dict(status='complete',**plan,**result,edge=edge,roi_labels=list(weekly[atlas,plan['first'][0],mice[0]]['matrix'].index))
            for i,m in enumerate(mice):
                rows.append(dict(atlas=atlas,comparison=plan['name'],mouse=m,n_mice=len(mice),
                    first_weeks=';'.join(plan['first']),second_weeks=';'.join(plan['second']),
                    self_z=result['self_z'][i],other_z=result['other_z'][i],advantage_z=result['advantage_z'][i],
                    self_match_percentile=result['individual_metrics']['self_match_percentile'][i],
                    relative_effect=result['relative_effect'][i],cohen_dz=result['effect']['cohen_dz'],
                    accuracy=result['accuracy'],p_permutation=result.get('p_permutation',np.nan),
                    first_run_count=sum(len(weekly[atlas,w,m]['runs']) for w in plan['first']),
                    second_run_count=sum(len(weekly[atlas,w,m]['runs']) for w in plan['second'])))
    return results,pd.DataFrame(rows),cohort


def trajectories(ctx, weekly):
    rows=[]
    metadata=ctx.section('study'); dates=metadata.get('acquisition_dates',{})
    for atlas in ['_'.join(a) for a in ctx.atlases]:
        for mouse in ctx.mice:
            observed=[w for w in ctx.weeks if (atlas,w,mouse) in weekly]
            if not observed:continue
            start=observed[0]
            for i,week in enumerate(observed):
                record=weekly[atlas,week,mouse];z=record['z_edges']
                previous=observed[i-1] if i else None
                row=dict(dataset=ctx.dataset,group=metadata.get('group','unspecified'),atlas=atlas,mouse=mouse,week=week,
                    week_index=ctx.weeks.index(week),date=dates.get(mouse,{}).get(week),
                    observation_start=start,n_runs=len(record['runs']),complete=record['complete'],
                    previous_week=previous,gap_configured_steps=ctx.weeks.index(week)-ctx.weeks.index(previous) if previous else None,
                    distance_from_start_z=float(np.sqrt(np.mean((z-weekly[atlas,start,mouse]['z_edges'])**2))),
                    adjacent_distance_z=float(np.sqrt(np.mean((z-weekly[atlas,previous,mouse]['z_edges'])**2))) if previous else np.nan,
                    self_pattern_from_start=np.nan,other_pattern_from_start=np.nan,identity_advantage=np.nan,
                    adjacent_pattern_r=np.nan)
                try:row['self_pattern_from_start']=pattern_similarity(weekly[atlas,start,mouse]['z_edges'],z)
                except ValueError:pass
                if previous:
                    try:row['adjacent_pattern_r']=pattern_similarity(weekly[atlas,previous,mouse]['z_edges'],z)
                    except ValueError:pass
                other=[]
                for m in ctx.mice:
                    if m!=mouse and (atlas,start,m) in weekly:
                        try:other.append(pattern_similarity(weekly[atlas,start,m]['z_edges'],z))
                        except ValueError:pass
                if other:
                    row['other_pattern_from_start']=float(np.mean(other));row['identity_advantage']=row['self_pattern_from_start']-row['other_pattern_from_start']
                row['n_reference_others']=len(other)
                rows.append(row)
    return pd.DataFrame(rows)


def control_deviation(ctx, weekly, controls):
    """Same-week reference only; leave out the same animal in control analyses."""
    rows=[]
    for (atlas,week,mouse),record in weekly.items():
        peers=[(m,c) for (a,w,m),c in controls.items() if a==atlas and w==week and m!=mouse and c.get('complete',False)]
        labels=list(record['matrix'].index)
        for _,c in peers:
            if list(c['matrix'].index)!=labels or list(c['matrix'].columns)!=labels:
                raise ValueError('Control ROI labels/order do not match target')
            if not np.isfinite(c['z_edges']).all() or np.asarray(c['z_edges']).shape!=np.asarray(record['z_edges']).shape:
                raise ValueError('Invalid control edge vector')
        row=dict(atlas=atlas,week=week,mouse=mouse,n_controls=len(peers),
                 status='available' if peers else 'no_complete_same_week_controls',
                 distance_from_control_z=np.nan,control_pattern_r=np.nan)
        if peers:
            reference=np.mean([c['z_edges'] for _,c in peers],axis=0)
            row['distance_from_control_z']=float(np.sqrt(np.mean((record['z_edges']-reference)**2)))
            try:row['control_pattern_r']=pattern_similarity(record['z_edges'],reference)
            except ValueError:pass
        rows.append(row)
    return pd.DataFrame(rows)
