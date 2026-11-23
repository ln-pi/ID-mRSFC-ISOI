"""Standalone saved-SeedFC supplement: variability, matching and paired-week ICC.

python plugins/seedfc/supplement.py --outputs A_OUTPUT N_OUTPUT
Optional YAML: `seedfc_supplement` mapping; see DEFAULTS for supported settings.
Writes a new folder exclusively; never edits inputs or joins the main stage cache.
"""
from __future__ import annotations

import argparse
import base64
from datetime import datetime, timezone
import hashlib
import html
from itertools import combinations_with_replacement
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml


DEFAULTS = dict(reference='03w', cohorts=['fixed_window', 'available_pair'],
    cv_mean_threshold=1e-6, numerical_tolerance=1e-12, icc_min_mice=3,
    dpi=150, formats=['png', 'pdf'], cmap='viridis', signed_cmap='coolwarm',
    missing_color='lightgray', series_colors=['#356caa', '#bd502f'],
    view_markers=['o', '^', 's'], panel_size=[5.2, 3.8], font_size=9,
    matrix_weeks=['03w', '07w', '13w'])


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def clean(value):
    if isinstance(value, dict): return {str(k): clean(v) for k,v in value.items()}
    if isinstance(value, (list, tuple)): return [clean(v) for v in value]
    if isinstance(value, np.ndarray): return clean(value.tolist())
    if isinstance(value, np.generic): return clean(value.item())
    if isinstance(value, float) and not np.isfinite(value): return None
    return value


def json_write(path, value):
    Path(path).write_text(json.dumps(clean(value), ensure_ascii=False, indent=2,
                                   allow_nan=False), encoding='utf-8')


def correlation(a, b):
    """Pattern Pearson r; undefined constant/nonfinite vectors remain NaN."""
    def normalize(x):
        finite = np.isfinite(x).all(1)
        varied = np.ptp(np.where(finite[:,None], x, 0), axis=1) > 0
        z = np.where(finite[:,None], x, 0)
        z = z-z.mean(1, keepdims=True)
        length = np.linalg.norm(z, axis=1, keepdims=True)
        return np.divide(z,length,out=np.zeros_like(z),where=length>0), finite & varied
    x, ax = normalize(a); y, by = normalize(b)
    result = np.clip(x@y.T, -1, 1)
    result[~(ax[:,None]&by[None,:])] = np.nan
    return result


def paired_icc(a, b, minimum=3, tolerance=1e-12):
    """Two-way mixed, single measurement: A absolute agreement, C consistency.

    Rows are mice, two columns are fixed acquisitions. No edge is a subject.
    ANOVA estimators may be negative; do not clip. Undefined denominator -> NaN.
    """
    n, e = a.shape
    empty = np.full(e, np.nan)
    if n < minimum: return empty.copy(), empty.copy(), ['insufficient_mice']*e
    y = np.stack([a,b],axis=1)
    grand = y.mean((0,1)); row = y.mean(1); col = y.mean(0)
    ms_mouse = 2*np.sum((row-grand)**2,axis=0)/(n-1)
    ms_visit = n*np.sum((col-grand)**2,axis=0)
    error = np.sum((y-row[:,None]-col[None]+grand)**2,axis=(0,1))/(n-1)
    numerator = ms_mouse-error
    da = ms_mouse+error+2*(ms_visit-error)/n
    dc = ms_mouse+error
    constant = np.ptp(y,axis=(0,1)) <= tolerance
    absolute = np.divide(numerator,da,out=empty.copy(),where=(da>tolerance)&~constant)
    consistency = np.divide(numerator,dc,out=empty.copy(),where=(dc>tolerance)&~constant)
    status = np.where(constant,'constant',np.where(np.isfinite(absolute)&np.isfinite(consistency),
                      'complete','undefined_denominator'))
    return absolute, consistency, status.tolist()


def matching(s):
    n = len(s)
    if n < 2 or not np.isfinite(s).all(): return None
    own = np.diag(s); mask = ~np.eye(n,dtype=bool)
    other = s[mask].reshape(n,n-1)
    reverse = s.T[mask].reshape(n,n-1)
    row = (own[:,None]>other).all(1); col = (own[:,None]>reverse).all(1)
    row_tie = (own[:,None]>=other).all(1) & (own[:,None]==other).any(1)
    col_tie = (own[:,None]>=reverse).all(1) & (own[:,None]==reverse).any(1)
    self_mean = own.mean(); other_mean = other.mean()
    return dict(fa_query=row.mean(),fa_reverse=col.mean(),fa_bidirectional=(row&col).mean(),
        self_similarity=self_mean,other_similarity=other_mean,
        differential_identifiability=100*(self_mean-other_mean),chance=1/n,
        query_ties=int(row_tie.sum()),reverse_ties=int(col_tie.sum()))


def load_output(root, atlas):
    root = Path(root).resolve(); path=root/'results/edges'/f'{atlas}.npz'
    marker=path.with_suffix('.json'); info=json.loads(marker.read_text(encoding='utf-8'))
    if sha(path)!=info['sha256']: raise ValueError(f'Changed edge store: {path}')
    with np.load(path,allow_pickle=False) as z:
        meta=json.loads(str(z['metadata']))
        a={k:z[k].copy() for k in ('F_z','F_r','valid_record','edge_index')}
    cfg_path=root/'metadata/analysis_config.json'
    cfg=json.loads(cfg_path.read_text(encoding='utf-8'))['config']
    identity=root/'metadata/data_identity.json'
    scope=json.loads(identity.read_text(encoding='utf-8'))['scope']
    if scope['source']['id']!='seedfc': raise ValueError('SeedFC outputs required')
    mapping=cfg['descriptive']['network_mapping']; labels=meta['roi_labels']
    if not set(labels)<=set(mapping): raise ValueError('Explicit ROI network mapping required')
    n,w,r,e=a['F_z'].shape
    if (n,w,r)!=(len(meta['mouse_ids']),len(meta['week_ids']),len(meta['run_ids'])):
        raise ValueError('Saved axes differ')
    if not np.array_equal(a['edge_index'],np.array(np.triu_indices(len(labels),1))):
        raise ValueError('Unexpected edge order')
    if not np.isfinite(a['F_z'][a['valid_record']]).all(): raise ValueError('Nonfinite valid FC')
    views=list(meta['run_ids'])+['mean_fc']
    values=np.concatenate([a['F_z'],a['F_z'].mean(2,keepdims=True)],axis=2)
    valid=np.concatenate([a['valid_record'],a['valid_record'].all(2,keepdims=True)],axis=2)
    # Fixed cohort is joint across all configured weeks and Runs, shared by every comparison.
    pattern_valid=valid & (np.ptp(values,axis=-1)>0) & np.isfinite(values).all(-1)
    fixed=pattern_valid.all((1,2))
    groups=[('whole','All',np.arange(e))]
    for roi in mapping:
        if roi in labels:
            j=labels.index(roi); groups.append(('roi',roi,np.flatnonzero((a['edge_index']==j).any(0))))
    networks=list(dict.fromkeys(mapping[x] for x in mapping if x in labels))
    for first,second in combinations_with_replacement(networks,2):
        ids=[k for k,(i,j) in enumerate(a['edge_index'].T)
             if {mapping[labels[i]],mapping[labels[j]]}=={first,second}]
        if ids: groups.append(('network_pair',f'{first} | {second}',np.array(ids)))
    return dict(root=root,meta=meta,values=values,valid=valid,pattern_valid=pattern_valid,
        fixed=fixed,views=views,groups=groups,edges=a['edge_index'],mapping=mapping,
        series=scope['series'],inputs={str(p):sha(p) for p in (path,marker,cfg_path,identity)})


def compute(source, cfg, tables, arrays, audit):
    meta=source['meta']; weeks=meta['week_ids']; mice=np.array(meta['mouse_ids'])
    ref=weeks.index(cfg['reference']); views=source['views']; x=source['values']
    base_series=source['series']; labels=meta['roi_labels']; edges=source['edges']
    comparisons=[(q,r) for q in range(len(views)-1) for r in range(len(views)-1)]+[(len(views)-1,)*2]
    for cohort in cfg['cohorts']:
        for w,week in enumerate(weeks):
            for v,view in enumerate(views):
                keep=source['valid'][:,w,v].copy()
                if cohort=='fixed_window': keep &= source['fixed']
                ids=np.flatnonzero(keep); n=len(ids)
                base=dict(series=base_series,cohort=cohort,week=week,view=view,n_mice=n)
                audit.append(dict(**base,kind='weekly',mice=mice[ids].tolist()))
                if n<2: continue
                data=x[ids,w,v]; mean=data.mean(0); sd=data.std(0,ddof=1)
                cv=np.divide(sd,mean,out=np.full_like(sd,np.nan),where=mean!=0)
                stable=np.abs(mean)>cfg['cv_mean_threshold']
                for e,(i,j) in enumerate(edges.T):
                    tables['weekly_edges'].append(dict(**base,edge=e,roi_1=labels[i],roi_2=labels[j],
                        mean=mean[e],sd=sd[e],cv=cv[e],cv_near_zero=not stable[e],cv_negative_mean=mean[e]<0))
                for level,region,es in source['groups']:
                    selected=cv[es][stable[es]&np.isfinite(cv[es])]
                    tables['weekly_regions'].append(dict(**base,level=level,region=region,n_edges=len(es),
                        mean_sd=sd[es].mean(),median_sd=np.median(sd[es]),
                        median_cv=np.median(selected) if len(selected) else np.nan,
                        cv_usable_edges=len(selected),cv_negative_mean_edges=int((mean[es]<0).sum())))
                # Constant patterns cannot enter pattern similarity; audited separately.
                pattern_ids=ids[source['pattern_valid'][ids,w,v]]
                sim=correlation(x[pattern_ids,w,v],x[pattern_ids,w,v])
                key=f'{base_series}_{cohort}_{week}_{view}_within'
                arrays[key]=sim; arrays[key+'_mice']=mice[pattern_ids]
                m=len(sim); avg=sim[np.triu_indices(m,1)].mean() if m>=2 else np.nan
                tables['within_summary'].append(dict(**base,n_pattern_mice=m,mean_similarity=avg,
                    excluded_constant_mice='|'.join(mice[ids[~source['pattern_valid'][ids,w,v]]])))
                for i in range(m):
                    for j in range(i+1,m):
                        tables['within_pairs'].append(dict(**base,mouse_i=mice[pattern_ids[i]],
                            mouse_j=mice[pattern_ids[j]],similarity=sim[i,j]))
            for q,r in comparisons:
                comparison=f'{views[q]}__to__{views[r]}'
                base=dict(series=base_series,cohort=cohort,week=week,query=views[q],
                          reference=views[r],comparison=comparison,reference_week=cfg['reference'])
                # Identical acquisition reference is never treated as evidence.
                if w==ref and q==r:
                    audit.append(dict(**base,kind='reference',status='same_acquisition_excluded',mice=[]))
                    continue
                keep=source['pattern_valid'][:,w,q]&source['pattern_valid'][:,ref,r]
                if cohort=='fixed_window': keep &= source['fixed']
                ids=np.flatnonzero(keep); n=len(ids); base['n_mice']=n
                audit.append(dict(**base,kind='reference',status='complete' if n>=2 else 'insufficient_mice',mice=mice[ids].tolist()))
                if n<2: continue
                a=x[ids,w,q]; b=x[ids,ref,r]; s=correlation(a,b)
                key=f'{base_series}_{cohort}_{week}_{comparison}'
                arrays[key]=s; arrays[key+'_mice']=mice[ids]
                result=matching(s); tables['fingerprint'].append(dict(**base,**result))
                for i in range(n):
                    foreign=np.delete(s[i],i); rev=np.delete(s[:,i],i)
                    tables['self_similarity'].append(dict(**base,mouse=mice[ids[i]],self_similarity=s[i,i],
                        other_similarity=(foreign.mean()+rev.mean())/2,
                        advantage=s[i,i]-(foreign.mean()+rev.mean())/2))
                aa,cc,status=paired_icc(a,b,cfg['icc_min_mice'],cfg['numerical_tolerance'])
                for e,(i,j) in enumerate(edges.T):
                    tables['icc_edges'].append(dict(**base,edge=e,roi_1=labels[i],roi_2=labels[j],
                        icc_a1=aa[e],icc_c1=cc[e],status=status[e]))
                for level,region,es in source['groups']:
                    row=dict(**base,level=level,region=region,n_edges=len(es))
                    for metric,values in (('icc_a1',aa[es]),('icc_c1',cc[es])):
                        finite=values[np.isfinite(values)]
                        row[metric]=np.median(finite) if len(finite) else np.nan
                        row[metric+'_valid_edges']=len(finite)
                    tables['icc_regions'].append(row)


def render(sources, cfg, tables, arrays, out):
    plt.rcParams.update({'font.size':cfg['font_size']})
    weeks=sources[0]['meta']['week_ids']; series=[s['series'] for s in sources]
    views=sources[0]['views']; files=[]; plots=out/'figures'; plots.mkdir()
    def save(fig,name):
        for ext in cfg['formats']:
            path=plots/f'{name}.{ext}'; fig.savefig(path,dpi=cfg['dpi'],bbox_inches='tight')
            if ext=='png': files.append(path)
        plt.close(fig)
    def heat(frame,metric,name,cohort,comparison=None):
        for level in ('network_pair','roi'):
            subset=frame[(frame.cohort==cohort)&(frame.level==level)]
            if comparison: subset=subset[subset.comparison==comparison]
            panels=list(subset.groupby(['series'] if comparison else ['series','view'],sort=False))
            if not panels: continue
            signed=metric!='mean_sd'; values=subset[metric].to_numpy(); finite=values[np.isfinite(values)]
            limit=max(float(np.max(np.abs(finite))) if len(finite) else 1.,cfg['numerical_tolerance'])
            if metric.startswith('icc'): limit=1.
            fig,axes=plt.subplots(len(panels),1,squeeze=False,
                figsize=(cfg['panel_size'][0]*2,cfg['panel_size'][1]*len(panels)))
            for ax,(key,part) in zip(axes.flat,panels):
                regions=list(dict.fromkeys(part.region))
                mat=part.pivot(index='region',columns='week',values=metric).reindex(index=regions,columns=weeks)
                cm=plt.get_cmap(cfg['signed_cmap'] if signed else cfg['cmap']).with_extremes(bad=cfg['missing_color'])
                im=ax.imshow(mat,aspect='auto',cmap=cm,vmin=-limit if signed else 0,vmax=limit)
                ax.set_xticks(range(len(weeks)),weeks,rotation=90); ax.set_yticks(range(len(regions)),regions)
                ax.set_title(f'{key} | {cohort} | {metric}'+(' | exploratory signed CV' if metric=='median_cv' else ''))
                fig.colorbar(im,ax=ax,extend='both' if metric.startswith('icc') else 'neither')
            fig.tight_layout();save(fig,f'{name}_{level}_{cohort}')
    for cohort in cfg['cohorts']:
        heat(tables['weekly_regions'],'mean_sd','SD',cohort)
        heat(tables['weekly_regions'],'median_cv','CV_diagnostic',cohort)
        fig,axes=plt.subplots(1,len(views),squeeze=False,figsize=(cfg['panel_size'][0]*len(views),cfg['panel_size'][1]))
        for ax,view in zip(axes.flat,views):
            for k,label in enumerate(series):
                d=tables['within_summary']; d=d[(d.series==label)&(d.view==view)&(d.cohort==cohort)].set_index('week').reindex(weeks)
                ax.plot(weeks,d.mean_similarity,marker=cfg['view_markers'][k%len(cfg['view_markers'])],color=cfg['series_colors'][k%len(cfg['series_colors'])],label=label)
            ax.set_title(view);ax.tick_params(axis='x',rotation=90);ax.set_ylabel('Mean nonself pattern Pearson r');ax.legend()
        fig.suptitle(cohort);fig.tight_layout();save(fig,'within_similarity_'+cohort)
        for comparison in tables['fingerprint'].comparison.unique():
            fig,axes=plt.subplots(2,2,figsize=(cfg['panel_size'][0]*2,cfg['panel_size'][1]*2))
            for k,label in enumerate(series):
                d=tables['fingerprint'];d=d[(d.series==label)&(d.cohort==cohort)&(d.comparison==comparison)].set_index('week').reindex(weeks)
                style=dict(color=cfg['series_colors'][k%len(cfg['series_colors'])],marker=cfg['view_markers'][k%len(cfg['view_markers'])],label=label)
                axes[0,0].plot(weeks,d.fa_query,**style)
                axes[0,0].plot(weeks,d.fa_bidirectional,color=style['color'],ls='--',label=label+' strict bidirectional')
                axes[0,0].plot(weeks,d.chance,color=style['color'],ls=':',label=label+' 1/n query chance')
                axes[0,1].plot(weeks,d.self_similarity,**style)
                axes[1,0].plot(weeks,d.other_similarity,**style)
                axes[1,1].plot(weeks,d.differential_identifiability,**style)
            for ax,title in zip(axes.flat,('Fingerprint accuracy','Self similarity (r)','Others similarity (r)','Differential identifiability (100 x r difference)')):
                ax.set_title(title);ax.tick_params(axis='x',rotation=90);ax.legend(fontsize=cfg['font_size']-1)
            axes[0,0].set_ylim(0,1);fig.suptitle(f'{cohort} | {comparison} | reference {cfg["reference"]}')
            fig.tight_layout();save(fig,'fingerprint_'+cohort+'_'+comparison)
            fig,axes=plt.subplots(1,len(series),squeeze=False,
                figsize=(cfg['panel_size'][0]*len(series),cfg['panel_size'][1]))
            for ax,label in zip(axes.flat,series):
                d=tables['self_similarity']
                d=d[(d.series==label)&(d.cohort==cohort)&(d.comparison==comparison)]
                for mouse,part in d.groupby('mouse',sort=False):
                    part=part.set_index('week').reindex(weeks)
                    ax.plot(weeks,part.self_similarity,marker=cfg['view_markers'][0],label=mouse)
                ax.set_title(label);ax.set_ylabel('Individual self similarity (r)')
                ax.tick_params(axis='x',rotation=90);ax.legend(fontsize=cfg['font_size']-1)
            fig.suptitle(f'{cohort} | {comparison} | reference {cfg["reference"]}')
            fig.tight_layout();save(fig,'individual_self_'+cohort+'_'+comparison)
            for metric in ('icc_a1','icc_c1'):
                heat(tables['icc_regions'],metric,metric+'_'+comparison,cohort,comparison)
        # Representative original matrices; every matrix remains in NPZ.
        for mode in ('within', 'reference'):
            chosen=[]
            for source in sources:
                for week in cfg['matrix_weeks']:
                    suffix=views[0]+'_within' if mode=='within' else views[0]+'__to__'+views[0]
                    key=f'{source["series"]}_{cohort}_{week}_{suffix}'
                    if key in arrays: chosen.append((key,source['series'],week))
            if chosen:
                fig,axes=plt.subplots(1,len(chosen),squeeze=False,figsize=(cfg['panel_size'][0]*len(chosen),cfg['panel_size'][1]))
                for ax,(key,label,week) in zip(axes.flat,chosen):
                    im=ax.imshow(arrays[key],vmin=-1,vmax=1,cmap=cfg['signed_cmap'])
                    mice=arrays[key+'_mice'];ax.set_xticks(range(len(mice)),mice,rotation=90);ax.set_yticks(range(len(mice)),mice)
                    ax.set_title(f'{label} {week} | {mode} | {views[0]}');fig.colorbar(im,ax=ax)
                fig.tight_layout();save(fig,'matrices_'+mode+'_'+cohort)
    return files


def run(args):
    cfg=dict(DEFAULTS)
    if args.config:
        loaded=yaml.safe_load(Path(args.config).read_text(encoding='utf-8'))
        override=loaded.get('seedfc_supplement',{})
        unknown=set(override)-set(cfg)
        if unknown: raise ValueError(f'Unknown supplemental settings: {unknown}')
        cfg.update(override)
    if args.reference: cfg['reference']=args.reference
    if not cfg['cohorts'] or not set(cfg['cohorts'])<={'fixed_window','available_pair'}: raise ValueError('Unknown cohort')
    if cfg['cv_mean_threshold']<0 or cfg['numerical_tolerance']<=0 or cfg['icc_min_mice']<3: raise ValueError('Invalid numerical settings')
    if 'png' not in cfg['formats'] or not set(cfg['formats'])<={'png','pdf'}: raise ValueError('Formats must include png and optionally pdf')
    sources=[load_output(p,args.atlas) for p in args.outputs]
    if len({s['series'] for s in sources})!=len(sources): raise ValueError('Distinct series labels required')
    for s in sources:
        if s['meta']['roi_labels']!=sources[0]['meta']['roi_labels'] or s['meta']['week_ids']!=sources[0]['meta']['week_ids'] or s['views']!=sources[0]['views'] or s['mapping']!=sources[0]['mapping']:
            raise ValueError('Incompatible ordered ROIs, weeks, Runs or networks')
        if cfg['reference'] not in s['meta']['week_ids']: raise ValueError('Reference week absent')
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out=Path(args.destination or Path(__file__).resolve().parents[2]/'data/output/_supplement')/('SeedFC_AN_'+cfg['reference']+'_'+stamp)
    out.mkdir(parents=True,exist_ok=False)
    tables={name:[] for name in ('weekly_edges','weekly_regions','within_summary','within_pairs','fingerprint','self_similarity','icc_edges','icc_regions')}
    arrays={};audit=[]
    for s in sources:
        print('Computing',s['series'],'fixed mice',np.array(s['meta']['mouse_ids'])[s['fixed']].tolist(),flush=True)
        compute(s,cfg,tables,arrays,audit)
    numeric=out/'numeric';numeric.mkdir()
    tables={k:pd.DataFrame(v) for k,v in tables.items()}
    for name,table in tables.items(): table.to_csv(numeric/(name+'.csv'),index=False)
    np.savez_compressed(numeric/'matrices.npz',**arrays)
    json_write(numeric/'cohorts.json',audit)
    json_write(out/'effective_config.json',cfg)
    print('Rendering saved metrics',flush=True)
    files=render(sources,cfg,tables,arrays,out)
    overview=tables['fingerprint']; overview=overview[(overview.cohort=='fixed_window')&(overview['query']==overview.reference)]
    notes='''<p>输入为已有SeedFC edge store；14 ROI/91边（以实际元数据为准），Fisher-z空间。固定队列为所有配置周/Run/view共同完整且模式非恒定的动物；available_pair在周内指标处表示当周有效，在参考比较处表示两个端点共同有效。</p>
<p>SD为逐边跨鼠样本SD(ddof=1)，区域汇总为平均逐边SD。CV=SD/有符号均值，Fisher-z不是常规CV适用的比率尺度；仅供诊断。均值为0时CV未定义，近零阈值只用于标记及区域中位数筛选，不修改逐边CV；负分母保留。CV图不作为变异程度的主结论。</p>
<p>组内相似性为异鼠模式Pearson r的算术平均，不含对角线。Self/Others均基于模式r，I_diff=100(Self−Others)。FA主结果为目标周查03w的唯一最高匹配，另保存反向和严格双向；并列不算成功。1/n只对应单向查询机会水平，不对应严格双向统计量。</p>
<p>ICC逐边以鼠为对象，两次指定采集为固定测量，采用双向混合单次测量ICC(A,1)/ICC(C,1)的ANOVA公式；负值保留。区域图为逐边ICC中位数，不是对区域均值重新估计ICC；色标±1外保留扩展色条。03w同Run和同mean_fc是同记录，全部排除；03w跨Run保留。结果为描述性首轮计算，无新增CI或组间显著性检验。</p>
<p>方法：<a href="https://www.nature.com/articles/nn.4135">Finn et al. 2015</a>（识别框架，本项目严格双向口径另列）；<a href="https://www.nature.com/articles/s41598-018-25089-1">Amico &amp; Goñi 2018</a>（I_diff）；McGraw &amp; Wong 1996, Psychological Methods 1:30–46, DOI <a href="https://doi.org/10.1037/1082-989X.1.1.30">10.1037/1082-989X.1.1.30</a>（ICC A/C）；<a href="https://itl.nist.gov/div898/software/dataplot/refman2/auxillar/coefvari.htm">NIST CV适用边界</a>。</p>'''
    notes += '''<details><summary>计算公式与保存范围</summary><pre>SD_e = sqrt(sum_i (x_i,e - mean_i x_i,e)^2 / (n-1))
CV_e = SD_e / mean_i x_i,e (signed, exploratory)
S_ij(w) = Corr_edges(x_i,w,query, x_j,reference_week,reference)
Self = mean_i S_ii; Others = mean_i!=j S_ij; I_diff = 100(Self-Others)
ICC(A,1) = (MS_mouse-MS_error)/(MS_mouse+(k-1)MS_error+k(MS_visit-MS_error)/n)
ICC(C,1) = (MS_mouse-MS_error)/(MS_mouse+(k-1)MS_error), k=2
MS_mouse: between-mouse ANOVA mean square; MS_visit: between-visit mean square;
MS_error: residual mouse-by-visit mean square. All are calculated per FC edge.
All weekly within-group and reference matching matrices, with ordered animal IDs,
are retained in numeric/matrices.npz. Only selected weeks/Run01 are previewed.</pre></details>'''
    body=['<!doctype html><meta charset="utf-8"><title>SeedFC supplementary metrics</title><style>body{font:16px/1.8 system-ui;max-width:1350px;margin:auto;padding:25px}img{max-width:100%}table{border-collapse:collapse;font-size:12px;display:block;overflow:auto}td,th{border:1px solid #ccc;padding:5px}figure{margin:30px 0}figure:target{outline:3px solid teal}a,pre{overflow-wrap:anywhere}pre{white-space:pre-wrap}</style><h1>SeedFC A/N：组内变异、矩阵相似性与03w参考指纹</h1>',notes,
          '<h2>固定队列同Run概览</h2>',overview.to_html(index=False,float_format=lambda x:f'{x:.4f}'),'<h2>图表目录</h2><ul>']
    body.extend(f'<li><a href="#fig{i}">{html.escape(p.stem)}</a></li>' for i,p in enumerate(files))
    body.append('</ul>')
    for i,p in enumerate(files):
        body.append(f'<figure id="fig{i}"><figcaption>{html.escape(p.stem)}</figcaption><a href="{p.as_uri()}"><img loading="lazy" src="data:image/png;base64,{base64.b64encode(p.read_bytes()).decode()}"></a></figure>')
    body.append('<h2>数值下载</h2><ul>'+''.join(f'<li><a href="{p.as_uri()}">{p.name}</a></li>' for p in numeric.iterdir())+'</ul>')
    (out/'report.html').write_text(''.join(body),encoding='utf-8')
    inputs={p:h for s in sources for p,h in s['inputs'].items()}
    if args.config: inputs[str(Path(args.config).resolve())]=sha(args.config)
    for p,h in inputs.items():
        if sha(p)!=h: raise RuntimeError('Input changed during computation')
    json_write(out/'manifest.json',dict(status='complete',config=cfg,inputs=inputs,code={str(Path(__file__).resolve()):sha(__file__)},
        outputs={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()},
        cohorts={s['series']:np.array(s['meta']['mouse_ids'])[s['fixed']].tolist() for s in sources},
        uncertainty='descriptive only; no confidence intervals or group p-values'))
    print('RESULT',out,flush=True)
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outputs',nargs='+',required=True,type=Path)
    parser.add_argument('--atlas',default='FourNetworks')
    parser.add_argument('--reference',help='Reference label; default from effective settings')
    parser.add_argument('--config',type=Path,help='YAML containing seedfc_supplement settings')
    parser.add_argument('--destination',type=Path,help='Parent directory for a new timestamped result')
    run(parser.parse_args())
