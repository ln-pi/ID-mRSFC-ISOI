"""Standalone saved-SeedFC supplement: variability, matching and paired-week ICC.

python -m plugins.function.seedfc_supplement --outputs A_OUTPUT N_OUTPUT
Optional YAML overrides analysis and visualization settings.
Verified numerical results and reports are cached separately.
"""
from __future__ import annotations

import argparse
import base64
import html
import json
from pathlib import Path
from tempfile import mkdtemp
from types import SimpleNamespace
import inspect
import sys

# Direct script execution must find the same project modules as the batch runner.
if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


from srcs.analysis.magnitude import select_views
from srcs.analysis.topology import edge_layout
from srcs.core.artifacts import sha256 as sha, write_json
from srcs.core.identity import digest
from srcs.core.runtime import exclusive_run
from srcs.utils.config import PROJECT_ROOT, load_config
from srcs.utils.visualization import displayed_weeks, merge, saved_config, settings, validate
from srcs.visualization.scales import resolve


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


def load_output(root, atlas, config=None):
    root = Path(root).resolve(); path=root/'results/edges'/f'{atlas}.npz'
    marker=path.with_suffix('.json'); info=json.loads(marker.read_text(encoding='utf-8'))
    if sha(path)!=info['sha256']: raise ValueError(f'Changed edge store: {path}')
    with np.load(path,allow_pickle=False) as z:
        meta=json.loads(str(z['metadata']))
        a={k:z[k].copy() for k in ('F_z','F_r','valid_record','edge_index')}
    cfg = saved_config(root) if config is None else config
    identity=root/'metadata/data_identity.json'
    scope=json.loads(identity.read_text(encoding='utf-8'))['scope']
    if scope['source']['id']!='seedfc': raise ValueError('SeedFC outputs required')
    descriptive = cfg.get("descriptive", {})
    labels = meta["roi_labels"]
    n,w,r,e=a['F_z'].shape
    if (n,w,r)!=(len(meta['mouse_ids']),len(meta['week_ids']),len(meta['run_ids'])):
        raise ValueError('Saved axes differ')
    if not np.array_equal(a['edge_index'],np.array(np.triu_indices(len(labels),1))):
        raise ValueError('Unexpected edge order')
    if not np.isfinite(a['F_z'][a['valid_record']]).all(): raise ValueError('Nonfinite valid FC')
    store = SimpleNamespace(weeks=meta["week_ids"], runs=meta["run_ids"],
                            F_z=a["F_z"], valid_record=a["valid_record"])
    weeks, views, values, valid, _ = select_views(store, descriptive)
    meta["week_ids"] = weeks
    # Fixed cohort is joint across all configured weeks and Runs, shared by every comparison.
    pattern_valid=valid & (np.ptp(values,axis=-1)>0) & np.isfinite(values).all(-1)
    fixed=pattern_valid.all((1,2))
    layout, _, topology = edge_layout(labels, a["edge_index"], descriptive.get("network_mapping"))
    groups = [("whole", "All", np.arange(e))]
    groups.extend((group["level"], group["name"], group["indices"]) for group in layout
                  if group["level"] in ("roi", "network_pair"))
    return dict(root=root,meta=meta,values=values,valid=valid,pattern_valid=pattern_valid,
        fixed=fixed,views=views,groups=groups,edges=a['edge_index'],mapping=topology["roi_network_mapping"],
        series=scope['series'], inputs={str(p):sha(p) for p in (path, marker, identity)})


def compute(source, cfg, tables, arrays, audit):
    meta=source['meta']; weeks=meta['week_ids']; mice=np.array(meta['mouse_ids'])
    ref=weeks.index(cfg['reference']); views=source['views']; x=source['values']
    base_series=source['series']; labels=meta['roi_labels']; edges=source['edges']
    run_indices = [i for i, view in enumerate(views) if view in meta["run_ids"]]
    comparisons = [(q, r) for q in run_indices for r in run_indices]
    if "mean_fc" in views:
        comparisons.append((views.index("mean_fc"),) * 2)
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


def render(sources, cfg, tables, arrays, out, viz):
    weeks = displayed_weeks(sources[0]['meta']['week_ids'], viz, cfg['reference'])
    series = [s['series'] for s in sources]
    scales = []
    tables = {k: v[v.week.isin(weeks)] if 'week' in v else v for k, v in tables.items()}
    cfg = {**cfg, **viz['supplement'], 'series_colors': viz['series_colors'], 'view_markers': viz['view_markers']}
    plt.rcParams.update({'font.size': cfg['font_size']})
    views=sources[0]['views']; files=[]; plots=out/'figures'; plots.mkdir()
    def save(fig,name):
        for ext in cfg['formats']:
            path=plots/f'{name}.{ext}'; fig.savefig(path,dpi=cfg['dpi'],bbox_inches='tight')
            if ext=='png': files.append(path)
        plt.close(fig)
    def heat(frame,metric,name,cohort,comparison=None):
        if frame.empty:
            return
        for level in ('network_pair','roi'):
            subset=frame[(frame.cohort==cohort)&(frame.level==level)]
            if comparison: subset=subset[subset.comparison==comparison]
            panels=list(subset.groupby(['series'] if comparison else ['series','view'],sort=False))
            if not panels: continue
            policy = {'mean_sd': 'sd', 'median_cv': 'cv'}.get(metric, 'icc')
            scale = resolve(subset[metric], viz, policy, anchor_zero=policy in ('sd', 'cv'))
            if scale is None:
                continue
            scales.append(scale.record(metric=metric, level=level, cohort=cohort, comparison=comparison,
                                       weeks=weeks, series=series, views=views))
            fig,axes=plt.subplots(len(panels),1,squeeze=False,
                figsize=(cfg['panel_size'][0]*2,cfg['panel_size'][1]*len(panels)))
            for ax,(key,part) in zip(axes.flat,panels):
                regions=list(dict.fromkeys(part.region))
                mat=part.pivot(index='region',columns='week',values=metric).reindex(index=regions,columns=weeks)
                im = ax.imshow(mat, aspect='auto', **scale.image)
                ax.set_xticks(range(len(weeks)),weeks,rotation=90); ax.set_yticks(range(len(regions)),regions)
                ax.set_title(f'{key} | {cohort} | {metric}'+(' | exploratory signed CV' if metric=='median_cv' else ''))
                fig.colorbar(im,ax=ax,extend=scale.extend)
            fig.tight_layout();save(fig,f'{name}_{level}_{cohort}')
    for cohort in cfg['cohorts']:
        heat(tables['weekly_regions'],'mean_sd','SD',cohort)
        heat(tables['weekly_regions'],'median_cv','CV_diagnostic',cohort)
        if tables['within_summary'].empty:
            continue
        fig,axes=plt.subplots(1,len(views),squeeze=False,figsize=(cfg['panel_size'][0]*len(views),cfg['panel_size'][1]))
        for ax,view in zip(axes.flat,views):
            for k,label in enumerate(series):
                d=tables['within_summary']; d=d[(d.series==label)&(d.view==view)&(d.cohort==cohort)].set_index('week').reindex(weeks)
                ax.plot(weeks,d.mean_similarity,marker=cfg['view_markers'][k%len(cfg['view_markers'])],color=cfg['series_colors'][k%len(cfg['series_colors'])],label=label)
            ax.set_title(view);ax.tick_params(axis='x',rotation=90);ax.set_ylabel('Mean nonself pattern Pearson r');ax.legend()
        fig.suptitle(cohort);fig.tight_layout();save(fig,'within_similarity_'+cohort)
        for comparison in tables['fingerprint'].get('comparison', pd.Series(dtype=str)).unique():
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
                    if week not in weeks:
                        continue
                    suffix=views[0]+'_within' if mode=='within' else views[0]+'__to__'+views[0]
                    key=f'{source["series"]}_{cohort}_{week}_{suffix}'
                    if key in arrays and (mode != 'within' or len(arrays[key]) > 1):
                        chosen.append((key, source['series'], week))
            if chosen:
                metric = mode + "_similarity"
                population = np.concatenate([arrays[key].ravel() for key, _, _ in chosen])
                scale = resolve(population, viz, metric, bounds=(-1, 1))
                if scale is None:
                    continue
                scales.append(scale.record(metric=metric, cohort=cohort, weeks=list(dict.fromkeys(w for _, _, w in chosen)),
                                           series=series, view=views[0]))
                fig,axes=plt.subplots(1,len(chosen),squeeze=False,figsize=(cfg['panel_size'][0]*len(chosen),cfg['panel_size'][1]))
                for ax,(key,label,week) in zip(axes.flat,chosen):
                    im = ax.imshow(arrays[key], **scale.image)
                    mice=arrays[key+'_mice'];ax.set_xticks(range(len(mice)),mice,rotation=90);ax.set_yticks(range(len(mice)),mice)
                    ax.set_title(f'{label} {week} | {mode} | {views[0]}');fig.colorbar(im,ax=ax,extend=scale.extend)
                fig.tight_layout();save(fig,'matrices_'+mode+'_'+cohort)
    json_write(out / 'color_scales.json', scales)
    return files


TABLE_NAMES = ("weekly_edges", "weekly_regions", "within_summary", "within_pairs", "fingerprint",
               "self_similarity", "icc_edges", "icc_regions")
NUMERIC_OPTIONS = ("reference", "cohorts", "cv_mean_threshold", "numerical_tolerance", "icc_min_mice")


def verified_cache(folder):
    manifest = folder / "manifest.json"
    try:
        record = json.loads(manifest.read_text(encoding="utf-8"))
        if record["status"] == "complete" and record["outputs"] and all(
                (folder / path).is_file() and sha(folder / path) == checksum
                for path, checksum in record["outputs"].items()):
            return record
    except (OSError, ValueError, KeyError, TypeError):
        pass
    return None


def publish_cache(pending, target, **metadata):
    write_json(pending / "manifest.json", dict(status="complete", **metadata,
        outputs={str(file.relative_to(pending)): sha(file) for file in pending.rglob("*") if file.is_file()}))
    # A corrupt previous generation is retained separately for inspection.
    if target.exists():
        from uuid import uuid4
        target.rename(target.with_name(f".{target.name}.{uuid4().hex}.invalid"))
    pending.rename(target)


def calculate(sources, cfg, destination, *, reuse, require_existing=False):
    # Rendering edits must not invalidate numerical results.
    numerical_code = digest([inspect.getsource(func) for func in
        (load_output, select_views, edge_layout, compute, correlation, matching, paired_icc)])
    inputs = {path: checksum for source in sources for path, checksum in source["inputs"].items()}
    numeric_options = {key: cfg[key] for key in NUMERIC_OPTIONS}
    selection = [dict(series=s["series"], weeks=s["meta"]["week_ids"], views=s["views"],
                      mapping=s["mapping"], regions=[(level, name) for level, name, _ in s["groups"]])
                 for s in sources]
    key = digest([inputs, numeric_options, selection, numerical_code])
    target = destination / "numeric" / key
    with exclusive_run(target.parent / f"{key}.lock"):
        cached = verified_cache(target) if reuse else None
        if cached:
            return target, cached, True
        if require_existing:
            raise ValueError("No verified numerical results for these inputs and analysis settings")
        tables = {name: [] for name in TABLE_NAMES}
        arrays, audit = {}, []
        for source in sources:
            compute(source, cfg, tables, arrays, audit)
        pending = Path(mkdtemp(prefix=f".{key}.", dir=target.parent))
        for name, rows in tables.items():
            pd.DataFrame(rows).to_csv(pending / f"{name}.csv", index=False)
        np.savez_compressed(pending / "matrices.npz", **arrays)
        json_write(pending / "cohorts.json", audit)
        if any(sha(path) != checksum for path, checksum in inputs.items()):
            raise RuntimeError("SeedFC inputs changed during supplementary calculation")
        publish_cache(pending, target, inputs=inputs, options=numeric_options,
                      selection=selection, numerical_code=numerical_code)
        return target, verified_cache(target), False


def read_table(path):
    try:
        return pd.read_csv(path)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def render_report(sources, cfg, viz, numeric, record, destination):
    from srcs.visualization import scales
    render_code = digest([inspect.getsource(func) for func in
                          (render, render_report, write_report, displayed_weeks, resolve, scales.step, scales.palette)])
    presentation = dict(visualization=viz, matrix_weeks=cfg["matrix_weeks"])
    key = digest([record, presentation, render_code])
    out = destination / "reports" / key
    with exclusive_run(out.parent / f"{key}.lock"):
        if verified_cache(out):
            return out, True
        # Render under the lock. Relative links stay valid after publication.
        pending = Path(mkdtemp(prefix=f".{key}.", dir=out.parent))
        tables = {name: read_table(numeric / f"{name}.csv") for name in TABLE_NAMES}
        with np.load(numeric / "matrices.npz", allow_pickle=False) as arrays, plt.rc_context():
            files = render(sources, cfg, tables, arrays, pending, viz)
        write_report(pending, files, tables, numeric, cfg, viz)
        json_write(pending / "effective_config.json", {**cfg, **presentation})
        publish_cache(pending, out, numeric=str(numeric), inputs=record["inputs"],
                      presentation=presentation, render_code=render_code)
    return out, False


def write_report(out, files, tables, numeric, cfg, viz):
    overview = tables["fingerprint"]
    if not overview.empty:
        weeks = displayed_weeks(cfg["weeks"], viz, cfg["reference"])
        overview = overview[(overview["query"] == overview.reference) & overview.week.isin(weeks)]
    notes='''<p>输入为已有SeedFC edge store；FC模式使用保存的ROI与边顺序，在Fisher-z空间计算。固定队列为所有配置周/Run/view共同完整且模式非恒定的动物；available_pair在周内指标处表示当周有效，在参考比较处表示两个端点共同有效。</p>
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
are retained in numeric/matrices.npz. The report previews the configured matrix weeks and the selected Run view.</pre></details>'''
    notes = notes.replace('03w', html.escape(cfg['reference']))
    body=['<!doctype html><meta charset="utf-8"><title>SeedFC supplementary metrics</title><style>body{font:16px/1.8 system-ui;max-width:1350px;margin:auto;padding:25px}img{max-width:100%}table{border-collapse:collapse;font-size:12px;display:block;overflow:auto}td,th{border:1px solid #ccc;padding:5px}figure{margin:30px 0}figure:target{outline:3px solid teal}a,pre{overflow-wrap:anywhere}pre{white-space:pre-wrap}</style><h1>SeedFC：组内变异、矩阵相似性与参考周指纹</h1>',notes,
          '<h2>同Run概览</h2>',overview.to_html(index=False,float_format=lambda x:f'{x:.4f}'),'<h2>图表目录</h2><ul>']
    body.extend(f'<li><a href="#fig{i}">{html.escape(p.stem)}</a></li>' for i,p in enumerate(files))
    body.append('</ul>')
    for i,p in enumerate(files):
        body.append(f'<figure id="fig{i}"><figcaption>{html.escape(p.stem)}</figcaption><a href="{p.relative_to(out).as_posix()}"><img loading="lazy" src="data:image/png;base64,{base64.b64encode(p.read_bytes()).decode()}"></a></figure>')
    body.append('<h2>数值下载</h2><ul>'+''.join(f'<li><a href="{p.as_uri()}">{p.name}</a></li>' for p in numeric.iterdir())+'</ul>')
    (out/'report.html').write_text(''.join(body),encoding='utf-8')
def run_outputs(outputs, atlas, *, configs=None, overrides=None, reference=None, destination=None, redraw=False):
    from srcs.workflow.plugins import function_options
    defaults = function_options(load_config(PROJECT_ROOT / "configs/default.yaml"), "seedfc_supplement")
    overrides = overrides or {}
    overrides = {key: overrides[key] for key in ("visualization", "descriptive", "seedfc_supplement", "plugins")
                 if key in overrides}
    supplement = dict(function_options(overrides, "seedfc_supplement"))
    legacy = {}
    for key in ("missing_color", "series_colors", "view_markers"):
        if key in supplement:
            legacy[key] = supplement.pop(key)
    for old, role in (("cmap", "sequential"), ("signed_cmap", "diverging")):
        if old in supplement:
            legacy.setdefault("colormaps", {})[role] = supplement.pop(old)
    for key in ("dpi", "formats", "panel_size", "font_size"):
        if key in supplement:
            legacy.setdefault("supplement", {})[key] = supplement.pop(key)
    overrides.pop("plugins", None)
    overrides.pop("seedfc_supplement", None)
    overrides["visualization"] = merge(legacy, overrides.get("visualization", {}))
    if configs is None:
        configs = [saved_config(output) for output in outputs]
    configs = [merge(config, overrides) for config in configs]
    if not configs or len(configs) != len(outputs):
        raise ValueError("Expected one configuration per SeedFC output")
    resolved = []
    for config in configs:
        cfg = merge(merge(defaults, function_options(config, "seedfc_supplement")), supplement)
        unknown = set(cfg) - set(defaults) - {"reference"}
        if unknown:
            raise ValueError(f"Unknown supplementary settings: {sorted(unknown)}")
        cfg.pop("enabled", None)
        cfg["reference"] = reference or cfg.get("reference") or config.get("descriptive", {}).get("reference_week")
        effective = settings(config)
        validate(effective)
        resolved.append((cfg, effective))
    cfg, viz = resolved[0]
    if any(other != cfg or style != viz for other, style in resolved[1:]):
        raise ValueError("Supplementary output configurations must agree on analysis and visualization settings")
    if (not cfg["cohorts"] or len(set(cfg["cohorts"])) != len(cfg["cohorts"])
            or not set(cfg["cohorts"]) <= {"fixed_window", "available_pair"}):
        raise ValueError("Invalid supplementary cohort selection")
    if (not np.isfinite([cfg["cv_mean_threshold"], cfg["numerical_tolerance"]]).all()
            or cfg["cv_mean_threshold"] < 0 or cfg["numerical_tolerance"] <= 0
            or type(cfg["icc_min_mice"]) is not int or cfg["icc_min_mice"] < 3):
        raise ValueError("Invalid supplementary numerical settings")
    formats = viz["supplement"]["formats"]
    if "png" not in formats or not set(formats) <= {"png", "pdf"}:
        raise ValueError("Supplementary formats must include PNG and optionally PDF")
    sources = [load_output(path, atlas, config) for path, config in zip(outputs, configs)]
    if len({source["series"] for source in sources}) != len(sources):
        raise ValueError("Distinct SeedFC series labels required")
    first = sources[0]
    for source in sources:
        if (any(source["meta"][key] != first["meta"][key] for key in ("roi_labels", "week_ids", "run_ids"))
                or source["views"] != first["views"] or source["mapping"] != first["mapping"]
                or [(level, name) for level, name, _ in source["groups"]] !=
                   [(level, name) for level, name, _ in first["groups"]]):
            raise ValueError("Incompatible SeedFC ROIs, weeks, Runs, views or networks")
        if cfg["reference"] not in source["meta"]["week_ids"]:
            raise ValueError("Set descriptive.reference_week or --reference to an analysis week")
    cfg["weeks"] = first["meta"]["week_ids"]
    root = Path(destination or (Path(configs[0]["paths"]["output"]) / "_supplement"))
    if not root.is_absolute():
        root = PROJECT_ROOT / root
    if redraw:
        # Enforce a verified numeric hit inside the same lock as regular calculations.
        numeric, record, reused = calculate(sources, cfg, root, reuse=True, require_existing=True)
    else:
        numeric, record, reused = calculate(sources, cfg, root, reuse=True)
    report, reused_report = render_report(sources, cfg, viz, numeric, record, root)
    return dict(atlas=atlas, path=str(report), numeric=str(numeric), reused_numeric=reused,
                reused_report=reused_report, manifest_sha256=sha(report / "manifest.json"))


def run_batch(outputs, contexts):
    atlases = contexts[0].atlases
    if any(ctx.atlases != atlases for ctx in contexts[1:]):
        raise ValueError("Supplementary reports require identical configured atlas partitions")
    for atlas in atlases:
        yield run_outputs(outputs, "_".join(atlas), configs=[ctx.config for ctx in contexts],
                          destination=Path(contexts[0].output_root) / "_supplement")


def run(args):
    overrides = load_config(args.config) if args.config else None
    atlas = args.atlas
    if atlas is None:
        atlas = "_".join(saved_config(args.outputs[0])["variables"]["Atlas"][0])
    result = run_outputs(args.outputs, atlas, overrides=overrides, reference=args.reference,
                         destination=args.destination, redraw=getattr(args, "redraw", False))
    print("RESULT", result["path"], flush=True)
    return Path(result["path"])


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outputs',nargs='+',required=True,type=Path)
    parser.add_argument('--atlas', help='Atlas partition; defaults to the first configured partition')
    parser.add_argument('--reference',help='Reference label; default from effective settings')
    parser.add_argument('--config',type=Path,help='YAML containing seedfc_supplement settings')
    parser.add_argument('--destination',type=Path,help='Root directory for cached supplementary results')
    parser.add_argument("--redraw", action="store_true", help="Require existing verified numerical results")
    run(parser.parse_args())
