"""Explicit staged workflow with verified outputs and durable record audits."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from srcs.utils.config import load_context
from srcs.utils.output_paths import fc_reports
from srcs.utils.progress import emit_progress, progress_items
from srcs.utils.io import save_json, save_pickle, save_table

STAGES=['preprocess','fc','reliability','weekly','longitudinal','report','phenotype']


def file_identity(path):
    p=Path(path)
    return dict(path=str(p.resolve()),sha256=hashlib.sha256(p.read_bytes()).hexdigest())


def code_identity():
    root=Path(__file__).resolve().parents[1]
    digest=hashlib.sha256()
    for p in sorted(root.rglob('*.py')):
        if 'snapshot' not in p.parts:
            digest.update(p.relative_to(root).as_posix().encode());digest.update(p.read_bytes())
    return digest.hexdigest()


class StageRunner:
    def __init__(self,ctx):
        self.ctx=ctx
        self.root=Path(ctx.output_dir('Manifests'))

    def run(self,name,action,outputs,dependencies=()):
        emit_progress(phase='checking saved stage')
        outputs=[Path(p) for p in outputs]
        marker=self.root/(name+'.json')
        identity=dict(config_tag=self.ctx.analysis_tag,code_sha256=code_identity(),
                      inputs=[file_identity(p) for p in dependencies])
        if marker.exists():
            try:
                manifest=json.loads(marker.read_text(encoding='utf-8'))
                if manifest.get('status')=='complete' and manifest.get('identity')==identity and manifest.get('outputs')==[file_identity(p) for p in outputs]:
                    self.ctx.logger.info('Skipping complete stage: %s',name)
                    return 'reused'
            except (OSError,ValueError,KeyError):pass
        save_json(marker,dict(status='running',identity=identity))
        try:
            action()
            emit_progress(phase='saving output manifest')
            committed=[file_identity(p) for p in outputs]
        except BaseException as exc:
            save_json(marker,dict(status='interrupted' if isinstance(exc,KeyboardInterrupt) else 'failed',identity=identity,reason=str(exc)))
            raise
        save_json(marker,dict(status='complete',identity=identity,outputs=committed))
        return 'success'


def run_preprocess(ctx):
    from srcs.preprocessing.workflow import process_recording
    from srcs.utils.batch import BatchReport
    report=BatchReport(ctx,'01_preprocess',total=len(ctx.weeks)*len(ctx.mice)*len(ctx.runs))
    for w in ctx.weeks:
        for m in ctx.mice:
            for r in ctx.runs:
                report.run(lambda:process_recording(ctx,w,m,r),week=w,mouse=m,run=r)
    report.finish()


def run_fc(ctx):
    from srcs.analysis.connectivity import CorrelationCache
    from srcs.utils.batch import BatchReport
    from srcs.utils.dataset import table_path
    report=BatchReport(ctx,'02_correlation',total=len(ctx.atlases)*len(ctx.weeks)*len(ctx.mice)*len(ctx.runs))
    for atlas in ctx.atlases:
        tag='_'.join(atlas)
        try:cache=CorrelationCache(ctx,atlas)
        except Exception as exc:
            cache=None;initialization_error=exc
        for w in ctx.weeks:
            frame=pd.DataFrame(index=ctx.mice,columns=ctx.runs,dtype=object)
            for m in ctx.mice:
                for r in ctx.runs:
                    def action():
                        if cache is None:raise initialization_error
                        matrix,status=cache.process(w,m,r)
                        frame.at[m,r]=matrix
                        return dict(status=status,input_path=matrix.attrs['correlation_signature']['source'][0])
                    report.run(action,week=w,mouse=m,run=r,atlas=tag)
            frame.attrs['analysis']=dict(config_tag=ctx.analysis_tag,week=w,atlas=tag)
            save_pickle(table_path(ctx,w,tag),frame)
    report.finish()
    # A missing or incomplete report is not accepted as a completed FC stage.
    if not report.persisted:raise OSError('FC batch report was not saved')


def load_records(ctx):
    """Cross-check every PKL cell with the current-config latest complete report."""
    from srcs.analysis.connectivity import valid_fc
    from srcs.utils.dataset import table_path
    reports=fc_reports(ctx.output_path)
    if not reports:raise FileNotFoundError('Current-config FC audit required; run --start fc after clearing outputs, or --import-fc for old outputs')
    report_path=max(reports,key=lambda p:(p.stat().st_mtime_ns,p.name))
    payload=json.loads(report_path.read_text(encoding='utf-8'))
    expected={(a,w,m,r) for a in ['_'.join(a) for a in ctx.atlases] for w in ctx.weeks for m in ctx.mice for r in ctx.runs}
    report={}
    for row in payload['records']:
        key=tuple(row[k] for k in ['atlas','week','mouse','run'])
        if key in report:raise ValueError('Duplicate FC report key: '+str(key))
        report[key]=row
    if set(report)!=expected:raise ValueError('Incomplete or mismatched FC report scope')
    records={};audit=[];inputs=[report_path]
    for atlas in ['_'.join(a) for a in ctx.atlases]:
        labels=None
        for w in progress_items(ctx.weeks, phase=f'{atlas}: checking FC inputs'):
            path=table_path(ctx,w,atlas)
            frame=pd.read_pickle(path);inputs.append(path)
            if not isinstance(frame,pd.DataFrame) or not frame.index.is_unique or not frame.columns.is_unique:raise ValueError('Invalid FC table: '+str(path))
            meta=frame.attrs.get('analysis',{})
            if meta.get('config_tag')!=ctx.analysis_tag:raise ValueError('FC table belongs to another configuration: '+str(path))
            for m in ctx.mice:
                for r in ctx.runs:
                    key=(atlas,w,m,r);entry=report[key]
                    cell=frame.at[m,r] if m in frame.index and r in frame.columns else None
                    reason=''
                    if entry['status'] not in ['success','reused']:reason=entry.get('reason') or entry['status']
                    elif not isinstance(cell,pd.DataFrame) or not valid_fc(cell,list(cell.index)):reason='missing/invalid PKL cell'
                    elif len(cell)<3:reason='fewer than three ROIs for fingerprint analysis'
                    elif labels is not None and list(cell.index)!=labels:reason='ROI labels/order mismatch'
                    if not reason:
                        labels=list(cell.index);records[key]=cell
                    audit.append(dict(atlas=atlas,week=w,mouse=m,run=r,status='invalid' if reason else 'valid',reason=reason,
                        missing_kind=ctx.section('study').get('missing_reasons',{}).get(f'{m}/{w}/{r}','unknown') if reason else '',input_path=entry.get('input_path')))
    return records,pd.DataFrame(audit),inputs


def import_fc(ctx,source):
    """Explicit, audited migration of old FC tables; never automatic cross-config reuse."""
    from srcs.utils.dataset import table_path
    from srcs.analysis.connectivity import valid_fc
    source=Path(source)
    config_path=source/'metadata/analysis_config.json'
    if not config_path.is_file():config_path=source/'analysis_config.json'
    if not config_path.is_file():raise FileNotFoundError('Import requires source analysis_config.json')
    source_config=json.loads(config_path.read_text(encoding='utf-8'))['config']
    # Old output config may differ in ID/longitudinal options; FC inputs must match.
    if source_config.get('variables')!=ctx.config.get('variables'):raise ValueError('Import mismatch: variables')
    def effective_preprocess(config):
        options=dict(config.get('preprocess',{}))
        for name in ['filter','bin','robust','resume']:options.pop(name,None)
        options.setdefault('signal',dict(input_fps=30,output_fps=5,highpass=0.009,lowpass=0.08))
        return options
    if effective_preprocess(source_config)!=effective_preprocess(ctx.config):raise ValueError('Import mismatch: effective preprocessing options')
    for key in ['mask','processed']:
        if source_config['paths'][key]!=ctx.paths[key]:raise ValueError('Import mismatch: paths.'+key)
    reports=fc_reports(source)
    if not reports:raise FileNotFoundError('Import requires source 02 BatchReport')
    report=max(reports,key=lambda p:(p.stat().st_mtime_ns,p.name))
    payload=json.loads(report.read_text(encoding='utf-8'))
    expected={(w,m,r,a) for w in ctx.weeks for m in ctx.mice for r in ctx.runs for a in ['_'.join(a) for a in ctx.atlases]}
    keys=[tuple(row[k] for k in ['week','mouse','run','atlas']) for row in payload['records']]
    if len(keys)!=len(set(keys)) or set(keys)!=expected:raise ValueError('Import report scope mismatch')
    provenance=[file_identity(config_path),file_identity(report)]
    prepared=[]
    for atlas in ['_'.join(a) for a in ctx.atlases]:
        for w in ctx.weeks:
            path=source/'Correlation'/f'{w}_correlation_matrices_{atlas}_excluded.pkl'
            frame=pd.read_pickle(path)
            if not isinstance(frame,pd.DataFrame) or list(frame.index)!=ctx.mice or list(frame.columns)!=ctx.runs:raise ValueError('Import table dimensions mismatch')
            for m in ctx.mice:
                for r in ctx.runs:
                    cell=frame.at[m,r]
                    row=payload['records'][keys.index((w,m,r,atlas))]
                    if row['status'] in ['success','reused'] and (not isinstance(cell,pd.DataFrame) or not valid_fc(cell,list(cell.index))):raise ValueError('Successful source cell is invalid')
            frame.attrs['analysis']=dict(config_tag=ctx.analysis_tag,week=w,atlas=atlas,imported=True)
            prepared.append((table_path(ctx,w,atlas),frame));provenance.append(file_identity(path))
    target_reports=Path(ctx.output_dir('BatchReports'))
    if list(target_reports.glob('02_correlation_*.json')):raise FileExistsError('Target already has FC reports; import into an unused configuration output')
    for path,frame in prepared:save_pickle(path,frame)
    save_json(target_reports/'02_correlation_imported.json',payload)
    save_json(Path(ctx.output_path,'metadata/fc_import_provenance.json'),dict(source=str(source.resolve()),files=provenance,verification='explicit matching acquisition/preprocessing config; source report and FC cells verified'))


def verified_report_inputs(ctx, fc_inputs):
    """Accept relocated saved results only when their config/input/output hashes match.

    Reporting is allowed to use results produced before a plotting-code change.
    No inference is recomputed here, and stale or altered numerical files fail closed.
    """
    root = Path(ctx.output_path)
    dependencies = []
    expected = sorted((Path(p).name, file_identity(p)['sha256']) for p in fc_inputs)
    for stage, files in (
        ('03_reliability', ['results/reliability/within_week.csv', 'results/reliability/comparisons.json']),
        ('05_longitudinal', ['results/longitudinal/' + ctx.statistics_file, 'results/longitudinal/comparisons.json',
                             'results/longitudinal/cohort.json', 'results/longitudinal/paired_mice.csv']),
    ):
        marker = root / 'metadata/manifests' / (stage + '.json')
        if not marker.is_file():
            raise FileNotFoundError(f'Report requires saved {stage} results; run --start reliability first')
        manifest = json.loads(marker.read_text(encoding='utf-8'))
        identity = manifest.get('identity', {})
        saved_inputs = sorted((str(row['path']).replace('\\', '/').rsplit('/', 1)[-1], row['sha256'])
                              for row in identity.get('inputs', []))
        if manifest.get('status') != 'complete' or identity.get('config_tag') != ctx.analysis_tag or saved_inputs != expected:
            raise ValueError(f'Report input provenance mismatch in {stage}; rerun --start reliability')
        saved_outputs = {str(row['path']).replace('\\', '/').rsplit('/', 1)[-1]: row['sha256']
                         for row in manifest.get('outputs', [])}
        for name in files:
            path = root / name
            if not path.is_file() or saved_outputs.get(path.name) != file_identity(path)['sha256']:
                raise ValueError(f'Report numerical output missing or changed: {path}; rerun {stage}')
            dependencies.append(path)
        dependencies.append(marker)
    return dependencies


def run_analysis(ctx,stage='all',start='reliability',plot_workers=1):
    from srcs.analysis.reliability import analyze_reliability
    from srcs.analysis.longitudinal import weekly_representatives,analyze_between_weeks,trajectories
    if stage not in STAGES+['all']:raise ValueError('Unknown stage')
    selected=STAGES[STAGES.index(start):] if stage=='all' else [stage]
    if 'preprocess' in selected:run_preprocess(ctx)
    if 'fc' in selected:run_fc(ctx)
    analyses=[s for s in selected if s not in ['preprocess','fc']]
    if not analyses:return
    records,audit,inputs=load_records(ctx)
    save_table(Path(ctx.output_path,'metadata/record_audit.csv'),audit,index=False)
    runner=StageRunner(ctx)
    weekly,qc=weekly_representatives(ctx,records)
    if 'reliability' in analyses:
        folder=Path(ctx.output_dir('Reliability'))
        def reliability():
            rows,result=analyze_reliability(ctx,records)
            save_table(folder/'within_week.csv',rows,index=False);save_json(folder/'comparisons.json',result)
        runner.run('03_reliability',reliability,[folder/'within_week.csv',folder/'comparisons.json'],inputs)
    if 'weekly' in analyses:
        folder=Path(ctx.output_dir('Weekly'))
        def representative():
            save_pickle(folder/'representatives.pkl',weekly);save_table(folder/'quality.csv',qc,index=False)
        runner.run('04_weekly',representative,[folder/'representatives.pkl',folder/'quality.csv'],inputs)
    if 'longitudinal' in analyses:
        folder=Path(ctx.output_dir('Identification'))
        def longitudinal():
            results,rows,cohort=analyze_between_weeks(ctx,weekly)
            save_pickle(folder/ctx.statistics_file,results);save_json(folder/'comparisons.json',results)
            save_json(folder/'cohort.json',cohort);save_table(folder/'paired_mice.csv',rows,index=False)
        runner.run('05_longitudinal',longitudinal,[folder/ctx.statistics_file,folder/'comparisons.json',folder/'cohort.json',folder/'paired_mice.csv'],inputs)
    if 'report' in analyses:
        folder=Path(ctx.output_dir('Statistics'))
        control_path=ctx.section('longitudinal').get('control_representatives')
        report_inputs=verified_report_inputs(ctx,inputs)
        from srcs.visualization.report_plots import build_jobs,render_jobs,structure_metrics
        results=json.loads(Path(ctx.output_path,'results/longitudinal/comparisons.json').read_text(encoding='utf-8'))
        reliability_results=json.loads(Path(ctx.output_path,'results/reliability/comparisons.json').read_text(encoding='utf-8'))
        frame=trajectories(ctx,weekly)
        jobs=build_jobs(ctx,records,results,reliability_results,frame)
        figure_index=folder/'figure_index.json'
        def report():
            save_table(folder/'tables/trajectories.csv',frame,index=False)
            save_table(folder/'tables/structure_metrics.csv',structure_metrics(ctx,results),index=False)
            render_jobs(ctx,jobs,figure_index,workers=plot_workers)
            summary=dict(status='complete' if records else 'no_valid_data',n_valid_recordings=len(records),n_invalid_recordings=int((audit.status!='valid').sum()),
                n_weekly_observations=len(weekly),trajectory_population='all valid observations',
                observation_start='first valid measurement per mouse; not assumed healthy baseline',
                inference='Descriptive; overlapping week pairs are not independent samples',
                phenotype='unavailable unless supplied',group=ctx.section('study').get('group','unspecified'))
            save_json(folder/'summary.json',summary)
            if control_path:
                from srcs.analysis.longitudinal import control_deviation
                controls=pd.read_pickle(control_path)
                if not isinstance(controls,dict):raise ValueError('Control representatives must be a weekly record dictionary')
                save_table(folder/'tables/control_deviation.csv',control_deviation(ctx,weekly,controls),index=False)
        report_outputs=[folder/'tables/trajectories.csv',folder/'tables/structure_metrics.csv',folder/'summary.json',figure_index]
        if ctx.render_options()['save']:
            report_outputs.extend(Path(ctx.output_path)/job.path for job in jobs)
        if control_path:report_outputs.append(folder/'tables/control_deviation.csv')
        runner.run('06_report',report,report_outputs,inputs+report_inputs+([control_path] if control_path else []))
    if 'phenotype' in analyses:
        folder=Path(ctx.output_dir('Phenotype'));options=ctx.section('phenotype')
        if not options.get('observations') and not options.get('events'):
            save_json(folder/'status.json',dict(status='unavailable',reason='No phenotype or event inputs supplied'))
        else:
            from srcs.analysis.phenotype import read_observations,align_observations,read_events
            paths=[options[k] for k in ['observations','events'] if options.get(k)]
            outputs=[folder/'status.json']
            if options.get('observations'):outputs.append(folder/'aligned_observations.csv')
            if options.get('events'):outputs.append(folder/'events.csv')
            def phenotype():
                if options.get('observations'):
                    obs=read_observations(options['observations'],ctx.mice)
                    save_table(folder/'aligned_observations.csv',align_observations(trajectories(ctx,weekly),obs,options.get('tolerance_days',7)),index=False)
                if options.get('events'):save_table(folder/'events.csv',read_events(options['events'],ctx.mice),index=False)
                save_json(folder/'status.json',dict(status='inputs_validated',model_status='not_fitted'))
            runner.run('07_phenotype',phenotype,outputs,inputs+paths)


def main(argv=None,default_stage='all'):
    parser=argparse.ArgumentParser(description='ISOI longitudinal workflow. Run from project root for relative paths.')
    parser.add_argument('-c','--config',default='configs/default.yaml')
    parser.add_argument('--stage',choices=STAGES+['all','legacy'],default=default_stage)
    parser.add_argument('--start',choices=STAGES,default='preprocess')
    parser.add_argument('--import-fc',metavar='SOURCE_OUTPUT')
    parser.add_argument('--plot-workers',type=int,default=4)
    args=parser.parse_args(argv)
    ctx=load_context(args.config,create_dirs=True)
    save_json(Path(ctx.output_path,'metadata/analysis_config.json'),dict(config=ctx.config,config_tag=ctx.analysis_tag,code_sha256=code_identity()))
    try:
        if args.import_fc:import_fc(ctx,args.import_fc)
        if args.stage=='legacy':
            from srcs.legacy.workflow import run_legacy
            run_legacy(ctx)
        else:
            run_analysis(ctx,args.stage,args.start,plot_workers=args.plot_workers)
    except (ValueError,OSError) as exc:
        ctx.logger.error('Workflow stopped: %s',exc)
        return 1
    return 0


if __name__=='__main__':sys.exit(main())
