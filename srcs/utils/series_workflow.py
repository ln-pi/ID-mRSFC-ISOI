"""Sequential A/N series workflow starting from existing FiltDRR recordings."""
from collections import Counter
from contextlib import contextmanager
from datetime import datetime
import json
import os
from pathlib import Path
import sys
from uuid import uuid4

from srcs.utils.config import load_context
from srcs.utils.dataset import resolve_recording
from srcs.utils.io import save_json
from srcs.utils.logger import setup_logger
from srcs.utils.live_console import StageDisplay, execute_stage

STAGES = ('preprocess', 'fc', 'reliability', 'weekly', 'longitudinal', 'report', 'phenotype')


@contextmanager
def project_directory(root):
    previous = Path.cwd()
    os.chdir(root)
    try:
        yield
    finally:
        os.chdir(previous)


def preflight(ctx, start):
    """Check source paths/masks without loading or calculating the movies."""
    from srcs.analysis.connectivity import CorrelationCache
    if start == 'preprocess':
        raw = Path(ctx.paths['raw'])
        if not raw.is_dir():
            raise FileNotFoundError(f'Raw input directory not found: {raw}')
        return dict(raw_available=True, raw_path=str(raw.resolve()), requires_manual_registration=True)
    if start != 'fc':
        from srcs.utils.pipeline import load_records
        records, audit, inputs = load_records(ctx)
        if start == 'report':
            from srcs.utils.pipeline import verified_report_inputs
            verified_report_inputs(ctx, inputs)
        return dict(valid_fc=len(records), invalid_fc=int((audit.status != 'valid').sum()))
    masks = []
    for atlas in ctx.atlases:
        cache = CorrelationCache(ctx, atlas)
        if len(cache.names) < 3:
            raise ValueError(f'Fewer than three retained ROIs: {atlas}')
        masks.append(dict(atlas='_'.join(atlas), n_rois=len(cache.names), shape=cache.shape))
    inputs = []
    for week in ctx.weeks:
        for mouse in ctx.mice:
            for run in ctx.runs:
                row = dict(week=week, mouse=mouse, run=run)
                try:
                    row.update(status='found', path=str(resolve_recording(ctx, week, mouse, run)))
                except (OSError, ValueError) as exc:
                    row.update(status='missing' if isinstance(exc, FileNotFoundError) else 'invalid', reason=str(exc))
                inputs.append(row)
    return dict(masks=masks, input_counts=dict(Counter(r['status'] for r in inputs)), inputs=inputs)


def inspect_stage(ctx, stage):
    """A zero process exit code does not imply scientifically available results."""
    root = Path(ctx.output_path)
    if stage == 'fc':
        from srcs.utils.pipeline import load_records
        records, audit, _ = load_records(ctx)
        invalid = int((audit.status != 'valid').sum())
        return dict(status='completed_with_missing' if records and invalid else 'complete' if records else 'no_valid_data',
                    valid_fc=len(records), invalid_fc=invalid)
    if stage == 'longitudinal':
        cohort = json.loads((root / 'results/longitudinal/cohort.json').read_text(encoding='utf-8'))
        results = json.loads((root / 'results/longitudinal/comparisons.json').read_text(encoding='utf-8'))
        complete = sum(r.get('status') == 'complete' for r in results.values())
        return dict(status='complete' if cohort['status'] == 'complete' and complete else 'unavailable',
                    cohort_status=cohort['status'], n_mice=len(cohort['mice']), weeks=cohort['weeks'],
                    complete_comparisons=complete, unavailable_comparisons=len(results)-complete)
    if stage == 'phenotype':
        result = json.loads((root / 'results/phenotype/status.json').read_text(encoding='utf-8'))
        return dict(status='skipped' if result['status'] == 'unavailable' else 'complete', detail=result)
    return dict(status='complete')


def run_batch(configs, *, project_root=None, start='fc', stage=None, dry_run=False,
              report_root='reports/automatic', execute=None, import_source=None, plot_workers=4):
    if start not in STAGES:
        raise ValueError('Unknown start stage: ' + start)
    project_root = Path(project_root or Path(__file__).resolve().parents[2]).resolve()
    execute = execute or execute_stage
    if stage is not None and stage not in (*STAGES, 'legacy'):
        raise ValueError('Unknown stage: ' + stage)
    stages = (stage,) if stage else STAGES[STAGES.index(start):]
    start = stages[0]
    with project_directory(project_root):
        stamp = datetime.now().strftime('%Y%m%d-%H%M%S')
        folder = Path(report_root).resolve() / f'{stamp}_{uuid4().hex[:8]}'
        folder.mkdir(parents=True, exist_ok=True)
        report_path = folder / 'batch_summary.json'
        logger = None
        logfile = None
        display = None
        report = dict(status='running', dry_run=dry_run, project_root=str(project_root),
                      start=start, stages=list(stages), python=sys.executable, series=[])

        def flush():
            save_json(report_path, report)

        def finish_batch():
            flush()
            if logger is not None:
                logger.info('Project status: %s; summary: %s', report['status'], report_path)
            print(f'Project ... {report["status"]}', flush=True)
            for item in report['series']:
                if item.get('log'):
                    print(f'Log ({Path(item["config"]).stem}): {item["log"]}', flush=True)
            close_log()

        def close_log():
            if logger is not None:
                for handler in list(logger.handlers):
                    handler.flush()
                    handler.close()
                    logger.removeHandler(handler)

        flush()
        try:
            for index, config in enumerate(configs, 1):
                close_log()
                logger = None
                logfile = None
                entry = dict(config=str(Path(config).resolve()), status='running', stages=[])
                report['series'].append(entry)
                flush()
                display = None
                try:
                    ctx = load_context(config)
                    logfile = Path(ctx.paths['log']).resolve()
                    entry['log'] = str(logfile)
                    logger = setup_logger(str(logfile), to_console=False, rotate=False)
                    logger.info('=== Series %s; batch %s; stages=%s; dry_run=%s ===',
                                Path(config).stem, folder.name, stages, dry_run)
                    if import_source:
                        from srcs.utils.pipeline import import_fc
                        import_fc(ctx, import_source)
                    entry.update(config_tag=ctx.analysis_tag, output_path=str(Path(ctx.output_path).resolve()))
                    entry['preflight'] = preflight(ctx, start)
                    available = (entry['preflight']['raw_available'] if start == 'preprocess' else
                                 entry['preflight']['input_counts'].get('found', 0) if start == 'fc'
                                 else entry['preflight']['valid_fc'])
                    if not available:
                        entry.update(status='blocked', reason='No available FiltDRR/valid FC inputs; downstream analysis not started')
                        logger.warning('%s: %s', Path(config).stem, entry['reason'])
                        print(f'{Path(config).stem} ... blocked: {entry["reason"]}', flush=True)
                        flush()
                        continue
                    print(f'[{index}/{len(configs)}] {Path(config).stem}', flush=True)
                    for stage in stages:
                        command = [sys.executable, '-u', str(Path(__file__).resolve().parents[2] / 'run.py'), '--worker',
                                   '-c', str(Path(config).resolve()), '--stage', stage, '--plot-workers', str(plot_workers)]
                        step = dict(stage=stage, command=command, log=str(logfile), status='planned' if dry_run else 'running')
                        entry['stages'].append(step)
                        flush()
                        if dry_run:
                            print(f'{stage} ... planned', flush=True)
                            continue
                        logger.info('--- %s / %s START ---', Path(config).stem, stage)
                        display = StageDisplay(stage)
                        result = execute(command, cwd=str(project_root), logger=logger, display=display,
                                         env=dict(os.environ, PYTHONUNBUFFERED='1', PYTHONIOENCODING='utf-8',
                                                  ISOI_MANAGED_WORKER='1'))
                        step['warning_count'] = getattr(result, 'warning_count', 0)
                        step['returncode'] = result.returncode
                        if result.returncode:
                            step['status'] = 'failed'
                            display.finish('failed')
                            display = None
                            logger.error('--- %s / %s FAILED (exit %s) ---', Path(config).stem, stage, result.returncode)
                            entry.update(status='failed', reason=f'{stage} failed; inspect {logfile}')
                            flush()
                            break
                        step.update(inspect_stage(ctx, stage))
                        if step['warning_count'] and step['status'] == 'complete':
                            step['status'] = 'completed_with_warnings'
                        display.finish(step['status'])
                        display = None
                        logger.info('--- %s / %s %s ---', Path(config).stem, stage, step['status'])
                        flush()
                        if step['status'] == 'no_valid_data':
                            entry.update(status='blocked', reason='FC produced no valid matrices; downstream analysis not started')
                            break
                    else:
                        warnings = (entry['preflight'].get('invalid_fc', 0) > 0
                                    or any(s['status'] in ('completed_with_missing', 'completed_with_warnings', 'unavailable') for s in entry['stages']))
                        entry['status'] = 'planned' if dry_run else 'completed_with_warnings' if warnings else 'complete'
                except Exception as exc:
                    if display is not None:
                        display.finish('failed')
                        display = None
                    if logger is None:
                        logfile = Path('logs').resolve() / (Path(config).stem + '.log')
                        entry['log'] = str(logfile)
                        logger = setup_logger(str(logfile), to_console=False, rotate=False)
                    logger.exception('%s failed', Path(config).stem)
                    entry.update(status='failed', reason=str(exc), error_type=type(exc).__name__)
                    if entry['stages'] and entry['stages'][-1]['status'] == 'running':
                        entry['stages'][-1].update(status='failed', reason=str(exc))
                    print(f'  {Path(config).stem} failed: {exc}', file=sys.stderr, flush=True)
                flush()
                logger.info('Series status: %s', entry['status'])
                print(f'{Path(config).stem} ... {entry["status"]}', flush=True)
        except KeyboardInterrupt:
            if display is not None:
                display.finish('interrupted')
            if logger is not None:
                logger.warning('Run interrupted by user')
            report['status'] = 'interrupted'
            if report['series'] and report['series'][-1]['status'] == 'running':
                entry = report['series'][-1]
                entry['status'] = 'interrupted'
                if entry['stages'] and entry['stages'][-1]['status'] == 'running':
                    entry['stages'][-1]['status'] = 'interrupted'
            finish_batch()
            return 130, report_path
        statuses = [s['status'] for s in report['series']]
        failed = any(s in ('failed', 'blocked') for s in statuses)
        warnings = 'completed_with_warnings' in statuses
        report['status'] = 'completed_with_errors' if failed else 'planned' if dry_run else 'completed_with_warnings' if warnings else 'complete'
        report['counts'] = dict(Counter(statuses))
        finish_batch()
        return (1 if failed else 0), report_path

