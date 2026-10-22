"""Stage orchestration; one worker process per config isolates library state."""
import argparse
from datetime import datetime, timezone
import importlib
from pathlib import Path

from srcs.utils.config import load_context
from srcs.utils.fc_import import import_fc
from srcs.utils.results import save_json
from srcs.workflow import STAGES, MODULES

def validate(ctx, stage, importing=False):
    _ = ctx.mice, ctx.weeks, ctx.runs, ctx.atlases, ctx.split_axis, ctx.selection_mode
    if stage in ('preprocess', 'fc'):
        from srcs.utils.dataset import resolve_recording
        found = 0
        for week in ctx.weeks:
            for mouse in ctx.mice:
                for run in ctx.runs:
                    try:
                        resolve_recording(ctx, week, mouse, run, raw=stage == 'preprocess')
                        found += 1
                    except FileNotFoundError:
                        pass
        print(f'Available recording inputs: {found}/{len(ctx.weeks)*len(ctx.mice)*len(ctx.runs)}')
        if not found:
            raise FileNotFoundError('No inputs found; configure paths or explicitly import FC with --start edges')
    elif stage == 'edges' and not importing:
        from srcs.utils.edge_store import source_identity
        for atlas in ctx.atlases:
            source_identity(ctx, '_'.join(atlas))
    elif stage in ('pairwise', 'descriptive', 'combined'):
        from srcs.utils.edge_store import load_store
        for atlas in ctx.atlases:
            store = load_store(ctx, '_'.join(atlas))
            print(f"{store.metadata['atlas']}: F {store.F_r.shape}; valid {store.valid_record.sum()}/{store.valid_record.size}")
        if stage == 'pairwise':
            from srcs.analysis.pairwise import selections, week_positions
            selections(ctx)
            week_positions(ctx)
        if stage == 'combined':
            from srcs.utils.dataset import select_complete_cohort
            mice, audit = select_complete_cohort(ctx)
            if len(mice) < 2:
                raise ValueError('Fewer than two complete mice remain')
            ctx.group_indices(audit['retained_weeks'], ctx.runs)
            print(f"Combined cohort: {len(mice)} mice x {len(audit['retained_weeks'])} weeks")
    elif stage == 'statistics':
        import pandas as pd
        from srcs.utils.dataset import statistics_metadata
        statistics_metadata(ctx, pd.read_pickle(ctx.statistics_path))
    elif stage == 'report':
        from srcs.workflow.report import report_inputs
        inputs = report_inputs(ctx)
        section = getattr(ctx.args, 'report_section', 'all')
        paths = list(inputs.values()) if section == 'all' else [inputs[section]]
        if not any(path.is_file() for path in paths):
            raise FileNotFoundError(f'No saved numerical input for report section {section}')


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', required=True)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--stage', choices=STAGES)
    selection.add_argument('--start', choices=STAGES, default='fc')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--import-fc')
    parser.add_argument('--report-section', choices=('all', 'pairwise', 'combined', 'statistics'), default='all')
    args = parser.parse_args(argv)
    ctx = load_context(args.config, create_dirs=False)
    ctx.args.report_section = args.report_section
    stages = [args.stage] if args.stage else list(STAGES[STAGES.index(args.start):])
    if args.import_fc and stages[0] != 'edges':
        parser.error('--import-fc requires --stage/--start edges')
    print(f'Output: {ctx.output_path}\nStages: {", ".join(stages)}', flush=True)
    if args.import_fc:
        details = import_fc(ctx, args.import_fc, dry_run=args.dry_run)
        print(f'FC import verified: {len(details["tables"])} tables; {details["audit_source"]}', flush=True)
    validate(ctx, stages[0], importing=bool(args.import_fc))
    if args.dry_run:
        return 0
    if not (Path(ctx.output_path) / "metadata/analysis_config.json").exists():
        ctx.save_analysis_metadata()
    for stage in stages:
        if stage == 'fc' and (Path(ctx.output_path) / 'metadata/fc_import_provenance.json').exists():
            raise ValueError('This output contains imported FC; start edges, or use a separate config to recompute FC')
        marker = Path(ctx.output_dir('Manifests')) / f'{stage}.json'
        state = dict(stage=stage, config_tag=ctx.analysis_tag, status='running',
                     started=datetime.now(timezone.utc).isoformat())
        save_json(marker, state)
        ctx.logger.info('Starting stage: %s', stage)
        try:
            importlib.import_module(MODULES[stage]).run(ctx)
        except BaseException as exc:
            state.update(status='interrupted' if isinstance(exc, KeyboardInterrupt) else 'failed', reason=str(exc))
            save_json(marker, state)
            raise
        state.update(status='complete', finished=datetime.now(timezone.utc).isoformat())
        save_json(marker, state)
        ctx.logger.info('Completed stage: %s', stage)
    return 0
