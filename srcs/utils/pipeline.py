"""Stage orchestration; one worker process per config isolates library state."""
import argparse
from datetime import datetime, timezone
import importlib
from pathlib import Path

from srcs.utils.config import load_context
from srcs.utils.fc_import import import_fc
from srcs.utils.results import save_json

STAGES = ('preprocess', 'fc', 'identification', 'statistics')
MODULES = dict(preprocess='preprocess', fc='correlation', identification='identification', statistics='statistics')


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', required=True)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--stage', choices=STAGES)
    selection.add_argument('--start', choices=STAGES, default='fc')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--import-fc')
    args = parser.parse_args(argv)
    ctx = load_context(args.config, create_dirs=False)
    stages = [args.stage] if args.stage else list(STAGES[STAGES.index(args.start):])
    print(f'Output: {ctx.output_path}\nStages: {", ".join(stages)}', flush=True)
    if args.import_fc:
        details = import_fc(ctx, args.import_fc, dry_run=args.dry_run)
        print(f'FC import verified: {len(details["tables"])} tables; {details["audit_source"]}', flush=True)
    if args.dry_run:
        if stages[0] == 'identification' and not args.import_fc:
            from srcs.utils.dataset import select_complete_cohort
            mice, audit = select_complete_cohort(ctx)
            if len(mice) < 2:
                raise ValueError('Fewer than two complete mice remain')
            ctx.group_indices(audit['retained_weeks'], ctx.runs)
            print(f'Cohort: {len(mice)} mice x {len(audit["retained_weeks"])} weeks x {len(ctx.runs)} runs')
        elif stages[0] == 'statistics':
            import pandas as pd
            from srcs.utils.dataset import statistics_metadata
            statistics_metadata(ctx, pd.read_pickle(ctx.statistics_path))
        elif stages[0] in ('fc', 'preprocess'):
            from srcs.utils.dataset import resolve_recording
            found = 0
            for w in ctx.weeks:
                for m in ctx.mice:
                    for r in ctx.runs:
                        try:
                            resolve_recording(ctx, w, m, r, raw=stages[0] == 'preprocess')
                            found += 1
                        except FileNotFoundError:
                            pass
            print(f'Available recording inputs: {found}/{len(ctx.weeks)*len(ctx.mice)*len(ctx.runs)}')
            if not found:
                raise FileNotFoundError('No recording inputs found; configure paths or import P2 FC')
        return 0
    if not (Path(ctx.output_path) / "metadata/analysis_config.json").exists():
        ctx.save_analysis_metadata()
    for stage in stages:
        if stage == 'fc' and (Path(ctx.output_path) / 'metadata/fc_import_provenance.json').exists():
            raise ValueError('This output contains imported FC; start identification, or use a separate config to recompute FC')
        marker = Path(ctx.output_dir('Manifests')) / f'{stage}.json'
        state = dict(stage=stage, config_tag=ctx.analysis_tag, status='running',
                     started=datetime.now(timezone.utc).isoformat())
        save_json(marker, state)
        ctx.logger.info('Starting stage: %s', stage)
        try:
            importlib.import_module('srcs.workflow.' + MODULES[stage]).run(ctx)
        except BaseException as exc:
            state.update(status='interrupted' if isinstance(exc, KeyboardInterrupt) else 'failed', reason=str(exc))
            save_json(marker, state)
            raise
        state.update(status='complete', finished=datetime.now(timezone.utc).isoformat())
        save_json(marker, state)
        ctx.logger.info('Completed stage: %s', stage)
    return 0
