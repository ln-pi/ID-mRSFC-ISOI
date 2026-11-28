"""Stage orchestration; one worker process per config isolates library state."""
import argparse
from datetime import datetime, timezone
import importlib
import json
from pathlib import Path

from srcs.utils.config import load_context
from srcs.utils.fc_import import import_fc
from srcs.utils.storage import save_json
from srcs.core.runtime import select_steps

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
    elif stage in ('pairwise', 'descriptive', 'stability', 'longitudinal', 'combined'):
        from srcs.utils.edge_store import load_store
        for atlas in ctx.atlases:
            store = load_store(ctx, '_'.join(atlas))
            print(f"{store.metadata['atlas']}: F {store.F_r.shape}; valid {store.valid_record.sum()}/{store.valid_record.size}")
        if stage == 'pairwise':
            from srcs.analysis.pairwise import selections, week_positions
            selections(ctx.weeks, ctx.section("pairwise"))
            week_positions(ctx.weeks, ctx.section("pairwise").get("week_positions"))
        if stage in ('stability', 'longitudinal'):
            from srcs.workflow.products import bootstrap_options, descriptive_products
            bootstrap_options(ctx, stage)
            for atlas in ctx.atlases:
                descriptive_products(ctx, '_'.join(atlas))
            if stage == 'longitudinal':
                from srcs.analysis.pairwise import week_positions
                week_positions(ctx.weeks, ctx.section("pairwise").get("week_positions"))
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
        from srcs.workflow.report import report_inputs, selected_sections
        inputs = report_inputs(ctx)
        section = getattr(ctx.args, 'report_section', 'all')
        paths = [inputs[name] for name in selected_sections(ctx, inputs)]
        if not any(path.is_file() for path in paths):
            print(f"No saved products for {section}; section will be marked skipped")


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-c', '--config', required=True)
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--stage')
    selection.add_argument('--start')
    parser.add_argument('--dry-run', action='store_true')
    parser.add_argument('--import-fc')
    parser.add_argument('--tag')
    parser.add_argument('--series-output', type=Path, required=True)
    parser.add_argument('--log-file', type=Path)
    parser.add_argument('--report-section', default='all')
    parser.add_argument('--function', action='append', default=[])
    args = parser.parse_args(argv)
    ctx = load_context(args.config, create_dirs=False, enabled_functions=args.function)
    ctx.args.tag = args.tag
    ctx.args.series_output = args.series_output.resolve()
    if args.dry_run:
        return execute(ctx,args,parser)
    from srcs.utils.logger import setup_logger
    if args.log_file is not None:
        setup_logger(str(args.log_file))
    from srcs.core.runtime import exclusive_run
    try:
        with exclusive_run(Path(ctx.output_path) / "metadata/workflow.lock"):
            return execute(ctx, args, parser)
    except Exception:
        ctx.logger.exception("Stage failed: %s", args.stage or args.start or "default")
        raise


def execute(ctx,args,parser):
    identity_path = Path(ctx.output_path) / "metadata/data_identity.json"
    if identity_path.exists() and json.loads(identity_path.read_text(encoding="utf-8")) != ctx.data_identity:
        raise ValueError("Data identity collision; use a distinct --tag")
    ctx.args.report_section = args.report_section
    workflow = ctx.workflow
    start = args.start or (workflow.source.spec["id"] if workflow.source else "fc")
    stages = select_steps(workflow.order, stage=args.stage, start=None if args.stage else start)
    workflow.validate_plan(stages)
    if args.import_fc and stages[0] != 'edges':
        parser.error('--import-fc requires --stage/--start edges')
    print(f'Output: {ctx.output_path}\nStages: {", ".join(stages)}', flush=True)
    if args.import_fc:
        if not args.dry_run and not identity_path.exists():
            save_json(identity_path, ctx.data_identity)
        details = import_fc(ctx, args.import_fc, dry_run=args.dry_run)
        print(f'FC import verified: {len(details["tables"])} tables; {details["audit_source"]}', flush=True)
    first = workflow.stages[stages[0]]
    if first.validate:
        getattr(importlib.import_module(first.module), first.validate)(ctx)
    else:
        validate(ctx, stages[0], importing=bool(args.import_fc))
    if args.dry_run:
        return 0
    if not identity_path.exists():
        save_json(identity_path, ctx.data_identity)
    if not (Path(ctx.output_path) / "metadata/analysis_config.json").exists():
        ctx.save_analysis_metadata()
    # Each stage validates when selected as a direct entry, or on demand in its adapter.
    for stage in stages:
        if stage == 'fc' and (Path(ctx.output_path) / 'metadata/fc_import_provenance.json').exists():
            raise ValueError('This output contains imported FC; start edges, or use a separate config to recompute FC')
        marker = Path(ctx.output_dir('Manifests')) / f'{stage}.json'
        state = dict(stage=stage, config_tag=ctx.analysis_tag, status='running',
                     started=datetime.now(timezone.utc).isoformat())
        save_json(marker, state)
        ctx.logger.info('Starting stage: %s', stage)
        try:
            action = importlib.import_module(workflow.stages[stage].module).run
            action(ctx)
        except BaseException as exc:
            state.update(status='interrupted' if isinstance(exc, KeyboardInterrupt) else 'failed', reason=str(exc))
            save_json(marker, state)
            raise
        state.update(status='complete', finished=datetime.now(timezone.utc).isoformat())
        save_json(marker, state)
        ctx.logger.info('Completed stage: %s', stage)
    return 0
