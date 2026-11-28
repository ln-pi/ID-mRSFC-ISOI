"""Run the ISOI workflow. Each new execution publishes one batch directory."""
import argparse
import os
from pathlib import Path
import sys


def main(argv=None):
    root = Path(__file__).resolve().parent
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == '--worker':
        os.chdir(root)
        from srcs.workflow.worker import main as worker
        return worker(argv[1:])
    parser = argparse.ArgumentParser(description=__doc__)
    configs = parser.add_mutually_exclusive_group()
    configs.add_argument('-c', '--config')
    configs.add_argument('--configs', nargs='+')
    stages = parser.add_mutually_exclusive_group()
    stages.add_argument('--stage', help='Run one configured stage')
    stages.add_argument('--start', help='Run from one configured stage')
    parser.add_argument('--dry-run', action='store_true', help='Print and validate a plan without creating output files')
    parser.add_argument('--import-fc', metavar='SOURCE_OUTPUT', help='Explicit source output/correlation directory; requires -c and edges stage')
    parser.add_argument('--from', dest='saved_output', metavar='DIRECTORY')
    parser.add_argument('--compare', nargs=2, metavar=('FIRST','SECOND'))
    parser.add_argument('--function', action='append', default=[], help='Enable a function plugin for this run')
    parser.add_argument('--visualization-config', help='Override visualization settings for a saved report or comparison')
    parser.add_argument('--report-section', default='all', help='Saved-result section to render (default: all available)')
    destinations = parser.add_mutually_exclusive_group()
    destinations.add_argument('--output-root', help='Parent directory for the new batch; defaults to paths.output')
    destinations.add_argument('--report-root', dest='output_root', help=argparse.SUPPRESS)
    destinations.add_argument('--resume-batch', metavar='BATCH_DIRECTORY', help='Reuse verified prerequisites in an existing batch')
    parser.add_argument('--tag', help='Append a suffix to the series output directory')
    args = parser.parse_args(argv)
    if args.compare:
        if any((args.stage, args.start, args.saved_output, args.config, args.configs, args.tag,
                args.import_fc, args.resume_batch, args.report_section != 'all')):
            parser.error('--compare is a standalone operation')
        from srcs.workflow.compare import compare_saved
        return compare_saved(*args.compare, project_root=root, enabled_functions=args.function,
                             visualization_config=args.visualization_config, output_root=args.output_root,
                             dry_run=args.dry_run)
    if args.saved_output:
        if args.stage != 'report' or args.config or args.configs or args.import_fc or args.tag or args.resume_batch:
            parser.error('--from accepts only saved-result report execution')
        from srcs.workflow.batch_layout import series_input
        from srcs.workflow.report import redraw_saved
        return redraw_saved(series_input(root, args.saved_output), args.report_section, args.dry_run,
                            args.visualization_config, args.function, output_root=args.output_root)
    if args.visualization_config:
        parser.error('--visualization-config requires --from or --compare')
    start = args.stage or args.start or 'fc'
    if args.import_fc and (not args.config or start != 'edges'):
        parser.error('--import-fc requires -c and --start/--stage edges')
    configs = [args.config] if args.config else args.configs or ['configs/Aseries.yaml', 'configs/Nseries.yaml']
    for config in configs:
        if not (root / config).resolve().is_file():
            parser.error(f'Configuration does not exist: {config}')
    from srcs.workflow.batch import run_batch
    code, _ = run_batch(configs, project_root=root, stage=args.stage, start=args.start,
        dry_run=args.dry_run, import_source=args.import_fc, report_section=args.report_section,
        output_root=args.output_root, resume_batch=args.resume_batch, tag=args.tag, enabled_functions=args.function)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
