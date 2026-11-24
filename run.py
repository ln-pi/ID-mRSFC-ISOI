"""Run the unified ISOI workflow from PowerShell or any working directory.

python run.py --dry-run
python run.py --configs configs/Aseries.yaml configs/Nseries.yaml --start edges
python run.py -c configs/Aseries.yaml --start edges --import-fc SOURCE_OUTPUT
"""
import argparse
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == '--worker':
        os.chdir(ROOT)
        from srcs.workflow.worker import main as worker
        return worker(argv[1:])
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    configs = parser.add_mutually_exclusive_group()
    configs.add_argument('-c', '--config')
    configs.add_argument('--configs', nargs='+')
    stages = parser.add_mutually_exclusive_group()
    stages.add_argument('--stage', help='Run one configured stage')
    stages.add_argument('--start', help='Run from one configured stage')
    parser.add_argument('--dry-run', action='store_true', help='Validate selected inputs and print output locations without running stages')
    parser.add_argument('--import-fc', metavar='SOURCE_OUTPUT', help='Explicit source output/correlation directory; requires -c and edges stage')
    parser.add_argument('--from', dest='saved_output', metavar='DIRECTORY')
    parser.add_argument('--compare', nargs=2, metavar=('FIRST','SECOND'))
    parser.add_argument('--function', action='append', default=[], help='Enable a function plugin for this run')
    parser.add_argument('--visualization-config', help='Override visualization settings for a saved report or comparison')
    parser.add_argument('--report-section', default='all',
                        help='Saved-result section to redraw during report (default: all available)')
    parser.add_argument('--report-root', help='Batch plan/status directory (defaults to _batches under the first output root)')
    parser.add_argument('--tag', help='Append a suffix to the data output directory')
    args = parser.parse_args(argv)
    if args.compare:
        if any((args.stage,args.start,args.saved_output,args.config,args.configs,args.tag,args.import_fc,args.dry_run)):
            parser.error('--compare is a standalone operation')
        from srcs.workflow.compare import compare
        roots = [ROOT / name if (ROOT / name).is_dir() else ROOT / "data/output" / name for name in args.compare]
        from srcs.utils.visualization import saved_config, merge
        from srcs.workflow.plugins import enable_functions
        configs = [enable_functions(saved_config(root), args.function) for root in roots]
        if args.visualization_config:
            from srcs.utils.config import load_config
            override = load_config(ROOT / args.visualization_config)
            configs = [merge(config, {"visualization": override.get("visualization", override)}) for config in configs]
        function_results = {}
        print(compare(*roots, configs=configs, function_results=function_results))
        for name, path in function_results.items():
            print(f"{name}: {path}")
        return 0
    if args.saved_output:
        if args.stage != 'report' or args.config or args.configs or args.import_fc or args.tag:
            parser.error('--from accepts only saved-result report execution')
        from srcs.workflow.report import redraw_saved
        return redraw_saved((ROOT / args.saved_output).resolve(), args.report_section, args.dry_run,
                            args.visualization_config, args.function)
    if args.visualization_config:
        parser.error('--visualization-config requires --from or --compare')
    start = args.stage or args.start or 'fc'
    if args.import_fc and (not args.config or start != 'edges'):
        parser.error('--import-fc requires -c and --start/--stage edges')
    configs = [args.config] if args.config else args.configs or ['configs/Aseries.yaml', 'configs/Nseries.yaml']
    # Absolute paths are accepted; relative CLI paths have the same project-root semantics as P2.
    for config in configs:
        config = (ROOT / config).resolve()
        if not config.is_file():
            parser.error(f'Configuration does not exist: {config}')
    from srcs.workflow.batch import run_batch
    code, _ = run_batch(configs, project_root=ROOT, stage=args.stage, start=args.start,
        dry_run=args.dry_run, import_source=args.import_fc, report_section=args.report_section,
        report_root=args.report_root, tag=args.tag, enabled_functions=args.function)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
