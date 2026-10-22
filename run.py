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
from srcs.workflow import STAGES


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if argv and argv[0] == '--worker':
        os.chdir(ROOT)
        from srcs.utils.pipeline import main as worker
        return worker(argv[1:])
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    configs = parser.add_mutually_exclusive_group()
    configs.add_argument('-c', '--config')
    configs.add_argument('--configs', nargs='+')
    stages = parser.add_mutually_exclusive_group()
    stages.add_argument('--stage', choices=STAGES)
    stages.add_argument('--start', choices=STAGES)
    parser.add_argument('--dry-run', action='store_true', help='Validate selected inputs and print output locations without running stages')
    parser.add_argument('--import-fc', metavar='SOURCE_OUTPUT', help='Explicit source output/correlation directory; requires -c and edges stage')
    parser.add_argument('--report-section', choices=('all', 'pairwise', 'combined', 'statistics'), default='all',
                        help='Saved-result section to redraw during report (default: all available)')
    parser.add_argument('--report-root', default='reports/automatic', help='Batch plan/status directory')
    args = parser.parse_args(argv)
    start = args.stage or args.start or 'fc'
    if args.import_fc and (not args.config or start != 'edges'):
        parser.error('--import-fc requires -c and --start/--stage edges')
    configs = [args.config] if args.config else args.configs or ['configs/Aseries.yaml', 'configs/Nseries.yaml']
    # Absolute paths are accepted; relative CLI paths have the same project-root semantics as P2.
    for config in configs:
        config = (ROOT / config).resolve()
        if not config.is_file():
            parser.error(f'Configuration does not exist: {config}')
    from srcs.utils.series_workflow import run_batch
    code, _ = run_batch(configs, project_root=ROOT, stage=args.stage, start=args.start or 'fc',
        dry_run=args.dry_run, import_source=args.import_fc, report_section=args.report_section,
        report_root=args.report_root)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
