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
        from srcs.workflow.worker import main as worker
        return worker(argv[1:])
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    configs = parser.add_mutually_exclusive_group()
    configs.add_argument('-c', '--config')
    configs.add_argument('--configs', nargs='+')
    stages = parser.add_mutually_exclusive_group()
    stages.add_argument('--stage', choices=(*STAGES,"seedfc"))
    stages.add_argument('--start', choices=(*STAGES,"seedfc"))
    parser.add_argument('--dry-run', action='store_true', help='Validate selected inputs and print output locations without running stages')
    parser.add_argument('--import-fc', metavar='SOURCE_OUTPUT', help='Explicit source output/correlation directory; requires -c and edges stage')
    parser.add_argument('--from', dest='saved_output', metavar='DIRECTORY')
    parser.add_argument('--compare', nargs=2, metavar=('FIRST','SECOND'))
    parser.add_argument('--report-section', choices=('all', 'geometry', 'pairwise', 'descriptive', 'stability', 'longitudinal', 'combined', 'statistics'), default='all',
                        help='Saved-result section to redraw during report (default: all available)')
    parser.add_argument('--report-root', default='reports/automatic', help='Batch plan/status directory')
    parser.add_argument('--tag', help='Append a suffix to the data output directory')
    args = parser.parse_args(argv)
    if args.compare:
        if any((args.stage,args.start,args.saved_output,args.config,args.configs,args.tag,args.import_fc,args.dry_run)):
            parser.error('--compare is a standalone operation')
        from srcs.workflow.compare import compare
        roots = [ROOT / name if (ROOT / name).is_dir() else ROOT / "data/output" / name for name in args.compare]
        print(compare(*roots))
        return 0
    if args.saved_output:
        if args.stage != 'report' or args.config or args.configs or args.import_fc or args.tag:
            parser.error('--from accepts only saved-result report execution')
        from srcs.workflow.report import redraw_saved
        return redraw_saved((ROOT / args.saved_output).resolve(), args.report_section, args.dry_run)
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
        report_root=args.report_root, tag=args.tag)
    return code


if __name__ == '__main__':
    raise SystemExit(main())
