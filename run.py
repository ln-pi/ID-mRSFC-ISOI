"""ISOI unified entry point (run from any working directory).

Stages:
  preprocess   raw recordings -> processed signals (manual registration may apply)
  fc           processed signals -> Correlation matrices and FC audit
  reliability  FC -> within-week Run reliability
  weekly       FC -> weekly representatives
  longitudinal FC -> between-week identification and edge contributions
  report       saved reliability/longitudinal results -> tables and figures
  phenotype    FC + supplied observations/events -> phenotype alignment

Examples:
  python run.py                           # Aseries and Nseries, fc onward
  python run.py -c configs/Aseries.yaml --stage preprocess
  python run.py -c configs/Aseries.yaml --stage fc
  python run.py --start reliability       # requires FC matrices AND FC audit
  python run.py --stage report            # requires saved numerical results
  python run.py --dry-run                 # inspect inputs and execution plan

--stage executes one stage; --start executes that stage and all later main stages.
After deleting everything except Correlation, restart with --start fc to recreate
FC audit records. Existing valid FC caches are reused under the usual checks.
Outputs: Correlation/ (unchanged), results/, report/, metadata/.
Terminal: live stage progress and final status; warnings remain in the summary.
Logs: each config's paths.log, all stages appended in order (Aseries/Nseries separate).
"""
import argparse
from pathlib import Path
import sys

from srcs.utils.series_workflow import run_batch, STAGES


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if argv and argv[0] == '--worker':
        from srcs.utils.pipeline import main as worker
        return worker(argv[1:])
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    configs = parser.add_mutually_exclusive_group()
    configs.add_argument('-c', '--config', help='One configuration')
    configs.add_argument('--configs', nargs='+', help='Configurations executed sequentially (default: Aseries and Nseries)')
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument('--stage', choices=STAGES, help='Execute only this stage')
    selection.add_argument('--start', choices=STAGES, help='Execute this and all following stages (default: fc)')
    parser.add_argument('--dry-run', action='store_true', help='Validate inputs and save a run plan')
    parser.add_argument('--report-root', default='reports/automatic')
    parser.add_argument('--plot-workers', type=int, default=4, help='Report rendering processes (default: 4; show=true uses one)')
    parser.add_argument('--import-fc', metavar='SOURCE_OUTPUT', help='Explicitly import FC before analysis; requires -c')
    args = parser.parse_args(argv)
    if args.plot_workers < 1:
        parser.error('--plot-workers must be positive')
    if args.import_fc and (not args.config or args.dry_run or (args.stage or args.start or 'fc') in ('fc', 'preprocess')):
        parser.error('--import-fc requires -c and an analysis --stage/--start; it cannot be combined with --dry-run')
    configs = [args.config] if args.config else args.configs or ['configs/Aseries.yaml', 'configs/Nseries.yaml']
    code, _ = run_batch(configs, project_root=Path(__file__).resolve().parent,
                        start=args.start or 'fc', stage=args.stage, dry_run=args.dry_run,
                        report_root=args.report_root, import_source=args.import_fc, plot_workers=args.plot_workers)
    return code


if __name__ == "__main__":
    raise SystemExit(main())
