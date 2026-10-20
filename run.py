"""Run P1 stages sequentially from PowerShell or any working directory.

python run.py --dry-run
python run.py --configs configs/Aseries.yaml configs/Nseries.yaml --start identification
python run.py -c configs/Aseries.yaml --start identification --import-fc SOURCE_OUTPUT
"""
import argparse
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
STAGES = ('preprocess', 'fc', 'identification', 'statistics')


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
    parser.add_argument('--import-fc', metavar='SOURCE_OUTPUT', help='P2 output root or Correlation folder; requires -c and --start/--stage identification')
    args = parser.parse_args(argv)
    start = args.stage or args.start or 'fc'
    if args.import_fc and (not args.config or start != 'identification'):
        parser.error('--import-fc requires -c and --start/--stage identification')
    configs = [args.config] if args.config else args.configs or ['configs/Aseries.yaml', 'configs/Nseries.yaml']
    # Absolute paths are accepted; relative CLI paths have the same project-root semantics as P2.
    for config in configs:
        config = (ROOT / config).resolve()
        if not config.is_file():
            parser.error(f'Configuration does not exist: {config}')
    for config in configs:
        cmd = [sys.executable, '-u', str(ROOT / 'run.py'), '--worker', '-c', str((ROOT / config).resolve()),
               '--stage' if args.stage else '--start', start]
        if args.dry_run:
            cmd.append('--dry-run')
        if args.import_fc:
            cmd.extend(['--import-fc', str((ROOT / args.import_fc).resolve())])
        print(f'\n=== {Path(config).name}: {start} ===', flush=True)
        result = subprocess.run(cmd, cwd=ROOT)
        if result.returncode:
            return result.returncode
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
