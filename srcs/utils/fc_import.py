"""Explicit transfer of labeled FC tables, without requiring raw recordings.

Copy table bytes unchanged. Configuration hashes identify analysis settings, not
FC compatibility. Source preprocessing is recorded, never claimed equivalent.
"""
from collections import Counter
import hashlib
import json
from pathlib import Path

import pandas as pd

from srcs.utils.output_paths import fc_reports, correlation_folder
from srcs.utils.results import _atomic, save_json


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_import(ctx):
    path = Path(ctx.output_path) / 'metadata/fc_import_provenance.json'
    if path.exists():
        saved = json.loads(path.read_text(encoding='utf-8'))
        if saved['config_tag'] != ctx.analysis_tag:
            raise ValueError('Imported FC belongs to another analysis configuration')
        for row in saved['tables']:
            target = Path(ctx.correlation_path) / row['name']
            if not target.is_file() or digest(target) != row['sha256']:
                raise ValueError(f'Imported FC is missing or changed: {target}')


def import_fc(ctx, source, dry_run=False):
    """Accept a P2 output root or its Correlation folder (including an in-place copy).

    A supplied audit must have exact scope. Without it, missing cells are audited
    from the tables, while malformed matrices fail. Nothing is written until all
    tables and existing destinations have passed validation.
    """
    from srcs.analysis.connectivity import valid_fc
    source = Path(source).resolve()
    folder = source if source.name.casefold() == 'correlation' else correlation_folder(source)
    root = folder.parent
    tags = ['_'.join(a) for a in ctx.atlases]
    expected = {(w, m, r, a) for w in ctx.weeks for m in ctx.mice for r in ctx.runs for a in tags}
    reports = fc_reports(root)
    audit = None
    provenance = []
    if reports:
        report = max(reports, key=lambda p: (p.stat().st_mtime_ns, p.name))
        rows = json.loads(report.read_text(encoding='utf-8'))['records']
        audit = {tuple(row[k] for k in ('week', 'mouse', 'run', 'atlas')): row for row in rows}
        if len(audit) != len(rows) or set(audit) != expected:
            raise ValueError('Source FC audit scope differs from configured weeks/mice/runs/atlases')
        provenance.append(dict(path=str(report), sha256=digest(report)))
    config_path = root / 'metadata/analysis_config.json'
    if not config_path.exists():
        config_path = root / 'analysis_config.json'
    source_metadata = json.loads(config_path.read_text(encoding='utf-8')) if config_path.exists() else None
    if source_metadata is not None:
        variables = source_metadata['config']['variables']
        for key in ('dataset', 'Week_ids', 'Mouse_ids', 'Run_ids', 'Atlas'):
            if variables.get(key) != ctx.variables.get(key):
                raise ValueError(f'Source configuration scope mismatch: {key}')
        provenance.append(dict(path=str(config_path), sha256=digest(config_path)))
    prepared, records = [], []
    for atlas in tags:
        labels = None
        for week in ctx.weeks:
            name = f'{week}_correlation_matrices_{atlas}_excluded.pkl'
            path = folder / name
            frame = pd.read_pickle(path)
            if not isinstance(frame, pd.DataFrame) or list(frame.index) != ctx.mice or list(frame.columns) != ctx.runs:
                raise ValueError(f'FC table dimensions/order differ from configuration: {path}')
            for mouse in ctx.mice:
                for run in ctx.runs:
                    cell = frame.at[mouse, run]
                    missing = cell is None or (pd.api.types.is_scalar(cell) and bool(pd.isna(cell)))
                    if not missing:
                        if not isinstance(cell, pd.DataFrame) or not cell.index.is_unique or len(cell) < 3 or not valid_fc(cell, list(cell.index)):
                            raise ValueError(f'Invalid FC matrix: {path}, {mouse}, {run}')
                        if labels is not None and labels != list(cell.index):
                            raise ValueError(f'Inconsistent ROI labels/order: {path}, {mouse}, {run}')
                        labels = list(cell.index)
                    original = audit[(week, mouse, run, atlas)] if audit is not None else None
                    if original and original['status'] in ('success', 'reused') and missing:
                        raise ValueError(f'Source audit reports success for missing FC: {path}, {mouse}, {run}')
                    status = original['status'] if original else ('missing' if missing else 'reused')
                    records.append(dict(week=week, mouse=mouse, run=run, atlas=atlas,
                        status=status, reason=original.get('reason', '') if original else ('missing FC cell' if missing else ''),
                        input_path=str(path)))
            target = Path(ctx.correlation_path) / name
            sha = digest(path)
            if target.exists() and digest(target) != sha:
                raise FileExistsError(f'Target FC differs; use a separate output configuration: {target}')
            prepared.append((path, target, sha))
    details = dict(source=str(folder), config_tag=ctx.analysis_tag,
        verification='table dimensions, finite symmetric FC, ROI order, available source audit; preprocessing equivalence not assumed',
        audit_source='source_report' if audit is not None else 'reconstructed_from_cells',
        source_metadata=source_metadata, source_files=provenance,
        tables=[dict(name=p.name, sha256=sha) for p, _, sha in prepared])
    if not dry_run:
        for path, target, _ in prepared:
            if not target.exists():
                data = path.read_bytes()
                _atomic(target, lambda stream, data=data: stream.write(data), binary=True)
        save_json(Path(ctx.output_dir('BatchReports')) / '02_correlation_imported.json',
                  dict(counts=dict(Counter(row['status'] for row in records)), records=records))
        save_json(Path(ctx.output_dir('Metadata')) / 'fc_import_provenance.json', details)
    return details
