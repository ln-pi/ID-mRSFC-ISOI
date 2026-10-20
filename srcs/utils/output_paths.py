"""P2-compatible top-level layout; combined and longitudinal results stay distinct."""
from pathlib import Path

DIRECTORIES = {
    'Correlation': 'correlation',
    'Identification': 'results/combined',
    'IdentificationFigures': 'report/combined',
    'Statistics': 'report',
    'BatchReports': 'metadata/batch_reports',
    'Manifests': 'metadata/manifests',
    'Metadata': 'metadata',
    'LegacyStructure': 'report/legacy_structure',
}


def output_parts(parts):
    if not parts:
        return parts
    first, *rest = parts
    if first == 'Identification' and len(rest) > 1:
        if rest[1] in ('Similarity', 'Stat'):
            first = 'IdentificationFigures'
        elif rest[1] == 'Numeric':
            rest[1] = 'numeric'
    if first == 'Statistics' and rest and rest[0] == 'Tables':
        rest[0] = 'tables'
    return (DIRECTORIES.get(first, first), *rest)


def fc_reports(root):
    root = Path(root)
    current = list((root / 'metadata/batch_reports').glob('02_correlation_*.json'))
    return current or list((root / 'BatchReports').glob('02_correlation_*.json'))
