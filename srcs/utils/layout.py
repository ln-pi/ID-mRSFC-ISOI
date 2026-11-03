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


def correlation_folder(root):
    """Write lowercase; read the old spelling if that is the only directory.

    Refuse two distinct directories on case-sensitive filesystems.
    """
    root = Path(root)
    lower, old = root / 'correlation', root / 'Correlation'
    if lower.exists() and old.exists() and not lower.samefile(old):
        raise ValueError(f'Ambiguous correlation/Correlation directories under {root}')
    return old if old.exists() and not lower.exists() else lower
