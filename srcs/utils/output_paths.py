"""Output directory layout; Correlation retains its established location."""
from pathlib import Path

DIRECTORIES = {
    'Correlation': 'correlation',
    'Reliability': 'results/reliability',
    'Weekly': 'results/weekly',
    'Identification': 'results/longitudinal',
    'Phenotype': 'results/phenotype',
    'Legacy': 'results/legacy',
    'Statistics': 'report',
    'Manifests': 'metadata/manifests',
    'BatchReports': 'metadata/batch_reports',
}


def fc_reports(root):
    """Read older FC audit locations until the user regenerates them."""
    root = Path(root)
    current = list((root / DIRECTORIES['BatchReports']).glob('02_correlation_*.json'))
    return current or list((root / 'BatchReports').glob('02_correlation_*.json'))
