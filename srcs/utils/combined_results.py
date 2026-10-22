"""Validation and export of saved per-mouse combined results."""
import numpy as np
import pickle
from pathlib import Path
import pandas as pd

from srcs.analysis.statistics import individual_identifiability
from srcs.utils.results import save_table, save_json




def export_subject_metrics(ctx, stat, analysis_meta):
    """Use full similarity matrices: the aggregate alone cannot recover ranks."""
    rows = []
    excluded = set(ctx.section('analysis').get('statistics_exclude', []))
    atlases = ['_'.join(a) for a in ctx.atlases if '_'.join(a) not in excluded]
    mice = list(analysis_meta['mice'])
    if len(mice) != len(set(mice)):
        raise ValueError('Duplicate mouse IDs in analysis metadata.')
    for atlas in atlases:
        for condition in stat.columns:
            path = Path(ctx.output_dir('Identification', str(condition), 'Numeric', create=False))/f'identification_{atlas}.pkl'
            with path.open('rb') as stream:
                payload = pickle.load(stream)
            if payload.get('analysis') != analysis_meta or payload.get('mice') != mice:
                raise ValueError(f'Similarity metadata/cohort differs from aggregate: {path}')
            metrics = individual_identifiability(payload['similarity'])
            if len(metrics['individual_similarity']) != len(mice):
                raise ValueError(f'Similarity size differs from mouse labels: {path}')
            cell = stat.at[atlas, condition]
            for actual, expected in [(metrics['individual_similarity'], cell[1]),
                                     (metrics['other_similarity'], cell[0])]:
                if np.shape(actual) != np.shape(expected) or not np.allclose(actual, expected, rtol=1e-9, atol=1e-10):
                    raise ValueError(f'Similarity values differ from aggregate: {path}')
            for i, mouse in enumerate(mice):
                rows.append(dict(mouse=mouse, atlas=atlas, analysis_condition=str(condition),
                                 n_mice=len(mice), **{k:float(v[i]) for k,v in metrics.items()}))
    columns = ['mouse','atlas','analysis_condition','n_mice','individual_similarity','other_similarity',
               'self_match_advantage','self_match_percentile','row_self_match_percentile','col_self_match_percentile']
    frame = pd.DataFrame(rows, columns=columns)
    save_table(Path(ctx.output_dir('results','statistics'))/'individual_identifiability.csv',frame,index=False)
    save_json(Path(ctx.output_dir('Metadata'))/'individual_identifiability_metadata.json',dict(
        analysis=analysis_meta, status='complete' if rows else 'empty_selection',
        inference='descriptive_only', advantage='self minus mean of row/column nonself similarities',
        percentile='mean row/column fraction of nonself candidates beaten; exact tie=0.5',
        percentile_range=[0,1], similarity_scale='saved_combined_similarity_z_capped_4',
        caveats=['Shared nonself references induce dependence.',
                 'Candidate count/composition affects interpretation; not a disease severity score.']))
    return frame
