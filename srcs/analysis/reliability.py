"""Within-week repeatability, retaining every valid paired observation."""
from itertools import combinations
import numpy as np
import pandas as pd
from srcs.analysis.identification import edge_vector, pattern_similarity, identify
from srcs.utils.progress import progress_items


def analyze_reliability(ctx, records):
    rows, comparisons = [], {}
    for atlas in ['_'.join(a) for a in ctx.atlases]:
        for week in progress_items(ctx.weeks, phase=f'{atlas}: paired Run weeks'):
            for r1, r2 in combinations(ctx.runs, 2):
                mice = [m for m in ctx.mice if all((atlas, week, m, r) in records for r in (r1,r2))]
                vectors = []
                for m in mice:
                    a = edge_vector(records[atlas,week,m,r1])
                    b = edge_vector(records[atlas,week,m,r2])
                    try:
                        r = pattern_similarity(a,b)
                    except ValueError:
                        r = np.nan
                    var = np.var(a)+np.var(b)+(a.mean()-b.mean())**2
                    rows.append(dict(atlas=atlas, week=week, mouse=m, first_run=r1, second_run=r2,
                                     pattern_r=r, rmse_z=float(np.sqrt(np.mean((a-b)**2))),
                                     concordance_z=float(2*np.mean((a-a.mean())*(b-b.mean()))/var) if var else np.nan,
                                     n_edges=len(a), status='valid' if np.isfinite(r) else 'undefined_pattern',
                                     self_z=np.nan, other_z=np.nan, advantage_z=np.nan, percentile=np.nan))
                    vectors.append((a,b))
                key = f'{atlas}/{week}/{r1}-{r2}'
                if len(mice) >= 2:
                    try:
                        result = identify([v[0] for v in vectors], [v[1] for v in vectors], mice)
                    except ValueError as exc:
                        comparisons[key] = dict(status='unavailable', reason=str(exc), mice=mice)
                    else:
                        comparisons[key] = dict(status='complete', **result)
                        for row, i in zip(rows[-len(mice):], range(len(mice))):
                            row.update(self_z=result['self_z'][i], other_z=result['other_z'][i],
                                       advantage_z=result['advantage_z'][i])
                            row['percentile']=result['individual_metrics']['self_match_percentile'][i]
                else:
                    comparisons[key] = dict(status='unavailable', reason='fewer than two paired mice', mice=mice)
    return pd.DataFrame(rows), comparisons
