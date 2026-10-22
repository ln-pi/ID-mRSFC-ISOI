"""Identification of two aligned sets of FC edge vectors; no I/O or plotting."""
import numpy as np
from scipy.stats import pearsonr, ttest_rel, wilcoxon
from statsmodels.stats.diagnostic import lilliefors

from srcs.analysis.connectome import pattern_pearson
from srcs.analysis.statistics import individual_identifiability, paired_effect_summary, monte_carlo_pvalue_greater


def identification_accuracy(similarity):
    matrix = np.asarray(similarity, dtype=float)
    # Validate square shape, aligned cohort size and finite observations.
    individual_identifiability(matrix)
    own = np.diag(matrix)
    other = ~np.eye(len(matrix), dtype=bool)
    rows = np.all(np.where(other, own[:, None] > matrix, True), axis=1)
    columns = np.all(np.where(other, own[:, None] > matrix.T, True), axis=1)
    return float(np.mean(rows & columns)), float(rows.mean()), float(columns.mean())


def identify(first, second):
    """Average FC upstream, then Pearson and legacy similarity-z cap +/-4.

    Inputs are aligned mouse x edge arrays already in FC Fisher-z space.
    The valid even-length P1 convention is preserved; no FC transform is repeated.
    """
    first, second = np.asarray(first, float), np.asarray(second, float)
    if (first.shape != second.shape or first.ndim != 2 or len(first) < 2
            or first.shape[1] < 3):
        raise ValueError("Identification requires matching (mouse, edge) arrays with >=2 mice and >=3 edges")
    raw = pattern_pearson(first, second)
    if not np.isfinite(raw).all():
        raise ValueError("Identification patterns must be finite and nonconstant")
    similarity = np.clip(np.arctanh(np.clip(raw, -np.tanh(4), np.tanh(4))), -4, 4)
    metrics = individual_identifiability(similarity)
    individual, group = metrics["individual_similarity"], metrics["other_similarity"]
    relative = 1 - np.divide(group, individual, out=np.full_like(individual, np.nan), where=individual != 0)
    accuracy, row_accuracy, col_accuracy = identification_accuracy(similarity)
    li = lilliefors(individual) if len(individual) >= 4 and np.ptp(individual) > 0 else (np.nan, np.nan)
    lg = lilliefors(group) if len(group) >= 4 and np.ptp(group) > 0 else (np.nan, np.nan)
    t, tp = ttest_rel(individual, group)
    w, wp = wilcoxon(individual, group) if np.any(individual != group) else (0., 1.)
    r, rp = pearsonr(individual, group)
    effect = paired_effect_summary(individual, group)
    statistics = dict(identification_accuracy=accuracy, row_identification_accuracy=row_accuracy,
        col_identification_accuracy=col_accuracy, Indiv_lilliefors_statistic=li[0], Indiv_lilliefors_p_value=li[1],
        Group_lilliefors_statistic=lg[0], Group_lilliefors_p_value=lg[1], t_test_statistic=t, t_test_p_value=tp,
        wilcoxon_statistic=w, wilcoxon_p_value=wp, pearson_correlation=r, pearson_p_value=rp,
        **effect)
    statistics["effect_size_direction"] = "individual_minus_group"
    return dict(similarity=similarity, individual=individual, group=group,
                relative_effect=relative, cohen_d=effect["cohen_dz"], statistics=statistics)


def permutation_test(similarity, iterations=1000, seed=None):
    """Relabel second-half mice; keep every mouse's observed edge pattern intact."""
    if isinstance(iterations, bool) or not isinstance(iterations, (int, np.integer)) or iterations < 1:
        raise ValueError("Permutation iterations must be a positive integer")
    matrix = np.asarray(similarity, float)
    real = identification_accuracy(matrix)
    rng = np.random.default_rng(seed)
    samples = np.asarray([identification_accuracy(matrix[:, rng.permutation(len(matrix))])
                          for _ in range(iterations)]).T
    primary, rows, columns = [monte_carlo_pvalue_greater(value, sample) for value, sample in zip(real, samples)]
    return dict(real_accuracy=real, shuffle_acc=tuple(samples), shuffle_mean=samples[0].mean(),
        shuffle_std=samples[0].std(), **primary, seed=seed, row_permutation=rows, col_permutation=columns)
