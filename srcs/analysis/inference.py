"""Existing within-mouse combined statistics; no new group-level inference."""
from itertools import combinations

import numpy as np
import pandas as pd
import pingouin as pg
from scipy.stats import ttest_rel
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multitest import multipletests

from srcs.analysis.statistics import paired_effect_summary


def observations(stat, atlas, metadata):
    axis = "week" if metadata["split_axis"] == "run" else "run"
    mice = metadata["mice"]
    rows = []
    for condition in stat.columns:
        cell = stat.at[atlas, condition]
        if len(cell[0]) != len(mice) or len(cell[1]) != len(mice):
            raise ValueError("Combined observations differ from the saved cohort")
        for i, mouse in enumerate(mice):
            for label, values in (("Group", cell[0]), ("Indiv", cell[1])):
                rows.append(dict(mouse=mouse, condition=label, **{axis: str(condition)}, similarity=values[i]))
    return pd.DataFrame(rows)


def repeated_measures(frame, within):
    fitted = AnovaRM(data=frame, depvar="similarity", subject="mouse", within=within).fit()
    corrected = pg.rm_anova(data=frame, dv="similarity", within=within, subject="mouse", detailed=True)
    return fitted.anova_table, corrected


def posthoc(stat, condition, exclude=()):
    atlases = [a for a in stat.index if a not in exclude]
    rows = []
    for metric, position in (("group", 0), ("indiv", 1), ("rel", 2)):
        for first, second in combinations(atlases, 2):
            a, b = np.asarray(stat.at[first, condition][position]), np.asarray(stat.at[second, condition][position])
            effect = paired_effect_summary(a, b)
            t, p = ttest_rel(a, b)
            rows.append(dict(metric=metric, comparison=f"{first} vs {second}", t=t, p_uncorrected=p, **effect))
    if rows:
        finite = np.isfinite([r["p_uncorrected"] for r in rows])
        adjusted = np.full(len(rows), np.nan)
        rejected = np.zeros(len(rows), bool)
        if finite.any():
            rejected[finite], adjusted[finite] = multipletests(
                np.asarray([r["p_uncorrected"] for r in rows])[finite], method="fdr_bh")[:2]
        for i, row in enumerate(rows):
            row.update(p_FDR=adjusted[i], significant=bool(rejected[i]),
                       test_status="complete" if finite[i] else "undefined_p_value")
    return pd.DataFrame(rows)
