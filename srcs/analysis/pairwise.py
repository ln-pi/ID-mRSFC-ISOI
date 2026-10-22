"""Phase selection and candidate-aware per-mouse metrics; no I/O or plotting."""
from itertools import combinations
import json
import re
import numpy as np
from srcs.analysis.statistics import individual_identifiability


def week_positions(ctx):
    supplied = ctx.section("pairwise").get("week_positions")
    if supplied is None:
        if not all(re.fullmatch(r"\d+w", w) for w in ctx.weeks):
            raise ValueError("Provide pairwise.week_positions for non-numeric week IDs")
        supplied = {w: int(w[:-1]) for w in ctx.weeks}
    if (set(supplied) != set(ctx.weeks) or any(isinstance(x, bool) or not isinstance(x, (int, float))
            or not np.isfinite(x) for x in supplied.values())
            or any(supplied[b] <= supplied[a] for a, b in zip(ctx.weeks, ctx.weeks[1:]))):
        raise ValueError("Week positions must be finite, strictly increasing, and cover all configured weeks")
    return supplied


def selections(ctx):
    """Phase labels select full arrays, never alter their construction."""
    cfg = ctx.section("pairwise")
    phases = cfg.get("phases", {})
    for name, ids in phases.items():
        if not isinstance(ids, list) or not ids or len(set(ids)) != len(ids) or not set(ids) <= set(ctx.weeks):
            raise ValueError(f"Invalid phase weeks: {name}")
    views = {"all_weeks": list(combinations(ctx.weeks, 2))}
    for name, ids in phases.items():
        ordered = [w for w in ctx.weeks if w in ids]
        views[f"within_{name}"] = list(combinations(ordered, 2))
    for request in cfg.get("phase_comparisons", []):
        if len(request) != 2 or any(p not in phases for p in request):
            raise ValueError("phase_comparisons must name two configured phases")
        a, b = request
        if a == b or set(phases[a]) & set(phases[b]):
            raise ValueError("Between-phase windows must be distinct and disjoint")
        views[f"{a}_to_{b}"] = [(u, v) for u in ctx.weeks if u in phases[a]
                               for v in ctx.weeks if v in phases[b]]
    return views


METRICS = ["individual_similarity", "other_similarity", "self_match_advantage",
           "self_match_percentile", "row_self_match_percentile", "col_self_match_percentile",
           "row_identified", "column_identified", "identified"]
COLUMNS = ["atlas", "run", "view", "cohort", "mouse", "week_u", "week_v", "week_gap",
           "n_candidates", "candidate_ids", "scale", *METRICS]


def panel_metrics(panel, candidates):
    selected = panel[np.ix_(candidates, candidates)]
    if len(candidates) < 2 or not np.isfinite(selected).all():
        return None
    result = individual_identifiability(selected)
    diag = np.diag(selected)
    others = ~np.eye(len(candidates), dtype=bool)
    result["row_identified"] = np.all(np.where(others, diag[:, None] > selected, True), axis=1)
    result["column_identified"] = np.all(np.where(others, diag[:, None] > selected.T, True), axis=1)
    result["identified"] = result["row_identified"] & result["column_identified"]
    return result


def analyze_views(store, run, arrays, views, positions):
    raw = arrays["S_r"]
    defined = arrays["pattern_reason"] == ""
    rows, audits = [], []
    for name, pairs in views.items():
        window = sorted({w for pair in pairs for w in pair}, key=store.weeks.index)
        fixed = [i for i in range(len(store.mice)) if window and all(defined[i, store.weeks.index(w)] for w in window)]
        audit = dict(view=name, run=run, atlas=store.metadata["atlas"],
                     pairs=[dict(week_u=u, week_v=v, week_gap=abs(positions[v]-positions[u])) for u, v in pairs],
                     fixed_mice=[store.mice[i] for i in fixed],
                     fixed_policy="all_pattern_valid_weeks_in_selected_window",
                     available_policy="common_pattern_valid_mice_at_both_endpoints_per_pair",
                     status="no_cross_week_pairs" if not pairs else "complete", panels=[])
        for u, v in pairs:
            a, b = store.weeks.index(u), store.weeks.index(v)
            available = np.flatnonzero(defined[:, a] & defined[:, b]).tolist()
            for policy, candidates in (("available_pair", available), ("fixed_window", fixed)):
                metrics = panel_metrics(raw[:, :, a, b], candidates)
                audit["panels"].append(dict(week_u=u, week_v=v, cohort=policy,
                    candidate_ids=[store.mice[i] for i in candidates], n_candidates=len(candidates),
                    status="complete" if metrics is not None else "fewer_than_two_candidates"))
                if metrics is None:
                    continue
                for k, i in enumerate(candidates):
                    rows.append(dict(atlas=store.metadata["atlas"], run=run, view=name, cohort=policy,
                        mouse=store.mice[i], week_u=u, week_v=v, week_gap=abs(positions[v]-positions[u]),
                        n_candidates=len(candidates), candidate_ids=json.dumps([store.mice[j] for j in candidates]),
                        scale="pattern_Pearson_r_of_F_z", **{metric: values[k] for metric, values in metrics.items()}))
        audits.append(audit)
    return rows, audits
