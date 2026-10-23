"""Run repeatability, residual persistence, edgewise variance components and ICC."""
from pathlib import Path

import numpy as np

from srcs.analysis.connectome import save_arrays
from srcs.analysis.stability import SCHEMA, residual_persistence, run_repeatability
from srcs.analysis.topology import edge_layout
from srcs.analysis.variance import balanced_reml, variance_intervals
from srcs.utils.analysis_products import bootstrap_options, descriptive_products, publish
from srcs.utils.edge_store import load_store, store_path
from srcs.utils.provenance import code_identity, sha256


def run(ctx):
    folder = Path(ctx.output_dir("results", "stability"))
    count, seed = bootstrap_options(ctx, "stability")
    tables = {key: [] for key in ("repeatability", "residual_pairs", "edge_icc", "regional_icc", "cohorts")}
    products, sources, audits = [], [], []
    code = code_identity()
    for atlas in ["_".join(a) for a in ctx.atlases]:
        store = load_store(ctx, atlas)
        descriptive, source = descriptive_products(ctx, atlas)
        sources.append(source)
        meta = descriptive["available_week"][0]
        weeks = meta["week_ids"]
        rows, audit = run_repeatability(store, weeks)
        tables["repeatability"].extend(dict(atlas=atlas, **row) for row in rows)
        audits.append(dict(atlas=atlas, same_day_templates=audit))
        for cohort, (m, arrays) in descriptive.items():
            persistence = residual_persistence(arrays["residual_group_loo"])
            for i, mouse in enumerate(m["mouse_ids"]):
                for r, view in enumerate(m["view_ids"]):
                    for u in range(len(weeks)):
                        for v in range(u+1, len(weeks)):
                            tables["residual_pairs"].append(dict(atlas=atlas, cohort=cohort, mouse=mouse,
                                view=view, week_u=weeks[u], week_v=weeks[v], pattern_r=persistence[i, r, u, v],
                                template_policy="fixed_window_members" if cohort == "fixed_window" else "dynamic_week_members"))
            path = folder/atlas/f"residual_{cohort}.npz"
            save_arrays(path, dict(schema=SCHEMA, atlas=atlas, cohort=cohort, mouse_ids=m["mouse_ids"],
                week_ids=weeks, view_ids=m["view_ids"], axes="mouse,view,week_u,week_v",
                same_record="diagonal_excluded", source=source, code=code),
                residual_similarity=persistence, valid_record=arrays["valid_record"])
            products.append(dict(kind="residual", atlas=atlas, cohort=cohort, path=str(path), sha256=sha256(path)))

        wi = [store.weeks.index(w) for w in weeks]
        complete = store.valid_record[:, wi, :].all(axis=(1, 2))
        ids = np.flatnonzero(complete)
        for i, mouse in enumerate(store.mice):
            missing = [f"{week}/{run}: {store.invalid_reason[i, u, r]}"
                       for week, u in zip(weeks, wi) for r, run in enumerate(store.runs)
                       if not store.valid_record[i, u, r]]
            tables["cohorts"].append(dict(atlas=atlas, mouse=mouse, included=bool(complete[i]),
                                          reason="; ".join(missing)))
        y = store.F_z[ids][:, wi]
        e = store.F_z.shape[-1]
        supported = len(ids) >= 3 and len(weeks) >= 2 and len(store.runs) >= 2
        model_status = "balanced_REML" if supported else "unsupported_need_3_complete_mice_2_weeks_2_runs"
        if supported:
            fit = balanced_reml(y)
            fit.update(variance_intervals(y, count, seed))
            fit["conditional_residual"] = (y-fit["week_mean"][None, :, None, :]
                -fit["run_effect"][None, None, :, :]-fit["stable_mouse_blup"][:, None, None, :]
                -fit["mouse_week_blup"][:, :, None, :])
        else:
            fit = {key: np.full(e, np.nan) for key in ("mouse_variance", "mouse_week_variance", "error_variance",
                                                       "stable_icc", "same_week_icc", "mean_run_stable_icc")}
            fit["status"] = np.full(e, "unsupported")
            fit["stable_mouse_blup"] = np.full((len(ids), e), np.nan)
        groups, edges, mapping = edge_layout(store.roi_labels, store.edge_index,
                                              ctx.section("descriptive").get("network_mapping"))
        for j, edge in enumerate(edges):
            row = dict(atlas=atlas, **edge, n_mice=len(ids), n_weeks=len(weeks), n_runs=len(store.runs),
                       model_status=model_status, status=fit["status"][j])
            for key in ("mouse_variance", "mouse_week_variance", "error_variance", "stable_icc",
                        "same_week_icc", "mean_run_stable_icc"):
                row[key] = fit[key][j]
                if key+"_ci" in fit:
                    row[key+"_ci_low"], row[key+"_ci_high"] = fit[key+"_ci"][:, j]
            tables["edge_icc"].append(row)
        for group in groups:
            indices = group["indices"]
            total = sum(fit[key][indices].sum() for key in ("mouse_variance", "mouse_week_variance", "error_variance"))
            values = fit["stable_icc"][indices]
            finite = values[np.isfinite(values)]
            tables["regional_icc"].append(dict(atlas=atlas, level=group["level"], region=group["name"],
                n_edges=len(indices), mean_edge_stable_icc=float(finite.mean()) if len(finite) else np.nan,
                ratio_of_summed_variances=float(fit["mouse_variance"][indices].sum()/total) if total > 0 else np.nan,
                mean_mouse_variance=float(fit["mouse_variance"][indices].mean())))
        path = folder/atlas/"variance_components.npz"
        save_arrays(path, dict(schema=SCHEMA, atlas=atlas, model_status=model_status,
            mouse_ids=[store.mice[i] for i in ids], excluded_mice=[m for i, m in enumerate(store.mice) if not complete[i]],
            week_ids=weeks, run_ids=store.runs, roi_labels=store.roi_labels, **mapping,
            scale="FC_edge_Fisher_z", model="week_fixed + run_fixed + mouse_random + mouse_week_random",
            estimator="balanced_orthogonal_constrained_REML_ordered_variance_pooling",
            bootstrap=dict(unit="whole_mouse", count=count, seed=seed, interval="percentile_95"),
            edge_store=dict(path=str(store_path(ctx, atlas)), sha256=sha256(store_path(ctx, atlas))), code=code),
            edge_index=store.edge_index, **fit)
        products.append(dict(kind="variance", atlas=atlas, path=str(path), sha256=sha256(path)))
        ctx.logger.info("Stability %s: ICC model %s, %s complete mice", atlas, model_status, len(ids))
    publish(folder, SCHEMA, products, tables, dict(schema=SCHEMA, code=code, sources=sources,
        formulas=dict(stable_icc="mouse_variance/(mouse_variance+mouse_week_variance+error_variance)",
            residual_similarity="Pearson_edges(DeltaFC_i,u,DeltaFC_i,v)",
            same_week_icc="(mouse_variance+mouse_week_variance)/total"),
        icc_target="single_Run_conditional_on_fixed_week_and_Run_effects_not_classical_unadjusted_ICC",
        stable_component="conditional_BLUP_of_mouse_random_intercept; not_pure_neural_or_disease_effect",
        bootstrap=dict(unit="whole_mouse", count=count, seed=seed, interval="percentile_95_not_boundary_test"),
        same_day_templates=audits, residual_templates="dynamic_available_and_fixed_window_kept_separate",
        interpretation="same_day_repeatability_not_long_term_stability; no_imputation; informative_missingness_not_resolved"))
