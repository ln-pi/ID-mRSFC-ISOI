"""Persist group/individual trajectories, dispersion trends and network slopes."""
from pathlib import Path

import numpy as np
import pandas as pd

from srcs.analysis.connectome import save_arrays
from srcs.analysis.longitudinal import SCHEMA, dispersion_trajectory, linear_slope, random_slopes, rank_persistence
from srcs.analysis.magnitude import select_views
from srcs.analysis.pairwise import week_positions
from srcs.analysis.topology import edge_layout
from srcs.utils.analysis_products import bootstrap_options, descriptive_products, publish
from srcs.utils.edge_store import load_store
from srcs.utils.provenance import code_identity, sha256
from srcs.utils.results import save_json


def run(ctx):
    folder = Path(ctx.output_dir("results", "longitudinal"))
    count, seed = bootstrap_options(ctx, "longitudinal")
    cfg = ctx.section("longitudinal")
    model_pairs = cfg.get("model_network_pairs", "within")
    if model_pairs not in ("within", "all"):
        raise ValueError("longitudinal.model_network_pairs must be within or all")
    maxiter = cfg.get("maxiter", 300)
    if isinstance(maxiter, bool) or not isinstance(maxiter, int) or maxiter < 1:
        raise ValueError("longitudinal.maxiter must be a positive integer")
    tables = {key: [] for key in ("dispersion", "trends", "subjects", "rank_pairs", "random_slopes", "mouse_slopes", "model_diagnostics")}
    products, sources = [], []
    code = code_identity()
    positions = week_positions(ctx)
    for atlas in ["_".join(a) for a in ctx.atlases]:
        store = load_store(ctx, atlas)
        descriptive, source = descriptive_products(ctx, atlas)
        sources.append(source)
        meta = descriptive["available_week"][0]
        weeks, views, values, valid, definitions = select_views(store, dict(ctx.section("descriptive"), weeks=meta["week_ids"]))
        time = np.array([positions[w] for w in weeks])
        if views != meta["view_ids"]:
            raise ValueError("Descriptive view definitions changed; rerun descriptive")
        for cohort, (m, arrays) in descriptive.items():
            mask = arrays["valid_record"]
            scores = np.sqrt(np.mean(arrays["residual_group_loo"]**2, axis=-1))
            observed, intervals = [], []
            for r, view in enumerate(views):
                result = dispersion_trajectory(values[:, :, r], mask[:, :, r], time, count, seed)
                observed.append(result["dispersion"])
                intervals.append(result["dispersion_ci"])
                tables["trends"].append(dict(atlas=atlas, cohort=cohort, view=view,
                    slope=result["slope"], slope_ci_low=result["slope_ci"][0], slope_ci_high=result["slope_ci"][1],
                    direction=result["direction"], bootstrap_valid=result["bootstrap_valid"],
                    interpretation="fixed_members_dispersion_trend" if cohort == "fixed_window" else "composition_sensitive_descriptive_trend"))
                for u, week in enumerate(weeks):
                    tables["dispersion"].append(dict(atlas=atlas, cohort=cohort, view=view, week=week,
                        time=time[u], dispersion=result["dispersion"][u],
                        ci_low=result["dispersion_ci"][0, u], ci_high=result["dispersion_ci"][1, u],
                        n_mice=result["n_valid_mice"][u], group_mean_fc=float(arrays["group_mean"][u, r].mean())))
                for i, mouse in enumerate(store.mice):
                    slope = linear_slope(time, scores[i, :, r])
                    for u, week in enumerate(weeks):
                        tables["subjects"].append(dict(atlas=atlas, cohort=cohort, view=view, mouse=mouse,
                            week=week, time=time[u], deviation_rms=scores[i, u, r],
                            mean_signed_deviation=float(arrays["residual_group_loo"][i, u, r].mean()),
                            deviation_linear_slope=slope, valid=bool(mask[i, u, r])))
                for row in rank_persistence(scores[:, :, r], mask[:, :, r]):
                    u, v = row.pop("u"), row.pop("v")
                    tables["rank_pairs"].append(dict(atlas=atlas, cohort=cohort, view=view,
                        week_u=weeks[u], week_v=weeks[v], time_gap=time[v]-time[u], **row))
            path = folder/atlas/f"trajectories_{cohort}.npz"
            save_arrays(path, dict(schema=SCHEMA, atlas=atlas, cohort=cohort, mouse_ids=store.mice,
                week_ids=weeks, view_ids=views, roi_labels=store.roi_labels, time_positions=time.tolist(),
                run_views=definitions, source=source, code=code),
                dispersion=np.stack(observed, axis=1), dispersion_ci=np.stack(intervals, axis=2),
                deviation_rms=scores, group_mean=arrays["group_mean"],
                group_step_change=np.diff(arrays["group_mean"], axis=0),
                valid_record=mask, edge_index=store.edge_index)
            products.append(dict(atlas=atlas, cohort=cohort, path=str(path), sha256=sha256(path)))

        groups, _, mapping = edge_layout(store.roi_labels, store.edge_index,
                                          ctx.section("descriptive").get("network_mapping"))
        targets = [("all_edges", np.arange(store.F_z.shape[-1]))]
        for group in groups:
            if group["level"] == "network_pair":
                a, b = group["name"].split(" | ")
                if model_pairs == "all" or a == b:
                    targets.append((group["name"], group["indices"]))
        wi = [store.weeks.index(w) for w in weeks]
        for target, edges in targets:
            observations = []
            for i, mouse in enumerate(store.mice):
                for u, week in enumerate(weeks):
                    for r, run_id in enumerate(store.runs):
                        if store.valid_record[i, wi[u], r]:
                            observations.append(dict(mouse=mouse, week=week, run=run_id, time=time[u],
                                value=float(store.F_z[i, wi[u], r, edges].mean())))
            frame = pd.DataFrame(observations, columns=["mouse", "week", "run", "time", "value"])
            summary, effects, diagnostics = random_slopes(frame, maxiter)
            base = dict(atlas=atlas, target=target, n_edges=len(edges), network_kind=mapping["network_kind"])
            tables["random_slopes"].append(dict(**base, **summary))
            tables["mouse_slopes"].extend(dict(**base, **row) for row in effects)
            tables["model_diagnostics"].extend(dict(**base, **row) for row in diagnostics)
            ctx.logger.info("Random slope %s/%s: %s", atlas, target, summary["status"])
    model_path = folder/"models.json"
    save_json(model_path, tables["random_slopes"])
    products.append(dict(kind="model_diagnostics", path=str(model_path), sha256=sha256(model_path)))
    publish(folder, SCHEMA, products, tables, dict(schema=SCHEMA, code=code, sources=sources,
        time_positions=positions, time_units="configured_week_positions_not_assumed_actual_acquisition_dates",
        trajectory_bootstrap=dict(unit="whole_mouse_trajectory", n_bootstrap=count, seed=seed, interval="percentile_95"),
        trend="OLS_slope_of_dispersion_over_time; increasing/decreasing_only_if_bootstrap_CI_excludes_zero",
        trend_scope="descriptive_convergence_divergence_not_disease_causality; compare_fixed_and_available_cohorts",
        rank="Spearman_across_common_mice_no_cell_level_p_values",
        random_slopes=dict(response="mean_raw_Fisher_z_FC_of_network_edges_not_dependent_LOO_scores",
            fixed="linear_time_and_Run", random="correlated_mouse_intercept_slope_plus_mouse_week_intercept",
            estimator="statsmodels_MixedLM_REML", model_network_pairs=model_pairs, maxiter=maxiter,
            missing="available_observations_no_imputation_MAR_assumption_not_correction_for_death",
            inference="no_automatic_significance_or_causal_claims; nonconverged_BLUPs_not_published"),
        reference="https://www.statsmodels.org/stable/generated/statsmodels.regression.mixed_linear_model.MixedLM.html"))
