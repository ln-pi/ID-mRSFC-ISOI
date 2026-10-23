"""Persist amplitude, LOO topology and transparent descriptive examples."""
from pathlib import Path

import numpy as np
import pandas as pd

from srcs.analysis.connectome import save_arrays
from srcs.analysis.descriptive import same_day_reliability
from srcs.analysis.magnitude import SCHEMA, magnitude, select_views
from srcs.analysis.topology import edge_layout, regional_summaries
from srcs.utils.results import save_json, save_table
from srcs.utils.edge_store import load_store, store_path
from srcs.utils.provenance import code_identity, sha256


def run(ctx):
    folder = Path(ctx.output_dir("results", "descriptive", "v2"))
    cfg = ctx.section("descriptive")
    reference = cfg.get("reference_week", "03w")
    products, inputs, selection = [], [], []
    tables = {name: [] for name in ("distances", "subjects", "dispersion", "regional", "examples", "reliability")}
    code = code_identity()
    for atlas in ["_".join(a) for a in ctx.atlases]:
        store = load_store(ctx, atlas)
        weeks, views, values, valid, definitions = select_views(store, cfg)
        groups, edges, mapping = edge_layout(store.roi_labels, store.edge_index, cfg.get("network_mapping"))
        save_table(folder / f"edges_{atlas}.csv", pd.DataFrame(edges), index=False)
        inputs.append(dict(path=str(store_path(ctx, atlas)), sha256=sha256(store_path(ctx, atlas))))
        tables["reliability"].extend(row for row in same_day_reliability(store) if row["week"] in weeks)
        all_weeks, _, all_values, all_valid, _ = select_views(store, dict(cfg, weeks=store.weeks))
        ri = all_weeks.index(reference) if reference in all_weeks else None
        for cohort in ("available_week", "fixed_window"):
            selected = valid.copy()
            fixed = valid.all(axis=1)
            if cohort == "fixed_window":
                selected &= fixed[:, None, :]
            n, w, v, e = values.shape
            residual, own_change = np.full_like(values, np.nan), np.full_like(values, np.nan)
            variance = np.full((w, v, e), np.nan)
            group_mean = np.full_like(variance, np.nan)
            distances = {key: np.full((w, v, n, n), np.nan)
                         for key in ("euclidean", "rms", "correlation")}
            audits = []
            for u, week in enumerate(weeks):
                for r, view in enumerate(views):
                    base = dict(atlas=atlas, cohort=cohort, week=week, view=view)
                    ids = np.flatnonzero(selected[:, u, r])
                    result = magnitude(values[:, u, r], selected[:, u, r])
                    residual[:, u, r] = result["residual_group_loo"]
                    variance[u, r], group_mean[u, r] = result["edge_variance"], result["group_mean"]
                    for key in distances:
                        distances[key][u, r] = result["distances"][key]
                    for k, i in enumerate(ids):
                        for j in ids[k+1:]:
                            tables["distances"].append(dict(**base, mouse_i=store.mice[i], mouse_j=store.mice[j],
                                n_mice=len(ids), euclidean=result["distances"]["euclidean"][i, j],
                                rms=result["distances"]["rms"][i, j],
                                correlation_distance=result["distances"]["correlation"][i, j],
                                correlation_status="defined" if np.isfinite(result["distances"]["correlation"][i, j])
                                else "constant_edge_vector"))
                    status = "complete" if len(ids) >= 2 else "fewer_than_two_mice"
                    audits.append(dict(**base, mouse_ids=[store.mice[i] for i in ids], status=status))
                    tables["dispersion"].append(dict(**base, n_mice=len(ids), n_unique_pairs=len(ids)*(len(ids)-1)//2,
                        dispersion=result["dispersion"], mean_edge_variance=float(np.mean(variance[u, r])),
                        status=status))
                    for row in regional_summaries(variance[u, r], groups, "edge_variance"):
                        tables["regional"].append(dict(**base, mouse="", **row))
                    scores = np.sqrt(np.mean(residual[ids, u, r] ** 2, axis=1))
                    for k, i in enumerate(ids):
                        vector = residual[i, u, r]
                        refs = [store.mice[j] for j in ids if j != i]
                        has_reference = ri is not None and all_valid[i, ri, r]
                        if has_reference:
                            own_change[i, u, r] = values[i, u, r] - all_values[i, ri, r]
                        rank = 1 + int(np.sum(scores > scores[k])) if np.isfinite(scores[k]) else np.nan
                        tables["subjects"].append(dict(**base, mouse=store.mice[i], n_reference_mice=len(refs),
                            reference_mice="|".join(refs), rms_from_group_loo=float(scores[k]),
                            euclidean_from_group_loo=float(scores[k]*np.sqrt(e)),
                            mean_signed_group_deviation=float(vector.mean()),
                            positive_fraction=float(np.mean(vector > 0)) if refs else np.nan,
                            negative_fraction=float(np.mean(vector < 0)) if refs else np.nan,
                            deviation_rank_descending=rank, group_status=status,
                            reference_week=reference, reference_status="available" if has_reference else "missing",
                            rms_from_reference=float(np.sqrt(np.mean(own_change[i, u, r] ** 2))),
                            mean_signed_reference_change=float(own_change[i, u, r].mean())))
                        for row in regional_summaries(vector, groups, "residual_group_loo"):
                            tables["regional"].append(dict(**base, mouse=store.mice[i], **row))
                    if len(ids) >= 2:
                        top = ids[scores == scores.max()]
                        tables["examples"].append(dict(**base, mouse=store.mice[top[0]],
                            tied_mouse_ids="|".join(store.mice[i] for i in top),
                            rms_from_group_loo=float(scores.max()), n_mice=len(ids),
                            selection="largest_LOO_RMS_config_order_tiebreak", interpretation="descriptive_example_not_outlier_test"))
            metadata = dict(schema=SCHEMA, atlas=atlas, cohort=cohort,
                mouse_ids=store.mice, week_ids=weeks, view_ids=views, roi_labels=store.roi_labels,
                run_views=definitions, scale="FC_edge_Fisher_z", edge_order="numpy.triu_indices(k=1)",
                reference_week=reference, reference_meaning="configured_own_reference_not_assumed_healthy",
                group_reference="same_cohort_week_view_leave_one_mouse_out", edge_variance_ddof=1,
                dispersion_denominator="N_mice", **mapping, edge_store=inputs[-1], code=code,
                axes=dict(distances="week,view,mouse_i,mouse_j", residual="mouse,week,view,edge",
                          edge_variance="week,view,edge", valid_record="mouse,week,view"))
            path = folder / atlas / f"{cohort}.npz"
            save_arrays(path, metadata, residual_group_loo=residual, delta_from_reference=own_change,
                edge_variance=variance, group_mean=group_mean, valid_record=selected,
                n_valid_mice=selected.sum(axis=0), edge_index=store.edge_index,
                **{f"distance_{key}": array for key, array in distances.items()})
            products.append(dict(atlas=atlas, cohort=cohort, path=str(path), sha256=sha256(path)))
            selection.extend(audits)
    for name, rows in tables.items():
        empty_columns = dict(
            distances=["atlas", "cohort", "week", "view", "mouse_i", "mouse_j", "n_mice",
                       "euclidean", "rms", "correlation_distance", "correlation_status"],
            subjects=["atlas", "cohort", "week", "view", "mouse", "rms_from_group_loo"],
            regional=["atlas", "cohort", "week", "view", "mouse", "quantity", "level", "region",
                      "mean_value", "rms", "positive_rms", "negative_rms"],
            examples=["atlas", "cohort", "week", "view", "mouse"],
            reliability=["atlas", "week", "mouse", "run_1", "run_2", "pattern_r", "pattern_status",
                         "rms_difference", "mean_signed_difference"],
            dispersion=["atlas", "cohort", "week", "view", "dispersion", "n_mice"])
        frame = pd.DataFrame(rows) if rows else pd.DataFrame(columns=empty_columns[name])
        save_table(folder / f"{name}.csv", frame, index=False)
    save_json(folder / "selection.json", dict(selections=selection,
        fixed_policy="all_selected_weeks_valid_separately_per_view", imputation="none"))
    save_json(folder / "methods.json", dict(schema=SCHEMA, scale="FC_edge_Fisher_z",
        formulas=dict(euclidean="sqrt(sum_e((f_i-f_j)^2))", rms="euclidean/sqrt(E)",
            correlation_distance="1-Pearson_over_edges(f_i,f_j)",
            residual_group_loo="f_i-mean(f_j,j!=i)",
            dispersion="mean_i(mean_e((f_i-group_mean)^2))", edge_variance="var_i(f,ddof=1)"),
        direction="positive=stronger_FC_than_same_week_LOO; negative=weaker; not_excitation_or_inhibition",
        roi_contribution="half_edge_squared_residual_or_variance_to_each_endpoint; normalize_by_total",
        network_contribution="unique_edges_per_unordered_network_pair; also_report_mean_per_edge",
        undefined="NaN; constant_pattern_correlation_undefined; N<2_variance_and_LOO_undefined",
        uncertainty="descriptive_only_no_confidence_intervals_or_outlier_hypothesis_test",
        icc="not_fitted", variance_model="not_fitted", inputs=inputs, code=code))
    files = sorted(p for p in folder.glob("*.csv")) + [folder / "selection.json", folder / "methods.json"]
    save_json(folder / "index.json", dict(schema=SCHEMA, products=products,
        files=[dict(path=str(p), sha256=sha256(p)) for p in files]))
    ctx.logger.info("Descriptive v2: %s atlas/cohort products", len(products))
