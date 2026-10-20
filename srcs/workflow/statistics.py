import numpy as np
import os
import pandas as pd
import pickle
import pingouin as pg

from scipy.stats import ttest_rel
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multitest import multipletests

from srcs.utils.config import load_context
from srcs.comp.statistics import paired_effect_summary
from srcs.utils.results import save_table, save_json
from srcs.utils.dataset import statistics_metadata
from srcs.io.writers import export_individual_identifiability
from srcs.vis.PlotSeries import PlotSeries


def run(ctx):
    logger = ctx.logger
    render_options = ctx.render_options()
    mean_error_cfg = ctx.viz_style("stat_mean_error")

    stat_df_path = ctx.statistics_path
    with open(stat_df_path, "rb") as f:
        stat = pickle.load(f)
    analysis_meta = statistics_metadata(ctx, stat)
    condition_axis = "week" if analysis_meta["split_axis"] == "run" else "run"
    within = ["condition", condition_axis] if len(stat.columns) > 1 else ["condition"]
    save_json(os.path.join(ctx.output_dir("Metadata"), "statistics_analysis.json"), analysis_meta)
    export_individual_identifiability(ctx, stat, analysis_meta)

    # Plot Group-Indiv-Diff curve
    for atlas in ctx.atlases:
        data, title = [], []
        ind = "_".join(atlas)
        logger.info(f"\n>>> Current class: {ind}")

        group = [info[0] for info in stat.loc[ind, :]]
        indiv = [info[1] for info in stat.loc[ind, :]]

        data.append((group, indiv))
        title.append(ind)

        plot_cfg = dict(mean_error_cfg)
        plot_cfg["title"] = title[0]
        plot_cfg["xlabel"] = ("Cumulative " if analysis_meta["selection_mode"] == "cumulative" else "") + condition_axis.capitalize()
        plot_cfg["xticklabels"] = list(stat.columns)
        PlotSeries(
            figsize=tuple(plot_cfg.pop("figsize", (8, 6))),
            dpi=plot_cfg.pop("dpi", ctx.figure_dpi()),
            save_path=ctx.output_dir("Statistics"),
        ).mean_error(
            data=(group, indiv),
            name=f"Line_Group-Indiv_{ind}_gradient",
            **plot_cfg,
        ).render(**render_options)

        rows = []
        for trial in stat.columns:

            cell = stat.loc[ind, trial]

            group = np.array(cell[0])
            indiv = np.array(cell[1])

            n_mouse = len(indiv) if len(indiv) == len(group) else min(len(indiv), len(group))

            for mouse in range(n_mouse):
                mouse_id = analysis_meta["mice"][mouse]

                rows.append({
                    "mouse": mouse_id,
                    "condition": "Indiv",
                    condition_axis: str(trial),
                    "similarity": indiv[mouse],
                })

                rows.append({
                    "mouse": mouse_id,
                    "condition": "Group",
                    condition_axis: str(trial),
                    "similarity": group[mouse]
                })

        df_long = pd.DataFrame(rows)
        save_table(os.path.join(ctx.output_dir("Statistics", "Tables"), f"observations_{ind}.csv"), df_long, index=False)

        anova = AnovaRM(
            data=df_long,
            depvar="similarity",
            subject="mouse",
            within=within
        )

        result = anova.fit()
        save_table(os.path.join(ctx.output_dir("Statistics", "Tables"), f"anova_{ind}.csv"), result.anova_table)

        logger.info("==================== ANOVA RESULT ====================")
        logger.info(result)

        aov = pg.rm_anova(
            data=df_long,
            dv="similarity",
            within=within,
            subject="mouse",
            detailed=True
        )

        save_table(os.path.join(ctx.output_dir("Statistics", "Tables"), f"rm_anova_{ind}.csv"), aov, index=False)
        logger.info("==================== RM_ANOVA ====================")
        logger.info(aov)

    # Post-hoc paired t-tests with FDR correction
    for trial in stat.columns:

        logger.info(f"\n>>> Current Trial: {trial}")

        exclude = set(ctx.section("analysis").get("statistics_exclude", []))
        systems = [s for s in stat.index.tolist() if s not in exclude]

        data_dict = {}
        for sys in systems:

            cell = stat.loc[sys, trial]

            data_dict[sys] = {
                "group": np.array(cell[0]),
                "indiv": np.array(cell[1]),
                "rel": np.array(cell[2]),
            }

        results, pvals = [], []
        metrics = ["group", "indiv", "rel"]

        for metric in metrics:
            for i in range(len(systems)):
                for j in range(i + 1, len(systems)):

                    sys1 = systems[i]
                    sys2 = systems[j]

                    a = data_dict[sys1][metric]
                    b = data_dict[sys2][metric]

                    t, p = ttest_rel(a, b)
                    effect = paired_effect_summary(a, b)

                    results.append({
                        "metric": metric,
                        "comparison": f"{sys1} vs {sys2}",
                        "t": t,
                        "p_uncorrected": p,
                        **effect,
                    })

                    pvals.append(p)

        reject, p_corr = (multipletests(pvals, method="fdr_bh")[0:2] if pvals else ([], []))

        for i in range(len(results)):
            results[i]["p_FDR"] = p_corr[i]
            results[i]["significant"] = reject[i]

        df_results = pd.DataFrame(results)
        save_table(os.path.join(ctx.output_dir("Statistics", "Tables"), f"posthoc_{trial}.csv"), df_results, index=False)

        logger.info("==================== Statistical Test Results ====================")
        logger.info(df_results)

        logger.info("==================== PAPER STYLE ====================")

        for _, row in df_results.iterrows():

            metric = row["metric"]
            comp = row["comparison"]
            t = row["t"]
            p = row["p_FDR"]
            d = row["cohen_d"]

            dfree = len(data_dict[systems[0]][metric]) - 1

            if p < 0.001:
                p_text = "P < 0.001"
            else:
                p_text = f"P = {p:.3f}"

            logger.info(
                f"{metric}: {comp}: "
                f"t({dfree}) = {t:.2f}, "
                f"{p_text}, "
                f"paired Cohen's dz = {d:.2f}"
            )
