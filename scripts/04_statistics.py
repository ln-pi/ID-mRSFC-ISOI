import numpy as np
import os
import pandas as pd
import pickle
import pingouin as pg

from scipy.stats import ttest_rel
from statsmodels.stats.anova import AnovaRM
from statsmodels.stats.multitest import multipletests

from srcs.io.readers import load_config
from srcs.utils.utils import parse_args
from srcs.utils.logger import setup_logger
from srcs.vis.PlotSeries import PlotSeries


args = parse_args()
config = load_config(args.config)
if config is None:
    raise ValueError("Config is empty (None). Exiting.")
logger = setup_logger(config["paths"]["log"])

processed_root = config["paths"]["processed"]
output_root = config["paths"]["output"]

dataset = config["variables"]["dataset"]
Session_ids = config["variables"]["Session_ids"]
Mouse_ids = config["variables"]["Mouse_ids"]
Trial_ids = config["variables"]["Trial_ids"]
Atlas = config["variables"]["Atlas"]

processed_path = os.path.join(processed_root, dataset)
output_path = os.path.join(output_root, dataset)
os.makedirs(output_path, exist_ok=True)

stat_df_path = os.path.join(output_path, "statistics.pkl")
with open(stat_df_path, "rb") as f:
    stat = pickle.load(f)

# Plot Group-Indiv-Diff curve
for atlas in Atlas:
    data, title = [], []
    ind = "_".join(atlas)
    logger.info(f"\n>>> Current class: {ind}")

    group = [info[0] for info in stat.loc[ind, :]]
    indiv = [info[1] for info in stat.loc[ind, :]]

    data.append((group, indiv))
    title.append(ind)

    PlotSeries(figsize=(8, 6), dpi=180, save_path=os.path.join(output_path, "Stat")).add_plot(
        data=(group, indiv),
        plot_type="mean_error",
        name=f"Line_Group-Indiv_{ind}_gradient",
        config={
            "error_type": "sem",
            "use_fill": False,
            "labels": ("Group", "Individual"),
            "show_diff": True,
            "diff_mode": "gradient",
            "diff_ylabel": "Gradient",
            "xlabel": "Trials",
            "ylabel": "Value",
            "title": title[0]
        }
    ).render(show=False, save=True, log=False)

    rows = []
    for trial in stat.columns:

        cell = stat.loc[ind, trial]

        group = np.array(cell[0])
        indiv = np.array(cell[1])

        n_mouse = len(indiv) if len(indiv) == len(group) else min(len(indiv), len(group))

        for mouse in range(n_mouse):

            rows.append({
                "mouse": mouse,
                "condition": "Indiv",
                "trial": str(trial),
                "similarity": indiv[mouse],
            })

            rows.append({
                "mouse": mouse,
                "condition": "Group",
                "trial": str(trial),
                "similarity": group[mouse]
            })
    
    df_long = pd.DataFrame(rows)

    anova = AnovaRM(
        data=df_long,
        depvar="similarity",
        subject="mouse",
        within=["condition", "trial"]
    )

    result = anova.fit()

    logger.info("==================== ANOVA RESULT ====================")
    logger.info(result)

    aov = pg.rm_anova(
        data=df_long,
        dv="similarity",
        within=["condition", "trial"],
        subject="mouse",
        detailed=True
    )

    logger.info("==================== RM_ANOVA ====================")
    logger.info(aov)

# Post-hoc paired t-tests with FDR correction
for trial in stat.columns:

    logger.info(f"\n>>> Current Trial: {trial}")

    exclude = "AUDv"
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
                diff = np.array(a) - np.array(b)
                d = np.mean(diff) / np.std(diff)

                results.append({
                    "metric": metric,
                    "comparison": f"{sys1} vs {sys2}",
                    "t": t,
                    "p_uncorrected": p,
                    "cohen_d": d,
                })

                pvals.append(p)
    
    reject, p_corr, _, _ = multipletests(pvals, method="fdr_bh")

    for i in range(len(results)):
        results[i]["p_FDR"] = p_corr[i]
        results[i]["significant"] = reject[i]
    
    df_results = pd.DataFrame(results)

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
            f"Cohen's d = {d:.2f}"
        )
