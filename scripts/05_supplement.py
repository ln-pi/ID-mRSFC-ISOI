import os
import pickle

import numpy as np
import pandas as pd

from statsmodels.multivariate.manova import MANOVA
from statsmodels.stats.anova import AnovaRM

from srcs.utils.config import load_context
from srcs.utils.results import save_pickle, save_table
from srcs.vis.PlotSeries import PlotSeries


class StructureAnalyzer:
    def __init__(self, DD, ND, group_names=("A", "B", "C"), logger=None):
        self.groups = list(group_names)
        self.logger = logger

        self.DD = {
            group: np.asarray(values)
            for group, values in zip(self.groups, DD)
        }
        self.ND = {
            group: np.asarray(values)
            for group, values in zip(self.groups, ND)
        }

        self.n_trials = len(DD[0])
        self.df_rm = self._build_rm_dataframe()
        self.df_manova = self._build_manova_dataframe()

    def _build_rm_dataframe(self):
        rows = []
        for trial in range(self.n_trials):
            for group in self.groups:
                rows.append({
                    "trial": trial,
                    "group": group,
                    "DD": self.DD[group][trial],
                    "ND": self.ND[group][trial],
                })
        return pd.DataFrame(rows)

    def _build_manova_dataframe(self):
        rows = []
        for group in self.groups:
            for trial in range(self.n_trials):
                rows.append({
                    "group": group,
                    "DD": self.DD[group][trial],
                    "ND": self.ND[group][trial],
                })
        return pd.DataFrame(rows)

    def rm_anova_dd(self):
        return AnovaRM(
            data=self.df_rm,
            depvar="DD",
            subject="trial",
            within=["group"],
        ).fit()

    def rm_anova_nd(self):
        return AnovaRM(
            data=self.df_rm,
            depvar="ND",
            subject="trial",
            within=["group"],
        ).fit()

    def manova(self):
        return MANOVA.from_formula("DD + ND ~ group", data=self.df_manova).mv_test()

    def plot_box(self, metric, save_path=None, figsize=(5, 4), dpi=180, show=False, save=True, log=False, config=None):
        metric = metric.upper()
        values = self.DD if metric == "DD" else self.ND
        cfg = dict(config or {})
        cfg.setdefault("xticklabel", self.groups)
        cfg.setdefault("ylabel", metric)

        return (
            PlotSeries(figsize=figsize, dpi=dpi, save_path=save_path)
            .box(
                [values[g] for g in self.groups],
                name=f"boxplot_{metric}",
                **cfg,
            )
            .render(show=show, save=save, log=log)
        )

    def plot_structure_space(self, save_path=None, figsize=(6, 5), dpi=180, show=False, save=True, log=False, config=None):
        data = [
            {"x": self.DD[group], "y": self.ND[group], "label": group}
            for group in self.groups
        ]
        return (
            PlotSeries(figsize=figsize, dpi=dpi, save_path=save_path)
            .scatter(data, name="structure_space", **dict(config or {}))
            .render(show=show, save=save, log=log)
        )

    def summary(self):
        outputs = {
            "RM-ANOVA DD": self.rm_anova_dd(),
            "RM-ANOVA ND": self.rm_anova_nd(),
            "MANOVA": self.manova(),
        }
        for title, result in outputs.items():
            message = f"\n========== {title} ==========\n{result}"
            if self.logger is None:
                print(message)
            else:
                self.logger.info(message)
        return outputs


def _atlas_tag(atlas):
    return "_".join(atlas)


def build_structure_inputs(stat, group_atlases, exclude=None):
    exclude = set(exclude or [])
    group_tags = [_atlas_tag(atlas) for atlas in group_atlases]
    DD = [[] for _ in group_tags]
    ND = [[] for _ in group_tags]

    for trial in stat.columns:
        systems = [system for system in stat.index.tolist() if system not in exclude]
        data = {
            system: {
                "DD": stat.loc[system, trial][3]["diagonal_dominance"],
                "ND": stat.loc[system, trial][3]["normalized_distance_to_diagonal"],
            }
            for system in systems
        }

        for idx, tag in enumerate(group_tags):
            if tag not in data:
                continue
            DD[idx].append(data[tag]["DD"])
            ND[idx].append(data[tag]["ND"])

    return DD, ND


def main():
    ctx = load_context()
    logger = ctx.logger
    supplement_cfg = ctx.section("supplement")
    render_options = ctx.render_options()

    stat_df_path = ctx.statistics_path
    with open(stat_df_path, "rb") as f:
        stat = pickle.load(f)

    group_atlases = supplement_cfg.get(
        "structure_atlases",
        [atlas for atlas in ctx.atlases if _atlas_tag(atlas) != "AUDv"],
    )
    group_names = supplement_cfg.get(
        "structure_group_names",
        [_atlas_tag(atlas) for atlas in group_atlases],
    )
    exclude = supplement_cfg.get("exclude", ["AUDv"])

    DD, ND = build_structure_inputs(stat, group_atlases, exclude=exclude)
    analyzer = StructureAnalyzer(DD, ND, group_names=group_names, logger=logger)
    summaries = analyzer.summary()
    table_path = ctx.output_dir("Stat", "Tables", "Supplement")
    save_table(os.path.join(table_path, "observations.csv"), analyzer.df_rm, index=False)
    save_table(os.path.join(table_path, "rm_anova_dd.csv"), summaries["RM-ANOVA DD"].anova_table)
    save_table(os.path.join(table_path, "rm_anova_nd.csv"), summaries["RM-ANOVA ND"].anova_table)
    save_pickle(os.path.join(table_path, "manova.pkl"), summaries["MANOVA"].results)
    for term, result in summaries["MANOVA"].results.items():
        save_table(os.path.join(table_path, f"manova_{term}.csv"), result["stat"])

    save_path = ctx.output_dir("Stat")

    dd_box_cfg = ctx.viz_style("stat_box", "supplement_dd_box")
    nd_box_cfg = ctx.viz_style("stat_box", "supplement_nd_box")
    space_cfg = ctx.viz_style("supplement_structure_space")

    analyzer.plot_box(
        "DD",
        save_path=save_path,
        figsize=tuple(dd_box_cfg.pop("figsize", (5, 4))),
        dpi=dd_box_cfg.pop("dpi", ctx.figure_dpi()),
        config=dd_box_cfg,
        **render_options,
    )
    analyzer.plot_box(
        "ND",
        save_path=save_path,
        figsize=tuple(nd_box_cfg.pop("figsize", (5, 4))),
        dpi=nd_box_cfg.pop("dpi", ctx.figure_dpi()),
        config=nd_box_cfg,
        **render_options,
    )
    analyzer.plot_structure_space(
        save_path=save_path,
        figsize=tuple(space_cfg.pop("figsize", (6, 5))),
        dpi=space_cfg.pop("dpi", ctx.figure_dpi()),
        config=space_cfg,
        **render_options,
    )


if __name__ == "__main__":
    main()
