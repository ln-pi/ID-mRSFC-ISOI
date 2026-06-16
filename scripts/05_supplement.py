import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from statsmodels.stats.anova import AnovaRM
from statsmodels.multivariate.manova import MANOVA

from srcs.vis.PlotSeries import PlotSeries


class StructureAnalyzer:

    def __init__(self, DD, ND, group_names=("A", "B", "C")):

        self.groups = group_names

        self.DD = {
            group_names[0]: np.asarray(DD[0]),
            group_names[1]: np.asarray(DD[1]),
            group_names[2]: np.asarray(DD[2])
        }

        self.ND = {
            group_names[0]: np.asarray(ND[0]),
            group_names[1]: np.asarray(ND[1]),
            group_names[2]: np.asarray(ND[2])
        }

        self.n_trials = len(DD[0])

        self.df_rm = self._build_rm_dataframe()
        self.df_manova = self._build_manova_dataframe()

    # =====================================================
    # dataframe
    # =====================================================

    def _build_rm_dataframe(self):

        rows = []

        for trial in range(self.n_trials):

            for group in self.groups:

                rows.append({
                    "trial": trial,
                    "group": group,
                    "DD": self.DD[group][trial],
                    "ND": self.ND[group][trial]
                })

        return pd.DataFrame(rows)

    def _build_manova_dataframe(self):

        rows = []

        for group in self.groups:

            for trial in range(self.n_trials):

                rows.append({
                    "group": group,
                    "DD": self.DD[group][trial],
                    "ND": self.ND[group][trial]
                })

        return pd.DataFrame(rows)

    # =====================================================
    # RM ANOVA
    # =====================================================

    def rm_anova_dd(self):

        aov = AnovaRM(
            data=self.df_rm,
            depvar="DD",
            subject="trial",
            within=["group"]
        )

        return aov.fit()

    def rm_anova_nd(self):

        aov = AnovaRM(
            data=self.df_rm,
            depvar="ND",
            subject="trial",
            within=["group"]
        )

        return aov.fit()

    # =====================================================
    # MANOVA
    # =====================================================

    def manova(self):

        return MANOVA.from_formula(
            "DD + ND ~ group",
            data=self.df_manova
        ).mv_test()

    # =====================================================
    # Boxplot
    # =====================================================

    def plot_dd_box(self,save_path=None):

        PlotSeries(figsize=(5, 4), dpi=180, save_path=save_path).add_plot(
            [self.DD[g] for g in self.groups],
            plot_type="box",
            name="boxplot_DD",
            config={
                "box_width": 0.1,
                "box_alpha": 0.7,
                "box_color": None,
                "showfliers": False,
                "fontsize": 7,
                "rotation": 0,
                "xticklabel": self.groups,
                "grid": False,
                "ylabel": "DD",
                "title": "Diagonal Dominance (DD)",
                "show_legend": False,
                "show_sub_label": True,
                "hide_ns": True,
            }
        ).render(show=False, save=True, log=False)

    def plot_nd_box(self,save_path=None):

        PlotSeries(figsize=(5, 4), dpi=180, save_path=save_path).add_plot(
            [self.ND[g] for g in self.groups],
            plot_type="box",
            name="boxplot_ND",
            config={
                "box_width": 0.1,
                "box_alpha": 0.7,
                "box_color": None,
                "showfliers": False,
                "fontsize": 7,
                "rotation": 0,
                "xticklabel": self.groups,
                "grid": False,
                "ylabel": "ND",
                "title": "Normalized Distance (ND)",
                "show_legend": False,
                "show_sub_label": True,
                "hide_ns": True,
            }
        ).render(show=False, save=True, log=False)

    # =====================================================
    # DD-ND Structure Space
    # =====================================================

    def plot_structure_space(self, save_path=None):

        plt.figure(figsize=(6, 5))

        for g in self.groups:

            plt.scatter(
                self.DD[g],
                self.ND[g],
                label=g,
                s=60
            )

        plt.xlabel("DD")
        plt.ylabel("ND")

        plt.title("Structure Space")

        plt.legend()

        plt.tight_layout()

        if save_path:
            plt.savefig(f"{save_path}/structure_space.png", dpi=180)
        else:
            plt.savefig("structure_space.png", dpi=180)

    # =====================================================
    # summary
    # =====================================================

    def summary(self):

        print("\n========== RM-ANOVA DD ==========")
        print(self.rm_anova_dd())

        print("\n========== RM-ANOVA ND ==========")
        print(self.rm_anova_nd())

        print("\n========== MANOVA ==========")
        print(self.manova())


if __name__ == "__main__":
    import pickle

    output_path = r"F:\Hou.T.N\ID-mRSFC-ISOI\data\output\WT\Stat"
    stat_df_path = r"F:\Hou.T.N\ID-mRSFC-ISOI\data\output\WT\statistics.pkl"
    with open(stat_df_path, "rb") as f:
        stat = pickle.load(f)
    
    Atlas = [["SM", "VIS", "AUD"], ["PF", "LAT", "MED"], ["None"]]

    DD = [[], [], []]
    ND = [[], [], []]
    for trial in stat.columns:

        exclude = "AUDv"
        systems = [s for s in stat.index.tolist() if s not in exclude]

        data_dict = {}
        for sys in systems:

            cell = stat.loc[sys, trial]

            data_dict[sys] = {
                "group": np.array(cell[0]),
                "indiv": np.array(cell[1]),
                "rel": np.array(cell[2]),
                "DD": cell[3]["diagonal_dominance"],
                "ND": cell[3]["normalized_distance_to_diagonal"],
            }

        for sys in systems:
            if sys in "_".join(Atlas[0]):
                idx = 0
            elif sys in "_".join(Atlas[1]):
                idx = 1
            else:
                idx = 2

            DD[idx].append(data_dict[sys]["DD"])
            ND[idx].append(data_dict[sys]["ND"])
    
    analyzer = StructureAnalyzer(DD, ND, group_names=["Association", "Sensory", "Whole Brain"])
    analyzer.summary()
    analyzer.plot_dd_box(save_path=output_path)
    analyzer.plot_nd_box(save_path=output_path)
    analyzer.plot_structure_space(save_path=output_path)
    
