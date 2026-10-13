import numpy as np
import os
import pandas as pd

from statsmodels.stats.diagnostic import lilliefors
from scipy.stats import ttest_rel, wilcoxon, pearsonr
from tqdm import tqdm

from srcs.utils.utils import df_to_array
from srcs.comp.statistics import paired_effect_summary, monte_carlo_pvalue_greater
from srcs.utils.logger import setup_logger
from srcs.vis.PlotSeries import PlotSeries


logger = setup_logger()


class Identify(object):

    name = "Identify correlation matrices"

    def __init__(self, df, atlas_exclude=None):
        self.df = df
        self.ticklabels = df.columns.tolist()
        self.labels = df.index.tolist()

        self.atlas = atlas_exclude

        self.sim_mat:np.ndarray
        self.Indiv:np.ndarray
        self.Group:np.ndarray
        self.relEffect:np.ndarray
        self.cohen_d:float
    
    @staticmethod
    def remove_nan(matrix):

        if matrix.ndim != 4:
            logger.error("Input should be (Run, Sub, N, N)")
            raise ValueError("Input should be (Run, Sub, N, N)")
        
        nan_row = np.all(np.isnan(matrix), axis=3) # (R, S, N)
        bad_nodes = np.any(nan_row, axis=(0, 1)) # (N,)

        good_nodes = ~bad_nodes
        kept_indices = np.where(good_nodes)[0]

        clean_matrix = matrix[:, :, good_nodes, :][:, :, :, good_nodes]

        return clean_matrix, kept_indices

    @staticmethod
    def Calc_NetworkSimilarity(series1, series2):
        """Compare aligned, symmetric half-average matrices (subject x ROI x ROI).

        The caller first Fisher-z transforms each FC edge and averages within
        each half (build_half_dataframe). Do not transform those edges again.
        Here each undirected off-diagonal edge is counted once, Pearson r is
        computed between edge vectors, then that *similarity* is Fisher-z
        transformed and capped at +/-4 for compatibility with existing outputs.
        ROI order must already match in both inputs. Diagonals are ignored.
        """
        series1 = np.asarray(series1, dtype=float)
        series2 = np.asarray(series2, dtype=float)
        if (series1.shape != series2.shape or series1.ndim != 3
                or series1.shape[1] != series1.shape[2]):
            raise ValueError("Inputs must have matching (subject, ROI, ROI) shapes.")
        sub_num, node_num, _ = series1.shape
        if sub_num < 2 or node_num < 3:
            raise ValueError("Similarity requires at least two subjects and three ROIs.")
        iu = np.triu_indices(node_num, k=1)
        vectors = []
        for series in (series1, series2):
            edges = series[:, iu[0], iu[1]]
            reverse = series[:, iu[1], iu[0]]
            if not np.isfinite(edges).all() or not np.isfinite(reverse).all():
                raise ValueError("Off-diagonal edges must be finite.")
            if not np.allclose(edges, reverse, rtol=1e-7, atol=1e-10):
                raise ValueError("Connectivity matrices must be symmetric.")
            if np.any(np.ptp(edges, axis=1) == 0):
                raise ValueError("Constant edge vectors have undefined Pearson similarity.")
            vectors.append(edges)
        sim_mat = np.zeros((sub_num, sub_num))
        r_limit = np.tanh(4.0)
        for i in range(sub_num):
            for j in range(sub_num):
                r = np.corrcoef(vectors[0][i], vectors[1][j])[0, 1]
                if not np.isfinite(r):
                    raise ValueError("Non-finite Pearson similarity; check edge variation.")
                # Equivalent cap, applied before arctanh to avoid infinities at +/-1.
                sim_mat[i, j] = np.clip(np.arctanh(np.clip(r, -r_limit, r_limit)), -4, 4)
        
        Indiv = np.diag(sim_mat)
        Group = np.zeros(sub_num)

        for i in range(sub_num):
            target = np.delete(np.arange(sub_num), i)
            Group[i] = np.nanmean(np.concatenate([sim_mat[i, target], sim_mat[target, i]]))
        
        relEffect = 1 - np.divide(Group, Indiv, out=np.full_like(Indiv, np.nan), where=Indiv != 0)
        cohen_d = paired_effect_summary(Indiv, Group)["cohen_dz"]

        return sim_mat, Indiv, Group, relEffect, cohen_d
    
    @staticmethod
    def Calc_IdentificationAccuracy(sim_mat):
        sub_num = sim_mat.shape[0]
        correct, row_correct, col_correct = 0, 0, 0

        for i in range(sub_num):
            group = np.delete(np.arange(sub_num), i)

            if np.all(sim_mat[i, i] > sim_mat[i, group]) and np.all(sim_mat[i, i] > sim_mat[group, i]):
                correct += 1
            
            if np.all(sim_mat[i, i] > sim_mat[i, group]):
                row_correct += 1
            
            if np.all(sim_mat[i, i] > sim_mat[group, i]):
                col_correct += 1
        
        accuracy = correct / sub_num
        row_accuracy = row_correct / sub_num
        col_accuracy = col_correct / sub_num

        return accuracy, row_accuracy, col_accuracy
    
    @staticmethod
    def Calc_Entropy(sim_mat, baseline=None):
        M = np.array(sim_mat, dtype=float)
        M -= baseline if baseline is not None else M.min()
        M[M <0] = 0

        total = M.sum()

        # Diagonal dominance
        diagonal_dominance = np.trace(M) / total

        # Distance to diagonal
        I, J = np.indices(M.shape)
        dist = np.abs(I - J)
        weighted_distance = np.sum(M * dist) / total
        max_dist = M.shape[0] - 1
        normalized_distance = weighted_distance / max_dist

        return diagonal_dominance, normalized_distance

    def ID_Test(self):
        """
        Input should be pandas.DataFrame (index = Runs, columns = Subjects)
        Filled in each cell should be corr_matrices sharing same shape (N x N)
        """

        corrMat = self.df
        matrices, ind, cols, labels = df_to_array(corrMat)

        matrices =matrices.squeeze()
        matrices, kept = self.remove_nan(matrices)
        labels = [labels[i] for i in kept]

        sim_mat, Indiv, Group, relEffect, cohen_d = self.Calc_NetworkSimilarity(matrices[0], matrices[1])

        self.sim_mat = sim_mat
        self.Indiv = Indiv
        self.Group = Group
        self.relEffect = relEffect
        self.cohen_d = cohen_d

        return [ind, cols, labels], sim_mat, Indiv, Group, relEffect, cohen_d
    
    def remove_bad_subject(self, obj:int):
        self.sim_mat = np.delete(self.sim_mat, obj, axis=0)
        self.sim_mat = np.delete(self.sim_mat, obj, axis=1)
        self.Indiv = np.delete(self.Indiv, obj)
        self.Group = np.delete(self.Group, obj)
        self.relEffect = np.delete(self.relEffect, obj)
    
    def statistics(self):
        
        accuracy, row_accuracy, col_accuracy = self.Calc_IdentificationAccuracy(self.sim_mat)
        logger.info(f"Identification Accuracy: {accuracy:.2%}")
        logger.info(f"Identification Accuracy in row: {row_accuracy:.2%}")
        logger.info(f"Identification Accuracy in column: {col_accuracy:.2%}")

        diag_dom, norm_dist = self.Calc_Entropy(self.sim_mat)
        logger.info(f"Diagonal Dominance: {diag_dom:.4f}")
        logger.info(f"Normalized Distance to Diagonal: {norm_dist:.4f}")

        stat_l_i, p_l_i = lilliefors(self.Indiv)
        stat_l_g, p_l_g = lilliefors(self.Group)
        logger.info(f"Lilliefors test for Indiv: stat = {stat_l_i}, p = {p_l_i}")
        logger.info(f"Lilliefors test for Group: stat = {stat_l_g}, p = {p_l_g}")

        stat_t, p_t = ttest_rel(self.Indiv, self.Group)
        logger.info(f"Paired t-test: t = {stat_t:.4f}, p = {p_t:.4e}")

        stat_w, p_w = wilcoxon(self.Indiv, self.Group)
        logger.info(f"Wilcoxon signed-rank test: stat = {stat_w:.4f}, p = {p_w:.4e}")

        r, p_r = pearsonr(self.Indiv, self.Group)
        logger.info(f"Pearson correlation: r = {r:.4f}, p = {p_r:.4e}")

        return {
            "identification_accuracy": accuracy,
            "row_identification_accuracy": row_accuracy,
            "col_identification_accuracy": col_accuracy,
            "diagonal_dominance": diag_dom,
            "normalized_distance_to_diagonal": norm_dist,
            "Indiv_lilliefors_statistic": stat_l_i,
            "Indiv_lilliefors_p_value": p_l_i,
            "Group_lilliefors_statistic": stat_l_g,
            "Group_lilliefors_p_value": p_l_g,
            "t_test_statistic": stat_t,
            "t_test_p_value": p_t,
            "wilcoxon_statistic": stat_w,
            "wilcoxon_p_value": p_w,
            "pearson_correlation": r,
            "pearson_p_value": p_r,
            **paired_effect_summary(self.Indiv, self.Group),
            "effect_size_direction": "individual_minus_group",
        }
    
    def shuffle(self, iter=1000, figsize=(5, 4), dpi=180, save_path=None, show=True, log=False, tag="shuffle", plot_config=None, seed=None):
        if isinstance(iter, (bool, np.bool_)) or not isinstance(iter, (int, np.integer)) or iter < 1:
            raise ValueError("Shuffle iter must be a positive integer.")
        plot_config = dict(plot_config or {})
        figsize = tuple(plot_config.pop("figsize", figsize))
        dpi = plot_config.pop("dpi", dpi)
        atlas_tag = "_".join(self.atlas)
        acc, r_acc, c_acc = [], [], []
        rng = np.random.default_rng(seed)
        real_acc, real_r_acc, real_c_acc = self.Calc_IdentificationAccuracy(self.sim_mat)

        for _ in tqdm(range(iter), total=iter, desc="Shuffle", leave=False):
            shuffled = self.df.copy()
            shuffled_mouse = rng.permutation(self.df.columns)
            shuffled.loc["2nd"] = self.df.loc["2nd", shuffled_mouse].values

            tmp = Identify(shuffled, atlas_exclude=self.atlas)
            tmp.ID_Test()
            a, ra, ca = tmp.Calc_IdentificationAccuracy(tmp.sim_mat)
            acc.append(a)
            r_acc.append(ra)
            c_acc.append(ca)
        
        acc = np.array(acc)
        r_acc = np.array(r_acc)
        c_acc = np.array(c_acc)

        permutation = monte_carlo_pvalue_greater(real_acc, acc)
        row_permutation = monte_carlo_pvalue_greater(real_r_acc, r_acc)
        col_permutation = monte_carlo_pvalue_greater(real_c_acc, c_acc)
        logger.info(f"Real Accuracy = {real_acc:.4f}")
        logger.info(f"Shuffle Accuracy = {np.mean(acc):.4f} +/- {np.std(acc):.4f}")
        logger.info(f"Permutation p-value (+1, greater) = {permutation['p_perm']:.6g}; "
                    f"extreme={permutation['n_extreme']}/{permutation['n_iter']}")

        if save_path is None:
            save = False
        else:
            save = True
        
        config = {
            "bins": 20,
            "range": (0, 1),
            "facecolor": "lightgray",
            "edgecolor": "dimgray",
            "arrow_x": real_acc,
            "arrow_color": "red",
            "show_mean_arrow": True,
            "mean_arrow_color": "gray",
            "xlabel": "Identification Accuracy",
            "ylabel": "Count",
            "title": "Shuffle Distribution",
            "arrow_text_format": "real={:.4f}",
            "mean_text_format": "shuffle={:.4f}",
        }
        config.update(plot_config)
        config["arrow_x"] = real_acc

        PlotSeries(figsize=figsize, dpi=dpi, save_path=save_path).add_plot(
            data=acc,
            plot_type="hist",
            name=f"ID_shuffle_acc_{atlas_tag}_excluded_{tag}",
            config=config,
        ).render(show=show, save=save, log=log)

        return {
            "real_accuracy": (real_acc, real_r_acc, real_c_acc),
            "shuffle_acc": (acc, r_acc, c_acc),
            "shuffle_mean": np.mean(acc),
            "shuffle_std": np.std(acc),
            **permutation,
            "seed": seed,
            "row_permutation": row_permutation,
            "col_permutation": col_permutation,
        }

    def output(self, figsize=None, dpi=180, tag="default", save_path=None, show=True, log=False, paired_config=None, matrix_config=None):
        paired_config = dict(paired_config or {})
        matrix_config = dict(matrix_config or {})
        atlas_tag = "_".join(self.atlas)
        if figsize is None:
            figsize = {
                "line": paired_config.pop("figsize", (4, 6)),
                "mat": matrix_config.pop("figsize", (6, 6)),
            }
        logger.info(f"Paired Cohen's dz (ddof=1): {self.cohen_d}")
        logger.info(f"Similarity Matrix: min = {np.min(self.sim_mat)}, max = {np.max(self.sim_mat)}")

        if save_path is None:
            save = False
        else:
            save = True

        line_figsize = tuple(paired_config.pop("figsize", figsize["line"]))
        line_dpi = paired_config.pop("dpi", dpi)
        line_config = {
            "linecolor": "gray",
            "alpha": 0.5,
            "linewidth": 1.5,
            "meancolor": "red",
            "meanwidth": 3,
            "xticklabels": ["Group", "Individual"],
            "ylabel": "Network Similarity"
        }
        line_config.update(paired_config)

        PlotSeries(figsize=line_figsize, dpi=line_dpi, save_path=save_path).add_plot(
            data=[self.Group, self.Indiv],
            plot_type="paired_vertical",
            name=f"IDtest_paired_{atlas_tag}_excluded_{tag}",
            config=line_config,
        ).render(show=show, save=save, log=log)

        mat_figsize = tuple(matrix_config.pop("figsize", figsize["mat"]))
        mat_dpi = matrix_config.pop("dpi", dpi)
        tick_fontsize = matrix_config.pop("tick_fontsize", 10)
        vmin_default = matrix_config.pop("vmin_default", 1.2)
        vmax_default = matrix_config.pop("vmax_default", 3.0)
        vmin_association = matrix_config.pop("vmin_association", 2.0)
        vmax_association = matrix_config.pop("vmax_association", 3.5)
        mat_config = {
            "title": "Similarity Matrix",
            "cmap": "inferno",
            "cbar": True,
            "axis": True,
            "vmin": vmin_association if atlas_tag == "SM_VIS_AUD" else vmin_default,
            "vmax": vmax_association if atlas_tag == "SM_VIS_AUD" else vmax_default,
            "xlabel": self.labels[0],
            "ylabel": self.labels[1],
            "xticks": {
                "ticks": np.arange(len(self.ticklabels)),
                "labels": self.ticklabels,
                "rotation": 90,
                "fontsize": tick_fontsize
            },
            "yticks": {
                "ticks": np.arange(len(self.ticklabels)),
                "labels": self.ticklabels,
                "rotation": 0,
                "fontsize": tick_fontsize
            }
        }
        mat_config.update(matrix_config)

        PlotSeries(figsize=mat_figsize, dpi=mat_dpi, save_path=save_path).add_plot(
            self.sim_mat,
            plot_type="origin",
            name=f"IDtest_simMat_{atlas_tag}_excluded_{tag}",
            config=mat_config,
        ).render(show=show, save=save, log=log)



class Edgewise_Contribution(object):

    name = "Quantifying edgewise contribution to identification"

    def __init__(self, path, mouse_ids, sessions, seed=None, atlas_exclude=None,
                 trial_cnt=0, trial_ids=None, frames=None, group_indices=None):
        self.path = path
        self.mouse_ids = mouse_ids
        self.sessions = sessions
        self.rng = np.random.default_rng(seed)
        self.atlas = atlas_exclude

        if isinstance(trial_cnt, bool) or not isinstance(trial_cnt, int) or trial_cnt < 0:
            raise ValueError("trial_cnt must be a nonnegative integer.")
        self.trial_cnt = trial_cnt
        self.trial_ids = trial_ids
        self.roi_labels = None
        self.data = {}
        self.frames = frames
        self.group_indices = group_indices
        self._load_data()
    
    def _load_data(self):
        atlas_tag = "_".join(self.atlas)
        
        for session in tqdm(self.sessions, desc="Sessions", position=0, leave=False):
            
            file_path = os.path.join(self.path, f"{session}_correlation_matrices_{atlas_tag}_excluded.pkl")
            if self.frames is not None:
                df = self.frames[session]
            else:
                if not os.path.exists(file_path):
                    logger.warning(f"No file for session: {session}")
                    continue
                df = pd.read_pickle(file_path)

            trials = list(self.trial_ids) if self.trial_ids is not None else list(df.columns)
            if self.trial_cnt > len(trials):
                raise ValueError("trial_cnt exceeds available trials.")
            trials = trials[:self.trial_cnt] if self.trial_cnt else trials
            if any(trial not in df.columns for trial in trials):
                raise ValueError(f"Requested trials missing in {session}.")
            self.data[session] = {}

            for mouse in tqdm(self.mouse_ids, desc="Mice", position=1, leave=False):

                if mouse not in df.index:
                    logger.warning(f"Mouse {mouse} not found in {session}")
                    continue

                mats = []

                for trial in tqdm(trials, desc="Trials", position=2, leave=False):

                    cell = df.loc[mouse, trial]

                    if cell is None:
                        continue

                    if isinstance(cell, dict):
                        mat = cell.get(0, None)
                    elif isinstance(cell, (list, tuple)) and len(cell) > 0:
                        mat = cell[0]
                    else:
                        mat = cell
                    
                    if mat is None:
                        continue

                    labels = mat.index.tolist() if isinstance(mat, pd.DataFrame) else list(range(np.asarray(mat).shape[0]))
                    if self.roi_labels is None:
                        self.roi_labels = labels
                    elif labels != self.roi_labels:
                        raise ValueError("Edge matrices must have matching ROI labels.")
                    mats.append(np.asarray(mat))
                
                if len(mats) == 0:
                    logger.warning(f"No valid matrices for {mouse} in {session}")
                    continue

                mats = np.stack(mats, axis=0)
                self.data[session][mouse] = mats
        
    @staticmethod
    def _flattrn(mat):
        iu = np.triu_indices_from(mat, k=1)
        return mat[iu], iu
    
    def _split_trials(self, mats):

        T = mats.shape[0]
        if self.group_indices is not None:
            first, second = self.group_indices
            return np.mean(mats[first], axis=0), np.mean(mats[second], axis=0)
        idx = self.rng.permutation(T)

        half = T // 2
        idx1 = idx[:half]
        idx2 = idx[half:]

        R1 = np.nanmean(mats[idx1], axis=0)
        R2 = np.nanmean(mats[idx2], axis=0)

        return R1, R2
    
    def build_pairs(self):

        pairs = {}

        for session, mice_data in tqdm(self.data.items(), total=len(self.data), desc="Sessions", position=0, leave=False):

            pairs[session] = {}

            for mouse, mats in tqdm(mice_data.items(), total=len(mice_data), desc=f"{session}", position=1, leave=False):
                if mats.shape[0] < 2:
                    continue

                R1, R2 = self._split_trials(mats)
                pairs[session][mouse] = (R1, R2)
        
        return pairs
    
    def compute_phi_dp(self, pairs):
        Phi, DP = {}, {}

        for session, mice_data in tqdm(pairs.items(), total=len(pairs), desc="Sessions", position=0, leave=False):
            mice = list(mice_data.keys())
            N = len(mice)

            if N < 2:
                continue

            R1, R2 = [], []

            for m in mice:
                r1, r2 = mice_data[m]
                r1, self.edge_index = self._flattrn(r1)
                r2, self.edge_index = self._flattrn(r2)
                R1.append(r1)
                R2.append(r2)
            
            R1 = np.stack(R1)
            R2 = np.stack(R2)

            # Φ
            phi_i = R1 * R2
            Phi[session] = np.nanmean(phi_i, axis=0)

            # DP
            M = R1.shape[1]
            dp = np.zeros(M)

            for i in range(N):
                phi_ii = R1[i] * R2[i]
                count = np.zeros(M)

                for j in range(N):
                    if i == j:
                        continue

                    phi_ij = R1[i] * R2[j]
                    phi_ji = R1[j] * R2[i]

                    count += (phi_ij > phi_ii)
                    count += (phi_ji > phi_ii)
                
                Pi = count / (2 * (N - 1))
                Pi = np.clip(Pi, 1e-10, None)

                dp += -np.log(Pi)
            
            DP[session] = dp
        
        return Phi, DP
    
    def run(self, repeat=1):

        if isinstance(repeat, bool) or not isinstance(repeat, int) or repeat < 1:
            raise ValueError("Edge.repeat must be a positive integer")
        # Shared groups are fixed; repeating identical pairs adds no information.
        if self.group_indices is not None:
            return self.compute_phi_dp(self.build_pairs())
        Phi_all, DP_all = [], []

        for r in range(repeat):

            pairs = self.build_pairs()
            Phi, DP = self.compute_phi_dp(pairs)

            Phi_all.append(Phi)
            DP_all.append(DP)
        
        Phi_mean, DP_mean = {}, {}

        sessions = Phi_all[0].keys()

        for s in sessions:
            Phi_mean[s] = np.mean([p[s] for p in Phi_all], axis=0)
            DP_mean[s] = np.mean([d[s] for d in DP_all], axis=0)
        
        return Phi_mean, DP_mean
