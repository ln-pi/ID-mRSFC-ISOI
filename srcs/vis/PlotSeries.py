import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns
import sys

from dataclasses import dataclass, field
from datetime import datetime
from matplotlib.patches import Patch
from scipy.stats import wilcoxon, ttest_rel
from typing import Any

from srcs.utils.logger import setup_logger


logger = setup_logger()


@dataclass
class PlotTask:
    data: Any
    plot_type: str = "heatmap"
    name: str | None = None
    config: dict[str, Any] = field(default_factory=dict)
    slice_dim: Any = None
    fmt: str | None = None


def roi_tick_config(names, fontsize_small=10, fontsize_scale=300):
    labels = [n.split("-", 1)[1] if "-" in n else n for n in names]
    fontsize = fontsize_scale / len(names) if len(names) >= 30 else fontsize_small
    ticks = np.arange(len(names)) + 0.5
    return {
        "xticks": {
            "ticks": ticks,
            "labels": labels,
            "rotation": 90,
            "fontsize": fontsize,
        },
        "yticks": {
            "ticks": ticks,
            "labels": labels,
            "rotation": 0,
            "fontsize": fontsize,
        },
    }


def fc_heatmap_config(names, **overrides):
    overrides = dict(overrides)
    fontsize_small = overrides.pop("label_fontsize_small", 10)
    fontsize_scale = overrides.pop("label_fontsize_scale", 300)
    config = {
        "center": 0,
        "vmin": -1,
        "vmax": 1,
        "square": True,
        "cmap": "coolwarm",
        "cbar": True,
    }
    config.update(roi_tick_config(names, fontsize_small=fontsize_small, fontsize_scale=fontsize_scale))
    config.update(overrides)
    return config


def save_fc_heatmap(data, names, name, save_path, figsize=(8, 8), dpi=180, show=False, save=True, log=False, config=None):
    config = dict(config or {})
    figsize = tuple(config.pop("figsize", figsize))
    dpi = config.pop("dpi", dpi)
    return (
        PlotSeries(figsize=figsize, dpi=dpi, save_path=save_path)
        .fc_heatmap(data, names=names, name=name, **config)
        .render(show=show, save=save, log=log)
    )


class PlotSeries:

    name = "A complex plotter, e.g., PlotSeries.add_plot().render()"

    def __init__(self, figsize=(10, 5), dpi=120, default_fmt="png", layout="auto", save_path=None, auto_close=True):
        """
        Params:
        - layout(str): "auto", "h": horizontal, "v": vertical
        - save_path(str or None): default save path
        - auto_close(bool): close figure after plot
        """
        self.figsize = figsize
        self.dpi = dpi
        self.default_fmt = default_fmt
        self.layout = layout
        self.save_path = save_path
        self.auto_close = auto_close
        self.plot_tasks = []
        self._fig = None
        self._axes = None
    
    # ====================
    # API
    # ====================
    def add_plot(self, data, slice_dim=None, plot_type="heatmap", config=None, name=None, fmt=None):
        self.plot_tasks.append(PlotTask(
            name=name,
            data=data,
            slice_dim=slice_dim,
            plot_type=plot_type,
            config=config or {},
            fmt=fmt,
        ))
        return self

    def add(self, data, plot_type="heatmap", name=None, slice_dim=None, fmt=None, **config):
        return self.add_plot(data, slice_dim=slice_dim, plot_type=plot_type, config=config, name=name, fmt=fmt)

    def heatmap(self, data, name=None, slice_dim=None, fmt=None, **config):
        return self.add(data, plot_type="heatmap", name=name, slice_dim=slice_dim, fmt=fmt, **config)

    def fc_heatmap(self, data, names, name=None, slice_dim=None, fmt=None, **overrides):
        return self.heatmap(
            data,
            name=name,
            slice_dim=slice_dim,
            fmt=fmt,
            **fc_heatmap_config(names, **overrides),
        )

    def origin(self, data, name=None, slice_dim=None, fmt=None, **config):
        return self.add(data, plot_type="origin", name=name, slice_dim=slice_dim, fmt=fmt, **config)

    def box(self, data, name=None, slice_dim=None, fmt=None, **config):
        return self.add(data, plot_type="box", name=name, slice_dim=slice_dim, fmt=fmt, **config)

    def hist(self, data, name=None, slice_dim=None, fmt=None, **config):
        return self.add(data, plot_type="hist", name=name, slice_dim=slice_dim, fmt=fmt, **config)

    def mean_error(self, data, name=None, slice_dim=None, fmt=None, **config):
        return self.add(data, plot_type="mean_error", name=name, slice_dim=slice_dim, fmt=fmt, **config)

    def paired_vertical(self, data, name=None, slice_dim=None, fmt=None, **config):
        return self.add(data, plot_type="paired_vertical", name=name, slice_dim=slice_dim, fmt=fmt, **config)

    def raw_save(self, data, name=None, slice_dim=None, fmt=None, **config):
        return self.add(data, plot_type="raw_save", name=name, slice_dim=slice_dim, fmt=fmt, **config)
    
    def render(self, show=False, save=True, log=True, dpi=300):
        if len(self.plot_tasks) == 0:
            logger.warning("No plot tasks to render.")
            return None, None
        
        raw_tasks = [t for t in self.plot_tasks if t.plot_type == "raw_save"]
        normal_tasks = [t for t in self.plot_tasks if t.plot_type != "raw_save"]

        fig, axes = None, None
        if len(normal_tasks) != 0:
            n = len(normal_tasks)
            nrows, ncols = self._compute_layout(n)

            fig, axes = plt.subplots(
                nrows, ncols,
                figsize=self.figsize,
                dpi=self.dpi,
            )
            axes = np.array(axes).reshape(-1)

            for ax, task in zip(axes, normal_tasks):
                data = self._slice_data(task.data, task.slice_dim)
                plot_type = task.plot_type
                cfg = task.config

                method = getattr(self, f"_plot_{plot_type}", None)
                if method is None:
                    logger.error(f"Unsupported plot type: {plot_type}")
                    continue
                method(ax, data, cfg)
            
            for ax in axes[len(normal_tasks):]:
                ax.axis("off")
            
            plt.tight_layout()
            self._fig = fig
            self._axes = axes

        if save:
            base_dir = self.save_path or "figs"

            if len(raw_tasks) != 0:
                for task in raw_tasks:
                    self._plot_raw_save(task, base_dir=base_dir, log=log, dpi=dpi)
            
            if len(normal_tasks) == 1:
                name = normal_tasks[0].name or "figure"
                base_path = os.path.join(base_dir, name)
                self.save(base_path, task_fmt=normal_tasks[0].fmt, log=log, dpi=dpi)
            elif len(normal_tasks) > 1:
                self._save_by_name(normal_tasks, base_dir=base_dir, log=log, dpi=dpi)

        if show:
            plt.show()

        if self.auto_close:
            plt.close(fig)
            plt.close("all")

        self.plot_tasks.clear()

        return fig, axes
    
    # ====================
    # Method
    # ====================
    def save(self, path, task_fmt=None, log=True, dpi=300):
        if self._fig is None:
            logger.error("No figure to save!")
            raise RuntimeError("No figure to save!")
        
        fmt = self._resolve_format(task_fmt, path)
        self._check_fmt_compatibility(fmt, mode="figure")

        root, ext = os.path.splitext(path)
        if not ext:
            path = f"{root}.{fmt}"

        folder = os.path.dirname(path)
        if folder and not os.path.exists(folder):
            os.makedirs(folder)
        
        self._fig.savefig(path, dpi=dpi, bbox_inches="tight", format=fmt)

        if log:
            logger.info(f"Figure saved to {path} in format {fmt}")
    
    def _resolve_format(self, fmt=None, path=None):
        if path is not None:
            _, ext = os.path.splitext(path)
            if ext:
                return ext.lstrip(".")
        if fmt is not None:
            return fmt
        return self.default_fmt
    
    def _check_fmt_compatibility(self, fmt, mode="figure"):
        fig_fmt = {"png", "pdf", "svg", "jpg", "jpeg", "tif", "tiff"}
        img_fmt = {"png", "jpg", "jpeg", "tif", "tiff"}
        if mode == "figure" and fmt not in fig_fmt:
            logger.error(f"Unsupported figure format: {fmt}")
            raise ValueError(f"Unsupported figure format: {fmt}")
        if mode == "image" and fmt not in img_fmt:
            logger.error(f"Unsupported image format: {fmt}")
            raise ValueError(f"Unsupported image format: {fmt}")
        
    def _auto_filename(self, fmt=None):
        fmt = fmt or self.default_fmt
        time_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        return f"figure_{time_stamp}.{fmt}"
    
    def _save_by_name(self, tasks, base_dir="figs", log=True, dpi=300):
        if self._fig is None:
            logger.error("No figure to save!")
            raise RuntimeError("No figure to save!")
        
        os.makedirs(base_dir, exist_ok=True)
        for idx, (ax, task) in enumerate(zip(self._axes, tasks)):
            name = task.name or f"plot_{idx}"
            fmt = self._resolve_format(task.fmt)
            self._check_fmt_compatibility(fmt, mode="figure")
            path = os.path.join(base_dir, f"{name}.{fmt}")
            fig = ax.figure
            fig.savefig(path, dpi=dpi, bbox_inches=ax.get_tightbbox(fig.canvas.get_renderer()), format=fmt)
            if log:
                logger.info(f"Saved subplot: {path}")

    def _slice_data(self, data, slice_dim):
        """
        slice_dim:
            None
            tuple(slice, slice)
            dict:
                {
                    "rows": ...,
                    "cols": ...,
                    "mode": "iloc" | "loc"
                }
        """
        # ----- None -----
        if slice_dim is None:
            return data
        
        # ----- NumPy -----
        if isinstance(data, np.ndarray):
            try:
                return data[slice_dim]
            except Exception as e:
                raise RuntimeError(
                    f"NumPy slicing failed: shape={data.shape}, slice_dim={slice_dim}"
                ) from e
        
        # ----- Pandas -----
        if isinstance(data, (pd.DataFrame, pd.Series)):
            if isinstance(slice_dim, tuple):
                try:
                    return data.iloc[slice_dim]
                except Exception as e:
                    raise RuntimeError(
                        f"Pandas iloc slicing failed: slice_dim={slice_dim}"
                    ) from e
            if isinstance(slice_dim, dict):
                rows = slice_dim.get("rows", slice(None))
                cols = slice_dim.get("cols", slice(None))
                mode = slice_dim.get("mode", "iloc")
                try:
                    if mode == "iloc":
                        return data.iloc[rows, cols]
                    elif mode == "loc":
                        return data.loc[rows,cols]
                    else:
                        raise ValueError(f"Unknown slicing mode: {mode}")
                except Exception as e:
                    raise RuntimeError(
                        f"Pandas {mode} slicing failed: rows={rows}, cols={cols}"
                    ) from e
            raise TypeError(f"Unsupported slice_dim type: {type(slice_dim)}")
        
        # ----- fallback -----
        raise TypeError(f"Unsupported data type for slicing: {type(data)}")

    def _compute_layout(self, n):
        if self.layout == "h":
            return 1, n
        elif self.layout == "v":
            return n, 1
        elif self.layout == "auto":
            ncols = math.ceil(math.sqrt(n))
            nrows = math.ceil(n / ncols)
            return nrows, ncols
        else:
            logger.error("Not a valid type! Check layout type plz!")
            sys.exit()
    
    @staticmethod
    def _get_star(p):
        if p < 0.001:
            return "***"
        elif p < 0.01:
            return "**"
        elif p < 0.05:
            return "*"
        else:
            return "n.s."
    
    @staticmethod
    def _compute_mean_error(data, error_type="sem"):
        means, errors = [], []

        for d in data:
            d = np.asarray(d)
            
            if len(d) == 0:
                means.append(np.nan)
                errors.append(0)
                continue

            m = np.mean(d)

            if error_type == "std":
                e = np.std(d)
            elif error_type == "sem":
                e = np.std(d) / np.sqrt(len(d))
            elif error_type == "ci":
                e = 1.96 * np.std(d) / np.sqrt(len(d))
            else:
                logger.raise_(ValueError(f"Unsupport error type: {error_type}"))
            
            means.append(m)
            errors.append(e)
        
        return np.array(means), np.array(errors)
    
    # ====================
    # Plot
    # ====================
    def _plot_raw_save(self, task, base_dir, log=True, dpi=300):
        # skip plotting
        data = self._slice_data(task.data, task.slice_dim)
        cfg = task.config
        name = task.name or "raw"
        fmt = self._resolve_format(task.fmt)
        self._check_fmt_compatibility(fmt, mode="image")
        os.makedirs(base_dir, exist_ok=True)
        path = os.path.join(base_dir, f"{name}.{fmt}")
        plt.imsave(
            path,
            data,
            cmap=cfg.get("cmap", "gray"),
            vmin=cfg.get("vmin", None),
            vmax=cfg.get("vmax", None)
        )
        if log:
            logger.info(f"Raw image saved to {path} in format {fmt}")

    def _plot_origin(self, ax, data, cfg):
        im = ax.imshow(
            data,
            cmap=cfg.get("cmap", "gray"),
            vmin=cfg.get("vmin", None),
            vmax=cfg.get("vmax", None),
            origin=cfg.get("origin", "upper"),
            interpolation=cfg.get("interpolation", "nearest")
        )
        if cfg.get("cbar", False):
            plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        
        xticks_cfg = cfg.get("xticks")
        yticks_cfg = cfg.get("yticks")

        if xticks_cfg == "index":
            if isinstance(data, pd.DataFrame):
                labels = list(data.columns)
            else:
                labels = list(range(data.shape[1]))
            ax.set_xticks(np.arange(len(labels)))
            ax.set_xticklabels(labels, rotation=90, fontsize=6)
        elif isinstance(xticks_cfg, dict):
            ax.set_xticks(xticks_cfg.get("ticks"))
            ax.set_xticklabels(
                xticks_cfg.get("labels", []),
                rotation=xticks_cfg.get("rotation", 0),
                fontsize=xticks_cfg.get("fontsize", 10),
            )
        
        if yticks_cfg == "index":
            if isinstance(data, pd.DataFrame):
                labels = list(data.index)
            else:
                labels = list(range(data.shape[0]))
            ax.set_yticks(np.arange(len(labels)))
            ax.set_yticklabels(labels, rotation=0, fontsize=6)
        elif isinstance(yticks_cfg, dict):
            ax.set_yticks(yticks_cfg.get("ticks"))
            ax.set_yticklabels(
                yticks_cfg.get("labels", []),
                rotation=yticks_cfg.get("rotation", 0),
                fontsize=yticks_cfg.get("fontsize", 10),
            )
        
        if  cfg.get("axis", True) is False:
            ax.axis("off")
        if "title" in cfg:
            ax.set_title(cfg["title"])
        if "xlabel" in cfg:
            ax.set_xlabel(cfg["xlabel"])
        if "ylabel" in cfg:
            ax.set_ylabel(cfg["ylabel"])
    
    def _plot_uint8(self, ax, data, cfg):
        vmin = cfg.get("vmin", np.nanmin(data))
        vmax = cfg.get("vmax", np.nanmax(data))
        if vmax <= vmin:
            logger.warning("Invalid vmin/vmax for uint8 normalization")
            norm = np.zeros_like(data, dtype=np.uint8)
        else:
            norm = (data - vmin) / (vmax - vmin)
            norm = np.clip(norm, 0, 1)
            norm = (norm * 255).astype(np.uint8)
        
        # plot origin
        self._plot_origin(ax, norm, cfg)
    
    def _plot_heatmap(self, ax, data, cfg):
        sns.heatmap(
            data=data,
            ax=ax,
            center=cfg.get("center"),
            vmin=cfg.get("vmin"),
            vmax=cfg.get("vmax"),
            square=cfg.get("square", False),
            cmap=cfg.get("cmap", "gray"),
            cbar=cfg.get("cbar", True),
            xticklabels=False,
            yticklabels=False,
        )

        xticks_cfg = cfg.get("xticks")
        if xticks_cfg == "index":
            assert isinstance(data, pd.DataFrame)
            ax.set_xticks(np.arange(len(data.columns)) + 0.5)
            ax.set_xticklabels(
                data.columns,
                rotation=90,
                fontsize=6,
            )
        elif isinstance(xticks_cfg, dict):
            ax.set_xticks(xticks_cfg.get("ticks"))
            ax.set_xticklabels(
                xticks_cfg.get("labels", []),
                rotation=xticks_cfg.get("rotation", 0),
                fontsize=xticks_cfg.get("fontsize", 10),
            )

        yticks_cfg = cfg.get("yticks")
        if yticks_cfg == "index":
            assert isinstance(data, pd.DataFrame)
            ax.set_yticks(np.arange(len(data.index)) + 0.5)
            ax.set_yticklabels(
                data.index,
                rotation=0,
                fontsize=6,
            )
        elif yticks_cfg:
            ax.set_yticks(yticks_cfg.get("ticks"))
            ax.set_yticklabels(
                yticks_cfg.get("labels", []),
                rotation=yticks_cfg.get("rotation", 0),
                fontsize=yticks_cfg.get("fontsize", 10),
            )

        if "title" in cfg:
            ax.set_title(cfg["title"])
    
    def _plot_violin(self, ax, data, cfg):
        sns.violinplot(
            y=data,
            ax=ax,
            inner=cfg.get("inner", "box"),
            cut=cfg.get("cut", 0),
            linewidth=cfg.get("linewidth", 1)
        )
        ax.set_xlabel(cfg.get("xlabel", ""))
        ax.set_ylabel(cfg.get("ylabel", "Value"))
        if "title" in cfg:
            ax.set_title(cfg["title"])

    def _plot_paired_line(self, ax, data, cfg):
        x1 = np.array(data[0])
        x2 = np.array(data[1])

        avg_x1 = np.mean(x1)
        avg_x2 = np.mean(x2)

        for p1, p2 in zip(x1, x2):
            ax.plot(
                [p1, p2], [0, 1],
                color=cfg.get("linecolor", "gray"),
                alpha=cfg.get("alpha", 0.6),
                linewidth=cfg.get("linewidth", 2),
            )
        
        ax.plot(
            [avg_x1, avg_x2], [0, 1],
            color=cfg.get("meancolor", "red"),
            linewidth=cfg.get("meanwidth", 3),
        )

        ax.set_yticks([0, 1])
        ax.set_yticklabels(cfg.get("labels", ["x1", "x2"]))
        ax.set_xlabel(cfg.get("xlabel", "Value"))

        if cfg.get("grid", True):
            ax.grid(True, alpha=0.3)
        
        if "title" in cfg:
            ax.set_title(cfg["title"])
    
    def _plot_paired_vertical(self, ax, data, cfg):
        x1 = np.array(data[0])
        x2 = np.array(data[1])

        pos1, pos2 = 0, 1

        for p1, p2 in zip(x1, x2):
            ax.plot(
                [pos1, pos2], [p1, p2],
                color=cfg.get("linecolor", "gray"),
                alpha=cfg.get("alpha", 0.5),
                linewidth=cfg.get("linewidth", 1.5)
            )
        ax.plot(
            [pos1, pos2], [np.mean(x1), np.mean(x2)],
            color=cfg.get("meancolor", "red"),
            linewidth=cfg.get("meanwidth", 3),
            zorder=4
        )

        ax.scatter(np.full_like(x1, pos1), x1, color=cfg.get("linecolor", "gray"), s=30, zorder=3)
        ax.scatter(np.full_like(x2, pos2), x2, color=cfg.get("linecolor", "gray"), s=30, zorder=3)

        ax.scatter(pos1, np.mean(x1), color=cfg.get("meancolor", "red"), s=120, zorder=5)
        ax.scatter(pos2, np.mean(x2), color=cfg.get("meancolor", "red"), s=120, zorder=5)

        stat, p = wilcoxon(x1, x2)
        star = self._get_star(p)

        y_max = max(np.max(x1), np.max(x2))
        y_bar = y_max + 0.1

        ax.plot([pos1, pos1, pos2, pos2],
                [y_bar, y_bar+0.05, y_bar+0.05, y_bar],
                color="black")

        ax.text((pos1+pos2)/2, y_bar+0.05, star,
                ha='center', va='bottom', fontsize=14)

        ax.set_xticks([pos1, pos2])
        ax.set_xticklabels(cfg.get("xticklabels", ["Group", "Individual"]))
        ax.set_ylabel(cfg.get("ylabel", "Network Similarity"))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    def _plot_box(self, ax, data, cfg):
        base_color = [
            '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
            '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
        ]
        is_2_group = len(data) == 2 and isinstance(data[0][0], (list, np.ndarray))
        
        if is_2_group:
            Group, Indiv = data
            all_data = Group + Indiv
            n_group = len(Group)
            colors = [base_color[i % len(base_color)] for i in range(n_group)]
            colors += colors
        else:
            all_data = data
            n_group = len(data)
            colors = [base_color[i % len(base_color)] for i in range(len(all_data))]
        
        gap = cfg.get("group_gap", 2)
        if is_2_group:
            pos1 = np.arange(n_group)
            pos2 = np.arange(n_group, len(all_data)) + gap
            positions = np.concatenate([pos1, pos2])
        else:
            positions = np.arange(len(all_data))
        
        box = ax.boxplot(
            all_data,
            positions=positions,
            widths=cfg.get("box_width", 0.1),
            patch_artist=True,
            showfliers=cfg.get("showfliers", False)
        )

        for i, (b, med) in enumerate(zip(box["boxes"], box["medians"])):
            color = colors[i]
            b.set(facecolor="none", edgecolor=color, linewidth=1.2)
            med.set(color=color)

        for w, c in zip(box["whiskers"], np.repeat(colors, 2)):
            w.set(color=c, linewidth=1)

        for cap, c in zip(box["caps"], np.repeat(colors, 2)):
            cap.set(color=c)

        all_values = np.concatenate(
            [np.asarray(d).ravel() for d in all_data]
        )

        y_min = np.nanmin(all_values)
        y_max = np.nanmax(all_values)

        y_range = y_max - y_min

        base_h = y_max + y_range * 0.05
        step_h = y_range * 0.08
        tick_h = y_range * 0.02
        
        def add_stats(pairs, y):
            results = []
            for i, j in pairs:
                _, p = ttest_rel(all_data[i], all_data[j])
                star = self._get_star(p)

                if cfg.get("hide_ns", False) and star == "n.s.":
                    continue

                span = j - i
                results.append((i, j, star, span))
            results.sort(key=lambda x: x[3])

            levels = []

            def overlap(a, b):
                return not (a[1] <= b[0] or b[1] <= a[0])
            
            for i, j, star, span in results:
                placed = False

                for level in levels:
                    if all(not overlap((i, j), (ii, jj)) for ii, jj, _, _ in level):
                        level.append((i, j, star, span))
                        placed = True
                        break
                
                if not placed:
                    levels.append([(i, j, star, span)])
            
            h = y
            for level in levels:
                for i, j, star, _ in level:
                    x1, x2 = positions[i]+0.05, positions[j]-0.05

                    color = "gray" if star == "n.s." else "black"

                    ax.plot([x1, x1], [h, h - tick_h], lw=0.8, c=color)
                    ax.plot([x2, x2], [h, h - tick_h], lw=0.8, c=color)
                    ax.plot([x1, x2], [h, h], lw=0.8, c=color)

                    ax.text((x1 + x2) / 2, h, star, ha="center", va="bottom", fontsize=7, color=color)

                h += step_h

            return h    

        if is_2_group:
            pairs1 = [(i, j) for i in range(n_group) for j in range(i + 1, n_group)]
            pairs2 = [(i, j) for i in range(n_group, len(all_data)) for j in range(i + 1, len(all_data))]
            h = add_stats(pairs1, base_h)
            current_h = add_stats(pairs2, h)
        else:
            pairs = [(i, j) for i in range(len(all_data)) for j in range(i + 1, len(all_data))]
            current_h = add_stats(pairs, base_h)
        
        ax.set_xticks(positions)
        sub_labels = cfg.get("xticklabel", [str(i + 1) for i in range(n_group)])
        if is_2_group and len(sub_labels) == n_group:
            sub_labels += sub_labels
        
        if cfg.get("show_legend", True):
            label_color = {}
            for lab, color in zip(sub_labels, colors):
                if lab not in label_color:
                    label_color[lab] = color
            handles = [
                Patch(facecolor="none", edgecolor=color, linewidth=1.5, label=lab)
                for lab, color in label_color.items()
            ]
            ax.legend(
                handles=handles,
                frameon=True,
                fontsize=8,
                loc=cfg.get("legend_loc", "upper left"),
                bbox_to_anchor=cfg.get("legend_box", (0.02, 0.98)),
                borderaxespad=0
            )

        if cfg.get("show_sub_label", False):
            ax.set_xticklabels(sub_labels, fontsize=7)
        else:
            ax.set_xticklabels([""] * len(positions))
        
        if is_2_group:
            mid1 = np.mean(positions[:n_group])
            mid2 = np.mean(positions[n_group:])

            ax.text(mid1, -0.08, "Group", ha="center", va="top", transform=ax.get_xaxis_transform(), fontsize=9)
            ax.text(mid2, -0.08, "Individual", ha="center", va="top", transform=ax.get_xaxis_transform(), fontsize=9)

        ax.set_ylim(
            y_min - y_range * 0.1,
            current_h + y_range * 0.15
        )
        if cfg.get("title", None) is not None:
            ax.set_title(cfg.get("title", None))
        ax.grid(cfg.get("grid", False))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    def _plot_hist(self, ax, data, cfg):

        # histogram config
        bins = cfg.get("bins", 20)
        hist_range = cfg.get("range", (0, 1))

        facecolor = cfg.get("facecolor", "lightgray")
        edgecolor = cfg.get("edgecolor", "dimgray")
        linewidth = cfg.get("linewidth", 1.5)

        # manual arrow
        arrow_x = cfg.get("arrow_x", None)
        arrow_color = cfg.get("arrow_color", "red")

        # mean arrow
        show_mean_arrow = cfg.get("show_mean_arrow", True)
        mean_arrow_color = cfg.get("mean_arrow_color", "gray")

        # shared arrow style
        arrow_height_ratio = cfg.get("arrow_height_ratio", 0.22)
        arrow_linewidth = cfg.get("arrow_linewidth", 2)

        # text
        show_text = cfg.get("show_text", True)
        arrow_text_format = cfg.get("arrow_text_format", "acc={:.4f}")
        mean_text_format = cfg.get("mean_text_format", "mean={:.4f}")
        text_offset_ratio = cfg.get("text_offset_ratio", 0.03)

        # labels
        xlabel = cfg.get("xlabel", "Accuracy")
        ylabel = cfg.get("ylabel", "Count")
        title = cfg.get("title", None)

        # histogram
        counts, edges, patches = ax.hist(
            data,
            bins=bins,
            range=hist_range,
            color=facecolor,
            edgecolor=edgecolor,
            linewidth=linewidth,
        )

        y_max = counts.max()

        arrow_top = y_max * arrow_height_ratio

        # manual arrow
        if arrow_x is not None:

            ax.annotate(
                "",
                xy=(arrow_x, 0),
                xytext=(arrow_x, arrow_top),

                arrowprops=dict(
                    arrowstyle="-|>",
                    color=arrow_color,
                    lw=arrow_linewidth,
                )
            )

            if show_text:

                ax.text(
                    arrow_x,
                    arrow_top + y_max * text_offset_ratio,
                    arrow_text_format.format(arrow_x),
                    color=arrow_color,
                    ha="center",
                    va="bottom",
                    fontsize=10,
                )
        
        # mean arrow
        if show_mean_arrow:

            mean_x = np.mean(data)
            ax.annotate(
                "",
                xy=(mean_x, 0),
                xytext=(mean_x, arrow_top),
                
                arrowprops=dict(
                    arrowstyle="-|>",
                    color=mean_arrow_color,
                    lw=arrow_linewidth,
                )
            )

            if show_text:

                ax.text(
                    mean_x,
                    arrow_top + y_max * text_offset_ratio,
                    mean_text_format.format(mean_x),
                    color=mean_arrow_color,
                    ha="center",
                    va="bottom",
                    fontsize=10,
                )
        
        # labels
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if title is not None:
            ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    
    def _plot_mean_error(self, ax, data, cfg):
        error_type = cfg.get("error_type", "sem")
        use_fill = cfg.get("use_fill", False)
        labels = cfg.get("labels", ("Group", "Indiv"))
        xlabel = cfg.get("xlabel", "Index")
        ylabel = cfg.get("ylabel", "Value")
        title = cfg.get("title", None)
        show_diff = cfg.get("show_diff", True)

        if isinstance(data, tuple):
            d1, d2 = data
        else:
            d1 = data
            d2 = None
        
        x = np.arange(len(d1))
        xticklabels = cfg.get("xticklabels", [str(i + 1) for i in range(len(d1))])

        # group 1
        m1, e1 = self._compute_mean_error(d1, error_type=error_type)
        if use_fill:
            ax.plot(
                x, m1,
                marker=cfg.get("marker1", "o"),
                label=labels[0]
            )
            ax.fill_between(
                x, m1 - e1, m1 + e1,
                alpha=cfg.get("fill_alpha", 0.3)
            )
        else:
            ax.errorbar(
                x, m1, yerr=e1,
                marker=cfg.get("marker1", "o"),
                capsize=cfg.get("capsize", 4),
                label=labels[0]
            )
        
        # group 2
        if d2 is not None:
            m2, e2 = self._compute_mean_error(d2, error_type=error_type)
            if use_fill:
                ax.plot(
                    x, m2,
                    marker=cfg.get("marker2", "s"),
                    label=labels[1]
                )
                ax.fill_between(
                    x, m2 - e2, m2 + e2,
                    alpha = cfg.get("fill_alpha", 0.3)
                )
            else:
                ax.errorbar(
                    x, m2, yerr=e2,
                    marker=cfg.get("marker2", "s"),
                    capsize=cfg.get("capsize", 4),
                    label=labels[1]
                )
        
        # difference curve
        if show_diff:
            
            diff_mode = cfg.get("diff_mode", "absolute")
            if diff_mode == "absolute":
                diff = m2 - m1
            elif diff_mode == "relative":
                diff = 1 - (m1 / m2)
            elif diff_mode == "gradient":
                diff = np.gradient(1 - (m1 / m2))

            ax2 = ax.twinx()
            ax2.plot(
                x, diff,
                marker=cfg.get("diff_marker", "^"),
                linestyle=cfg.get("diff_linestyle", ":"),
                linewidth=cfg.get("diff_linewidth", 2),
                color=cfg.get("diff_color", "tab:green"),
                label=cfg.get("diff_ylabel", "Diff"),
            )
            ax2.set_ylabel(
                cfg.get("diff_ylabel", "Diff"),
                color=cfg.get("diff_color", "tab:green")
            )
            ax2.tick_params(axis="y", colors=cfg.get("diff_color", "tab:green"))

            lines1, labels1 = ax.get_legend_handles_labels()
            lines2, labels2 = ax2.get_legend_handles_labels()

            ax.legend(
                lines1 + lines2, labels1 + labels2,
                loc=cfg.get("legend_loc", "upper center"),
                bbox_to_anchor=cfg.get("legend_box", (0.5, -0.1)),
                ncol=cfg.get("legend_ncol", 3),
                frameon=cfg.get("legend_frameon", False)
            )
        else:
            ax.legend(frameon=False)
        
        # labels
        ax.set_xticks(x)
        ax.set_xticklabels(
            xticklabels,
            rotation=cfg.get("xtick_rotation", 0),
            fontsize=cfg.get("xtick_fontsize", 10)
        )
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if title is not None:
            ax.set_title(title)
        ax.grid(cfg.get("grid", False))
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
