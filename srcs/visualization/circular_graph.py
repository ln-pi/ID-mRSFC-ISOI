"""Circular connectivity graph visualization.

Adapted from the MATLAB circularGraph project:
https://github.com/paul-kassebaum-mathworks/circularGraph
"""

import matplotlib.pyplot as plt
import numpy as np

from matplotlib import colormaps
from matplotlib.collections import LineCollection
from matplotlib.colors import LinearSegmentedColormap, to_rgb
from matplotlib.patches import FancyArrowPatch

from srcs.utils.logger import setup_logger


logger = setup_logger()


class Node:
    def __init__(self, idx, pos, label, color):
        self.idx = idx
        self.pos = pos
        self.label = label
        self.color = color

        self.artist = None
        self.text = None
        self.active = True

    def draw(self, ax):
        x, y = self.pos

        self.artist, = ax.plot(
            x, y,
            marker="o",
            markersize=8,
            color=self.color,
            linestyle="None",
            zorder=3,
            picker=5
        )

        r = 1.1
        lx, ly = x * r, y * r

        angle = np.degrees(np.arctan2(y, x))
        rotation = angle
        ha = "left"

        if angle > 90 or angle < -90:
            rotation += 180
            ha = "right"

        self.text = ax.text(
            lx, ly,
            self.label,
            fontsize=9,
            ha=ha,
            va="center",
            rotation=rotation,
            rotation_mode="anchor"
        )

    def toggle(self):
        self.active = not self.active

        if self.active:
            self.artist.set_marker("o")
            self.artist.set_color(self.color)
        else:
            self.artist.set_marker("x")
            self.artist.set_color("gray")


class circularGraph:
    """
    Generate a circular dot-edge graph based on the input adjacency matrix A.
    Nodes are arranged in a circle, and edges are drawn as arcs connecting the nodes.

    Parameters:
    - A: 2D array-like adjacency matrix (N x N)
    - labels: List of node labels (default: ["0", "1", ..., "N - 1"])
    - node_color_mode: "gradient" (default) or "fixed"
    - node_color_spec: List of node names for "fixed" mode (e.g., ["group1_L", "group1_R", "group2_L", "group2_R", ...])
    - edge_color_mode: "source" (default), "gradient", or "fixed"
    - edge_fixed_color: Color for edges in "fixed" mode (default: "black")
    - edge_width_mode: "weight" (default) or "fixed"
    - threshold: Minimum weight involved in visualization
    - gap: Distance between left&right half of the nodes circle
    - cross_only: Block cross-left-right edges (default: True)
    - show: Show graph after generated (default: False)
    - save_path: Path to save the graph (default: None, noSave)
    - save_mat: Permission for generating an additonal matrix-version plot
    """

    def __init__(
            self, A, labels=None,
            node_color_mode="gradient",
            node_color_spec=None,
            edge_color_mode="source",
            edge_fixed_color="black",
            edge_width_mode="weight",
            threshold=0.0,
            gap=0.2,
            cross_only=True,
            show=False,
            save_path=None,
            save_mat=False
    ):
        self.A = np.array(A)
        self.N = self.A.shape[0]

        self.labels = labels if labels else [str(i) for i in range(self.N)]

        self.node_color_mode = node_color_mode
        self.node_color_spec = node_color_spec

        self.edge_color_mode = edge_color_mode
        self.edge_fixed_color = edge_fixed_color
        self.edge_width_mode = edge_width_mode

        self.threshold = threshold
        self.gap = gap
        self.cross_only = cross_only

        self.nodes = []
        self.edges = []

        self.node_colors = self._build_node_colors()

        self.fig, self.ax = plt.subplots(figsize=(8, 8))
        self.ax.axis("off")

        self._create_nodes()
        self._draw_edges()
        self._draw_nodes()

        self.fig.canvas.mpl_connect("button_press_event", self._on_click)

        if show:
            plt.show()

        if save_path is not None:
            if isinstance(save_path, str):
                plt.savefig(save_path + "_cir.png")

            if save_mat:

                self._draw_mat(save_path + "_mat.png")

    # ====================
    # Node color system
    # ====================
    def _build_node_colors(self):
        half = self.N // 2

        if self.node_color_mode == "gradient":
            cmap = colormaps["viridis"].resampled(half)
            left_color = [cmap(i) for i in range(half)]
            right_color = left_color[::-1]
            return left_color + right_color

        elif self.node_color_mode == "fixed":
            if self.node_color_spec is None:
                logger.raise_(ValueError("'node_color_spec' required!"))

            left_color = [None] * half
            cmap = colormaps["tab10"]

            left_nodes = [name for name in self.node_color_spec if name.endswith("_L")]

            prefix_groups = {}
            for idx, name in enumerate(left_nodes):
                prefix = name.split("-")[0]
                if prefix not in prefix_groups:
                    prefix_groups[prefix] = []
                prefix_groups[prefix].append(idx)

            gid = 0
            for prefix, indices in prefix_groups.items():
                color = cmap(gid % 10)

                for idx in indices:
                    left_color[idx] = color

                gid += 1

            right_color = left_color[::-1]

            return left_color + right_color

    # ====================
    # Node layout
    # ====================
    def _create_nodes(self):
        half = self.N // 2

        theta_left = np.linspace(
            np.pi / 2 + self.gap,
            3 * np.pi / 2 - self.gap,
            half,
            endpoint=False
        )

        theta_right = np.pi - theta_left[::-1]

        theta = np.concatenate([theta_left, theta_right])

        for i in range(self.N):
            x, y = np.cos(theta[i]), np.sin(theta[i])
            node = Node(i, (x, y), self.labels[i], self.node_colors[i])
            self.nodes.append(node)

    # ====================
    # Edge drawing
    # ====================
    def _draw_edges(self):
        for i in range(self.N):
            for j in range(i + 1, self.N):

                if self.cross_only:
                    if (i < self.N / 2 and j < self.N / 2) or (i >= self.N / 2 and j >= self.N / 2):
                        continue

                w = self.A[i, j]
                if abs(w) < self.threshold or w == 0:
                    continue

                p1 = np.array(self.nodes[i].pos)
                p2 = np.array(self.nodes[j].pos)

                rad = self._compute_rad(i, j)

                lw = abs(w) * 2 if self.edge_width_mode == "weight" else 1.5

                if self.edge_color_mode == "gradient":
                    self._draw_gradient_edge(
                        p1, p2, rad,
                        self.node_colors[i], self.node_colors[j],
                        lw
                    )
                else:
                    if self.edge_color_mode == "source":
                        color = self.node_colors[i]
                    else:
                        color = self.edge_fixed_color

                    patch = FancyArrowPatch(
                        tuple(p1), tuple(p2),
                        connectionstyle=f"arc3,rad={rad}",
                        arrowstyle="-",
                        linewidth=lw,
                        color=color,
                        alpha=0.7,
                        zorder=1
                    )
                    self.ax.add_patch(patch)
                    self.edges.append((i, j, patch))

    def _compute_rad(self, i, j):
        p1 = np.array(self.nodes[i].pos)
        p2 = np.array(self.nodes[j].pos)

        chord = p2 - p1
        mid = (p1 + p2) / 2
        to_center = -mid
        cross = chord[0] * to_center[1] - chord[1] * to_center[0]

        diff = j - i
        if diff < 0:
            diff += self.N
        dist = min(diff, self.N - diff)
        norm = dist / (self.N / 2)

        max_rad = 0.6
        base = max_rad * (1 - norm)

        sign = np.sign(cross)
        if sign == 0:
            sign = 1

        return -sign * base

    def _draw_gradient_edge(self, p1, p2, rad, c1, c2, lw):
        n = 100

        xs = np.linspace(p1[0], p2[0], n)
        ys = np.linspace(p1[1], p2[1], n)

        points = np.array([xs, ys]).T.reshape(-1, 1, 2)
        segs = np.concatenate([points[:-1], points[1:]], axis=1)

        c1 = np.array(c1)
        c2 = np.array(c2)

        if len(c1) == 3:
            c1 = np.append(c1, 1.0)
        if len(c2) == 3:
            c2 = np.append(c2, 1.0)

        colors = [(1 - t) * c1 + t * c2 for t in np.linspace(0, 1, n-1)]

        lc = LineCollection(segs, colors=colors, linewidth=lw, alpha=0.8)
        lc._base_colors = np.array(colors)

        self.ax.add_collection(lc)

        return lc

    # ====================
    # Draw nodes
    # ====================
    def _draw_nodes(self):
        for node in self.nodes:
            node.draw(self.ax)

    # ====================
    # Interaction
    # ====================
    def _on_click(self, event):
        if event.inaxes != self.ax:
            return

        for node in self.nodes:
            x, y = node.pos
            if np.hypot(event.xdata - x, event.ydata - y) < 0.05:
                node.toggle()
                self._update_edges()
                self.fig.canvas.draw_idle()
                break

    def _update_edges(self):
        for i, j, obj in self.edges:

            active = self.nodes[i].active and self.nodes[j].active

            if hasattr(obj, "_base_colors"):

                if active:
                    new_colors = obj._base_colors.copy()
                    new_colors[:, 3] = 0.8
                else:
                    n = len(obj._base_colors)
                    new_colors = np.tile([0.7, 0.7, 0.7, 0.05], (n, 1))

                obj.set_color(new_colors)

            else:
                obj.set_alpha(0.7 if active else 0.05)

    # ====================
    # (Optional) Matrix of circularGraph
    #=====================
    def _draw_mat(self, save_path):
        N = self.N
        group_spec = self.node_color_spec

        prefixes = [s.split("-")[0] for s in group_spec]

        unique_groups = sorted(set(prefixes))
        G = len(unique_groups)

        prefix2id = {p: i for i, p in enumerate(unique_groups)}

        node_to_group = np.array([prefix2id[p] for p in prefixes], dtype=int)

        iu = np.triu_indices(N, 1)
        values = self.A[iu]
        mask = values > self.threshold

        edge_i = iu[0][mask]
        edge_j = iu[1][mask]

        density = np.zeros((G, G))

        for i, j in zip(edge_i, edge_j):
            gi = node_to_group[i]
            gj = node_to_group[j]
            density[gi, gj] += 1
            density[gj, gi] += 1

        total_edges = len(edge_i)
        proportion = density / (2 * total_edges) if total_edges > 0 else density

        tri_mask = np.triu(np.ones_like(proportion, bool), k=1)
        proportion_plot = proportion.copy()
        proportion_plot[tri_mask] = np.nan

        try:
            color = to_rgb(self.edge_fixed_color)
        except ValueError:
            raise ValueError(f"Invalid color name: {self.edge_fixed_color}")

        cmap = LinearSegmentedColormap.from_list("custom", [(1, 1, 1), color])
        cmap.set_bad(color="white")

        fig, ax = plt.subplots(figsize=(6, 6))

        for spine in ax.spines.values():
            spine.set_visible(False)

        im = ax.imshow(proportion_plot, cmap=cmap, interpolation="none")

        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label("Edge proportion")

        ax.set_title("")
        ax.set_xticks(np.arange(G))
        ax.set_yticks(np.arange(G))
        ax.set_xticklabels(unique_groups)
        ax.set_yticklabels(unique_groups)

        ax.grid(False)
        for i in range(G):
            ax.hlines(y=i - 0.5, xmin=-0.5, xmax=i + 0.5, color="black", linewidth=0.5)
            ax.vlines(x=i + 0.5, ymin=i - 0.5, ymax=G - 0.5, color="black", linewidth=0.5)
        ax.hlines(y=G - 0.5, xmin=-0.5, xmax=G - 0.5, color="black", linewidth=0.5)
        ax.vlines(x=-0.5, ymin=-0.5, ymax=G - 0.5, color="black", linewidth=0.5)

        plt.savefig(save_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
