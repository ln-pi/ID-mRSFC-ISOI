"""
These functions will align frames captured to brain atlas.

User need to provide Bregma. and Lambda. location handfully.
"""

import copy
import cv2
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm

from srcs.atlas.Atlas_Mask import Mask_Generator
from srcs.utils.logger import setup_logger
from srcs.io.readers import load_config
from srcs.utils.utils import parse_args


args = parse_args()
config = load_config(args.config)
if config is None:
    raise ValueError("Config is empty (None). Exiting.")
logger = setup_logger(config["paths"]["log"])


class Match_Pattern(Mask_Generator):

    name = "Align Stack to Pattern"

    def __init__(self, avr, bin=4, lbd=4.2):
        super().__init__(bin)

        self.height, self.width = self.idmap.shape
        logger.info(f"After Align Resolution: {self.height} x {self.width}")

        self.avr = ((avr / avr.max()) * 255).astype("u1")
        self.lbd = lbd * 100 / bin # pix distance between lambda and bregma: 4.2mm
        self.pad_num = int(800 / bin) # number of pad used for frame rotation

        self.idmap_sym = copy.deepcopy(self.idmap)
        self.idmap_sym[self.idmap_sym > 32] -= 31

        self.preview_fig = None
        self.retry_cnt = 0
        self.anchor_names = ["Bregma", "Lambda", "P3", "P4", "P5"]
    
    def Select_Anchor(self):
        while True:
            self.anchor = [] # anchor points, in seq. Y, X
            confirmed = False

            fig, ax = plt.subplots(figsize=(6, 6))
            ax.imshow(self.avr, cmap="gray")
            ax.set_title(
                "Left click: select | Backspace/u: undo | Enter: confirm | Esc: retry",
                fontsize=10
            )
            ax.axis("off")

            scatter = ax.scatter([], [], s=30)
            texts = []

            def redraw():
                nonlocal texts
                for t in texts:
                    t.remove()
                texts = []

                if len(self.anchor) == 0:
                    scatter.set_offsets(np.empty((0, 2)))
                    fig.canvas.draw_idle()
                    return
                
                pts = np.array(self.anchor) # (N, 2) [y, x]
                xy = np.stack([pts[:, 1], pts[:, 0]], axis=1)

                colors = []
                for i in range(len(pts)):
                    if i == 0:
                        colors.append("r") # Bregma
                    elif i == 1:
                        colors.append("g") # Lambda
                    else:
                        colors.append("b") # Others
                    
                    texts.append(
                        ax.text(
                            pts[i, 1] + 3,
                            pts[i, 0] + 3,
                            self.anchor_names[i],
                            color=colors[-1],
                            fontsize=9
                        )
                    )
                
                scatter.set_offsets(xy)
                scatter.set_color(colors)
                fig.canvas.draw_idle()

            def on_click(event):
                if event.inaxes != ax:
                    return
                if event.button == 1 and len(self.anchor) < 5:
                    y, x = int(event.ydata), int(event.xdata)
                    self.anchor.append([y, x])

                    if len(self.anchor) == 1:
                        logger.info("Selecting Bregma...")
                    elif len(self.anchor) == 2:
                        logger.info("Selecting Lambda...")
                    elif len(self.anchor) == 3:
                        logger.info("Selecting other anchor points (3 points in total)...")

                    redraw()

            def on_key(event):
                nonlocal confirmed

                if event.key in ["backspace", "u"] and len(self.anchor) > 0:
                    removed = self.anchor.pop()
                    logger.info(f"Undo last point: {removed}")
                    redraw()

                elif event.key == "enter":
                    if len(self.anchor) == 5:
                        confirmed = True
                        plt.close(fig)
                    else:
                        logger.info("Need exactly 5 anchor points to confirm")

                elif event.key == "escape":
                    confirmed = False
                    plt.close(fig)

            fig.canvas.mpl_connect("button_press_event", on_click)
            fig.canvas.mpl_connect("key_press_event", on_key)

            plt.show()

            if confirmed:
                logger.info("Anchors confirmed")
                plt.close()
                break

        self.anchor = np.array(self.anchor)
        self.realbreg = self.anchor[0]
        self.reallamb = self.anchor[1]
        self.point_demo = self.avr.copy()

        logger.info(
            f"Select points saved in self.anchor "
            f"(reselected {self.retry_cnt} times)"
        )
    
    def Fit_Align_Matrix(self):
        # calculate rescale index
        real_lbd = np.linalg.norm(self.realbreg - self.reallamb)
        self.scale = self.lbd / real_lbd

        # rotation and rescale
        self.slope, _ = np.polyfit(self.anchor[:, 1], self.anchor[:, 0], 1)
        if self.slope < 0: # negative angle
            self.rot_angle = np.degrees(np.arctan(self.slope)) + 90
        else:
            self.rot_angle = np.degrees(np.arctan(self.slope)) - 90
        if self.realbreg[0] > self.reallamb[0]:
            logger.warning("Graph is upside down! BE AWARE!")
            self.rot_angle += 180
        
        # transform
        padded_graph = np.pad(self.avr, ((self.pad_num, self.pad_num), (self.pad_num, self.pad_num)), mode="constant", constant_values=0)
        self.breg_pad = (int(self.realbreg[0] + self.pad_num), int(self.realbreg[1] + self.pad_num))
        # rotation matrix
        self.rot_mat = cv2.getRotationMatrix2D((self.breg_pad[1], self.breg_pad[0]), self.rot_angle, self.scale)
        
        result_raw = cv2.warpAffine(padded_graph, self.rot_mat, padded_graph.shape[1::-1], flags=cv2.INTER_NEAREST)
        self.LU_point = [self.breg_pad[0] - self.breg[0], self.breg_pad[1] - self.breg[1]]
        self.result = result_raw[self.LU_point[0]:self.LU_point[0] + self.height, self.LU_point[1]:self.LU_point[1] + self.width]
        show_img = cv2.circle(self.result, (self.breg[1], self.breg[0]), 3, (0, 0, 255), -1)
        
        if self.preview_fig is not None:
            plt.close(self.preview_fig)
        
        self.preview_fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(show_img, cmap="gray")
        ax.imshow(self.idmap_sym, alpha=0.15, cmap="jet")
        ax.set_title("Alignment preview")
        ax.axis("off")
        plt.show(block=False)

        while True:
            key = input("Press to continue (y/n): ")
            if key.lower() == 'y':
                break
            elif key.lower() == 'n':
                plt.close(self.preview_fig)
                self.preview_fig = None
                self.retry_cnt += 1
                logger.info(f"Reselecting anchors (retry {self.retry_cnt})")
                self.Select_Anchor()
                return self.Fit_Align_Matrix()
    
    def Transform_Series(self, stacks, intra=cv2.INTER_NEAREST):
        graph_num = len(stacks)
        transformed_stacks = np.zeros((graph_num, self.height, self.width), dtype="u2")

        # transform each frame
        for i in tqdm(range(graph_num), total=graph_num, desc="Transforming", leave=False):
            c_img = stacks[i, :, :].astype("f8")
            padded_img = np.pad(c_img, ((self.pad_num, self.pad_num), (self.pad_num, self.pad_num)), mode="constant", constant_values=0)
            rotted_img = cv2.warpAffine(padded_img, self.rot_mat, padded_img.shape[1::-1], flags=intra)
            cutted_img = rotted_img[self.LU_point[0]:self.LU_point[0] + self.height, self.LU_point[1]:self.LU_point[1] + self.width]
            transformed_stacks[i, :, :] = cutted_img
        
        return transformed_stacks


if __name__ == "__main__":
    avr = cv2.imread(r"", 0)
    MP = Match_Pattern(avr, bin=4)
    MP.Select_Anchor()
    MP.Fit_Align_Matrix()
