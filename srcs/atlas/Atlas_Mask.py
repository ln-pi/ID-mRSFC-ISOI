import numpy as np
import os

from srcs.preprocessing.readers import Load_Variable
from srcs.utils.logger import setup_logger


logger = setup_logger()


class Mask_Generator(object):

    name = "Mask of visable parts"

    def __init__(self, bin=4):
        base_path = os.path.dirname(__file__)
        breg = Load_Variable(os.path.join(base_path, "Bregma.pkl"))

        if bin == 1:
            self.breg = breg
            self.idmap = np.load(os.path.join(base_path, "Raw_Area_ID.npy"))
            self.masks = Load_Variable(os.path.join(base_path, "Brain_Area_Masks.pkl"))
        elif bin == 2:
            self.breg = (int(breg[0] / 2), int(breg[1] / 2))
            self.idmap = np.load(os.path.join(base_path, "Raw_Area_ID_bin2.npy"))
            self.masks = Load_Variable(os.path.join(base_path, "Brain_Area_Masks_bin2.pkl"))
        elif bin == 4:
            self.breg = (int(breg[0] / 4), int(breg[1] / 4))
            self.idmap = np.load(os.path.join(base_path, "Raw_Area_ID_bin4.npy"))
            self.masks = Load_Variable(os.path.join(base_path, "Brain_Area_Masks_bin4.pkl"))
        else:
            raise IOError(f"Bin {bin} not supported!")
        self.all_area = list(set(self.masks["Area"]))
