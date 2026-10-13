import copy
import numpy as np
import pandas as pd

from scipy.stats import pearsonr
from tqdm import tqdm

from srcs.atlas.Atlas_Mask import Mask_Generator
from srcs.utils.logger import setup_logger


logger = setup_logger()


class Atlas_Data_Tools(object):

    name = "Atlas Data Easy Tool"

    def __init__(self, series, bin=4, min_pix=20):
        self.series = series
        self.std = series.std()
        self.MG = Mask_Generator(bin=4)
        self.all_area_name = self.MG.all_area
        self.Area_Function = {}
        self.min_pix = min_pix
        for i, c_area in enumerate(self.all_area_name):
            if "VI" in c_area:
                self.Area_Function[c_area] = "Visual"
            elif "SS" in c_area:
                self.Area_Function[c_area] = "Somatosensory"
            elif "MO" in c_area:
                self.Area_Function[c_area] = "Motion"
            elif "AUD" in c_area:
                self.Area_Function[c_area] = "Auditory"
            elif "RSP" in c_area:
                self.Area_Function[c_area] = "Retro"
            else:
                self.Area_Function[c_area] = "Other"
    
    def Get_All_Area_Response(self, keep_unilateral=False):
        """Average response by area"""
        self.Area_Response = pd.DataFrame(columns=["Area", "Function", "LR", "pixnum", "Series"])
        logger.info("Calculating All Brain Area Response...")
        for i, cloc in tqdm(enumerate(self.all_area_name), desc="GettingAllAreaResponse", leave=False):
            for j, hemi in enumerate(["L", "R"]):
                c_mask = self.MG.Get_Mask(cloc, hemi)
                combined_mask = c_mask & (~np.isnan(self.series).all(axis=0))
                # combined_mask = ((self.std * c_mask) != 0)
                if combined_mask.sum() > self.min_pix:
                    c_response = np.nanmean(self.series[:, combined_mask], axis=1)
                    c_func = self.Area_Function[cloc]
                    self.Area_Response.loc[len(self.Area_Response), :] = [cloc, c_func, hemi, combined_mask.sum(), c_response]

        on_area_name = list(set(self.Area_Response["Area"]))
        if keep_unilateral is False:
            for i, c_name in enumerate(on_area_name):
                if (self.Area_Response["Area"] == c_name).sum() < 2:
                    self.Area_Response = self.Area_Response[self.Area_Response["Area"] != c_name]

        part_L = self.Area_Response.groupby("LR").get_group("L").sort_values(by=['Function','Area'],ascending=False)
        part_R = self.Area_Response.groupby('LR').get_group('R').sort_values(by=['Function','Area'],ascending=True)
        self.Area_Response = pd.concat((part_L,part_R),axis=0).reset_index(drop=True)
    
    def Combine_Response_Heatmap(self):
        """Return heatmap of area reaponse matrix"""
        try:
            self.Area_Response
        except NameError:
            logger.warning("Area response not generated yet!")
            self.Get_All_Area_Response()
        
        response_heatmap = pd.DataFrame(columns=np.arange(len(self.Area_Response.iloc[0, -1])))
        logger.info("Combining heatmap of brain area response...")
        for i in tqdm(range(len(self.Area_Response)), total=len(self.Area_Response), desc="CombiningResponseHeatmap", leave=False):
            c_slice = self.Area_Response.iloc[i, :]
            c_name = c_slice["Area"] + "_" + c_slice["LR"]
            response_heatmap.loc[c_name, :] = c_slice["Series"]
        self.Area_Response_Heatmap = response_heatmap.astype("f8")
        return self.Area_Response_Heatmap
    
    def Get_Corr_Matrix(self, win_size=600, win_step=150):
        """Calculate corr matrix for area pair"""
        used_area_response = copy.deepcopy(self.Area_Response)
        used_area_response["Fullname"] = used_area_response["Area"] + "_" + used_area_response["LR"]
        used_area_fullname = list(used_area_response["Fullname"])

        frame_num = len(self.Area_Response.iloc[0, -1])
        win_num = 1 + (frame_num - win_size) // win_step
        self.Corr_Matrix = {}
        logger.info("Calculating correlation matrix of given data frames...")
        for i in tqdm(range(win_num), total=win_num, desc="GettingCorrMatrix", leave=False):
            c_start = i * win_step
            c_end = i * win_step + win_size
            c_corr_frame = pd.DataFrame(0.0, columns=used_area_fullname, index=(used_area_fullname))
            for j, name_A in enumerate(used_area_fullname):
                series_A = used_area_response[used_area_response["Fullname"] == name_A]["Series"].iloc[0][c_start:c_end]
                for k, name_B in enumerate(used_area_fullname):
                    series_B = used_area_response[used_area_response["Fullname"] == name_B]["Series"].iloc[0][c_start:c_end]
                    
                    valid_mask = (~np.isnan(series_A)) & (~np.isnan(series_B))

                    if valid_mask.sum() > 2:
                        c_r, _ = pearsonr(series_A[valid_mask], series_B[valid_mask])
                    else:
                        c_r = np.nan
                        
                    c_corr_frame.loc[name_A, name_B] = c_r
            self.Corr_Matrix[i] = c_corr_frame
        return self.Corr_Matrix


def Contra_Similar(series, bin=4):
    """
    Calculate contralateral similarity of given series.
    MASK SHOULD BE APPLIED BEFORE
    """
    _, height, width = series.shape
    mirrored_series = series[:, :, ::-1]

    MG = Mask_Generator(bin=bin)
    area_mask = MG.idmap > 0
    series_mask = series.std(0) > 0
    joint_mask = area_mask * series_mask
    joint_mask = joint_mask * joint_mask[:, ::-1]
    similar_map = np.zeros((height, width), dtype="f8")

    for i in range(height):
        for j in range(width):
            if joint_mask[i, j] == True:
                c_r, _ = pearsonr(series[:, i, j], mirrored_series[:, i, j])
                similar_map[i, j] = c_r
    
    return similar_map
