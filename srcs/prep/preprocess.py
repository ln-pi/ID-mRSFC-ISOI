import cv2
import glob
import numpy as np
import os
import re
import tifffile

from scipy import signal

from srcs.utils.logger import setup_logger
from srcs.utils.utils import normalize_data, load_chamber_mask
from srcs.io.readers import Info_Reader, Graph_Reader_All
from srcs.comp.Align_Tools import Match_Pattern
from srcs.comp.Filters import Filter_2D


logger = setup_logger()


def Bin2Tiff(folder, info_path, output_tiff="output.tiff", keepna=False):
    if isinstance(keepna, str):
        keepna = keepna.strip().lower()
        if keepna in {"false", "0", "no", "off"}:
            keepna = False
        elif keepna in {"true", "1", "yes", "on"}:
            keepna = True
        else:
            logger.raise_(ValueError(f"Invalid keepna value: {keepna!r}"))

    # --- Step 1: get all *.bin files ---
    bin_files = glob.glob(os.path.join(folder, "img*.bin"))
    bin_files = sorted(bin_files, key=lambda x: int(os.path.basename(x).split(".")[0][4:]))

    if len(bin_files) == 0:
        logger.raise_(FileNotFoundError(f"No img*.bin files found in the folder: {folder}"))
    
    # --- Step 2: get channel names from info.txt ---
    info = Info_Reader(info_path)
    channel_names = info["Channel_Names"]

    # --- Step 3: get all frames with reader ---
    _, graph_all = Graph_Reader_All(
        bin_files,
        channel_names=channel_names,
        nan=keepna,
        log=True,
    )

    # --- Step 4: merge channels into TIFF stack ---
    stacks = []
    for ch in channel_names:
        stacks.append(graph_all[ch])
    
    final_stack = np.stack(stacks, axis=1).squeeze()
    num_frames = final_stack.shape[0]
    final_stack = final_stack.reshape(num_frames, final_stack.shape[1], final_stack.shape[2])
    final_stack = _sanitize_tiff_stack(final_stack)

    # --- Step 5: write multi-page TIFF ---
    tifffile.imwrite(output_tiff, final_stack, photometric="minisblack")

    logger.info(f"Tiff saved: {output_tiff}")

    return final_stack


def _sanitize_tiff_stack(stack):
    stack = np.asarray(stack)
    if stack.ndim != 3:
        logger.raise_(ValueError(f"Expected a 3-D TIFF stack, got shape {stack.shape}"))

    if not np.isfinite(stack).all():
        logger.warning("Non-finite values found before Suite2p; replacing them with finite values.")
        stack = np.nan_to_num(stack, nan=0, posinf=0, neginf=0)

    frame_sums = stack.reshape(stack.shape[0], -1).sum(axis=1)
    bad_frames = np.where(frame_sums == 0)[0]
    if len(bad_frames) > 0:
        logger.warning(f"Found {len(bad_frames)} zero frames before Suite2p; filling from nearby frames.")
        good_frames = np.where(frame_sums != 0)[0]
        if len(good_frames) == 0:
            logger.raise_(ValueError("All frames are zero; Suite2p cannot build a registration template."))
        for idx in bad_frames:
            previous = good_frames[good_frames < idx]
            fill_idx = previous[-1] if len(previous) else good_frames[0]
            stack[idx] = stack[fill_idx]

    return stack


def merge_tiffs(input_dir, merged_path):
    
    tiff_dir = os.path.join(input_dir, "plane0", "reg_tif")
    if not os.path.exists(tiff_dir):
        logger.raise_(FileNotFoundError(f"No such directory: {tiff_dir}"))
    
    tiff_files = []
    files = glob.glob(os.path.join(tiff_dir, "*.[tT][iI][fF]*"))
    for f in files:
        name = os.path.basename(f)
        m = re.search(r"(\d+)", name)
        idx = int(m.group(1)) if m else -1
        tiff_files.append((idx, f))
    
    if len(tiff_files) == 0:
        logger.raise_(FileNotFoundError(f"No TIFF files found in: {tiff_dir}"))
    
    tiff_files.sort(key=lambda x: x[0])
    tiff_path = [p for _, p in tiff_files]

    all_frames = []
    for path in tiff_path:
        with tifffile.TiffFile(path) as tif:
            all_frames.append(tif.asarray())
    
    merged_video = all_frames[0] if len(all_frames) == 1 else np.concatenate(all_frames, axis=0)
    tifffile.imwrite(merged_path, merged_video)

    return merged_video


def Affine_Atlas(processed_path, mouse, trial):
    tiff_path = os.path.join(processed_path, f"{mouse}_{trial}_merged.tif")
    with tifffile.TiffFile(tiff_path) as tif:
        video_data = tif.asarray()
    tiff_data = video_data

    avr = tiff_data.mean(axis=0)
    MP = Match_Pattern(avr, bin=4, lbd=4.2)
    MP.Select_Anchor()
    MP.Fit_Align_Matrix()

    transformed_stacks = MP.Transform_Series(tiff_data)
    os.makedirs(processed_path, exist_ok=True)
    affine_path = os.path.join(processed_path, f"{mouse}_{trial}_Affine.npy")
    np.save(affine_path, transformed_stacks)
    logger.info(f"Affined stacks saved: {affine_path}")
    logger.info(f"Affined stacks shape: {transformed_stacks.shape}")
    cv2.imwrite(os.path.join(processed_path, f"{mouse}_{trial}_Average_After.png"), transformed_stacks.mean(0).astype("u2"))


def drr_narrow(processed_path, mouse, trial, week, mask_root=os.path.join("data", "mask")):
    affine_data = np.load(os.path.join(processed_path, f"{mouse}_{trial}_Affine.npy"))
    chambermask = load_chamber_mask(mask_root, expected_shape=affine_data.shape[1:])
    mask_data = affine_data * chambermask

    drr = np.nan_to_num((mask_data - mask_data.mean(0)) / mask_data.mean(0))

    detrended_data = signal.detrend(drr, axis=0, type="linear")

    # GSR
    GSR_data = np.transpose(detrended_data * chambermask, (1, 2, 0))
    global_signal = (GSR_data[chambermask > 0]).mean(0)
    reshaped_data = GSR_data.reshape(-1, GSR_data.shape[-1])
    X = global_signal[:, np.newaxis]
    Y = reshaped_data.T
    beta = np.linalg.lstsq(X, Y, rcond=None)[0]
    residuals = Y - X @ beta
    regress_data = residuals.T.reshape(GSR_data.shape)
    NormalizedSeries = normalize_data(regress_data, mask=chambermask)
    ResidualGSR = np.transpose(NormalizedSeries, (2, 0, 1))

    #filtering
    filt_data = Filter_2D(ResidualGSR, chambermask, fs=30)

    step = int(30 / 5)
    data_ds = signal.decimate(filt_data, step, axis=0, zero_phase=True)

    np.save(os.path.join(processed_path, f"{mouse}_{trial}_{week}_FiltDRR.npy"), data_ds)
