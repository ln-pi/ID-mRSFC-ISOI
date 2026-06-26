import os

from tqdm import tqdm

from srcs.utils.config import load_context
from srcs.prep.preprocess import Bin2Tiff, Affine_Atlas, merge_tiffs, drr_narrow
from srcs.prep.suite2p import run_isolated_process


ctx = load_context()
logger = ctx.logger
keepna = ctx.section("preprocess")["keepna"]

total = len(ctx.sessions) * len(ctx.mice) * len(ctx.trials)
pbar = tqdm(total=total, desc="Processing", leave=False)
try:
    for session in ctx.sessions:
        for mouse in ctx.mice:
            for trial in ctx.trials:
                pbar.set_postfix({"Session": session, "Mouse": mouse, "Trial": trial})
                logger.info(f"\n>>> Processing Session: {session}, Mouse: {mouse}, Trail: {trial}")

                raw_path = os.path.join(ctx.raw_root, ctx.dataset, session, mouse, trial)
                info_path = os.path.join(raw_path, "info.txt")

                processed_folder = ctx.processed_dir(session, mouse, trial)
                tiff_path = os.path.join(processed_folder, f"{mouse}_{trial}.tif")

                # --- Step 1: Bin to TIFF ---
                tiff = Bin2Tiff(
                    folder=raw_path,
                    info_path=info_path,
                    output_tiff=tiff_path,
                    keepna=keepna,
                )

                # --- Step 2: Suite2p motion correction (isolated process) ---
                run_isolated_process(
                    tiff_file=tiff_path,
                    mouse=mouse,
                    trial=trial,
                    save_path=os.path.join(ctx.processed_path, "OIS_Corrected"),
                )

                # --- Step 3: Merge corrected TIFFs ---
                input_dir = os.path.join(ctx.processed_path, "OIS_Corrected", mouse, trial, "suite2p")
                merged_path = os.path.join(processed_folder, f"{mouse}_{trial}_merged.tif")

                merged_video = merge_tiffs(input_dir, merged_path)
                logger.info(f"Merged video saved at: {merged_path}, shape: {merged_video.shape}")

                # --- Step 4: Affine atlas registration ---
                Affine_Atlas(
                    processed_path=processed_folder,
                    mouse=mouse,
                    trial=trial,
                )

                # --- Step 5: dRR ---
                drr_narrow(
                    processed_path=processed_folder,
                    mouse=mouse,
                    trial=trial
                )

                pbar.update(1)
finally:
    pbar.close()
