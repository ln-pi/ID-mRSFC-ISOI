import os

from tqdm import tqdm

from srcs.utils.logger import setup_logger
from srcs.utils.utils import parse_args
from srcs.io.readers import load_config
from srcs.prep.preprocess import Bin2Tiff, Affine_Atlas, merge_tiffs, drr_narrow
from srcs.prep.suite2p import run_isolated_process


args = parse_args()
config = load_config(args.config)
if config is None:
    raise ValueError("Config is empty (None). Exiting.")
logger = setup_logger(config["paths"]["log"])

raw_root = config["paths"]["raw"]
processed_root = config["paths"]["processed"]
output_root = config["paths"]["output"]
mask_root = config["paths"]["mask"]

dataset = config["variables"]["dataset"]
Session_ids = config["variables"]["Session_ids"]
Mouse_ids = config["variables"]["Mouse_ids"]
Trial_ids = config["variables"]["Trial_ids"]
Atlas = config["variables"]["Atlas"]
keepna = config["preprocess"]["keepna"]

processed_path = os.path.join(processed_root, dataset)
os.makedirs(processed_path, exist_ok=True)
output_path = os.path.join(output_root, dataset)
os.makedirs(output_path, exist_ok=True)

total = len(Session_ids) * len(Mouse_ids) * len(Trial_ids)
pbar = tqdm(total=total, desc="Processing", leave=False)
for session in Session_ids:
    for mouse in Mouse_ids:
        for trial in Trial_ids:

            raw_path = os.path.join(raw_root, dataset, session, mouse, trial)
            info_path = os.path.join(raw_path, "info.txt")

            processed_folder = os.path.join(processed_path, session, mouse, trial)
            os.makedirs(processed_folder, exist_ok=True)
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
                save_path=os.path.join(processed_path, "OIS_Corrected"),
            )

            # --- Step 3: Merge corrected TIFFs ---
            input_dir = os.path.join(processed_path, "OIS_Corrected", mouse, trial, "suite2p")
            merged_path = os.path.join(processed_folder, f"{mouse}_{trial}_merged.tif")

            merged_video = merge_tiffs(input_dir, merged_path)
            logger.info(f"Merged video saved at: {merged_path}, shape: {merged_video.shape}")

            # --- Step 4: Affine atlas registration ---
            Affine_Atlas(
                processed_path=processed_folder,
                session=session,
                mouse=mouse,
                trial=trial,
            )

            # --- Ster 5: dRR ---
            drr_narrow(
                processed_path=processed_folder,
                session=session,
                mouse=mouse,
                trial=trial
            )
