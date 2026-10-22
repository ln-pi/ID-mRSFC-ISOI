import os
from pathlib import Path

from tqdm import tqdm

from srcs.utils.dataset import resolve_recording, preprocessing_folder
from srcs.utils.batch import BatchReport
from srcs.preprocessing.signals import load_chamber_mask
from srcs.preprocessing.preprocess import Bin2Tiff, Affine_Atlas, merge_tiffs, drr_narrow
from srcs.preprocessing.suite2p import run_isolated_process, suite2p_runtime_signature
from srcs.preprocessing.resume import (
    StepRunner, tiff_shape, check_tiff, check_registration, check_alignment,
    check_filtered, registered_files, preserve_suite2p,
)


def process_recording(ctx, week, mouse, run):
    logger = ctx.logger
    settings = ctx.section("preprocess")
    raw_path = str(resolve_recording(ctx, week, mouse, run, raw=True))
    info_path = os.path.join(raw_path, "info.txt")
    bins = sorted(Path(raw_path).glob("img*.bin"))
    if not bins or not os.path.isfile(info_path):
        raise FileNotFoundError(f"Missing info.txt or img*.bin in {raw_path}")

    spatial_shape = load_chamber_mask(ctx.mask_root).shape
    chamber_path = next(p for p in Path(ctx.mask_root).iterdir()
                        if p.name.casefold() == "chamber_mask.png")
    folder = str(preprocessing_folder(ctx, week, mouse, run))
    tiff_path = os.path.join(folder, f"{mouse}_{run}.tif")
    suite_path = os.path.join(folder, "suite2p")
    reg_path = os.path.join(suite_path, "plane0", "reg_tif")
    merged_path = os.path.join(folder, f"{mouse}_{run}_merged.tif")
    affine_path = os.path.join(folder, f"{mouse}_{run}_Affine.npy")
    preview_path = os.path.join(folder, f"{mouse}_{run}_Average_After.png")
    filtered_path = os.path.join(folder, f"{mouse}_{run}_{week}_FiltDRR.npy")
    runner = StepRunner(folder, logger, resume=settings.get("resume", True))

    # --- Step 1: Bin to TIFF ---
    runner.run("step1_bin_to_tiff",
        lambda: Bin2Tiff(folder=raw_path, info_path=info_path, output_tiff=tiff_path,
                        keepna=settings.get("keepna", False)),
        [tiff_path], lambda: check_tiff(tiff_path),
        dependencies=[info_path, *bins], parameters={"keepna": settings.get("keepna", False)})

    # --- Step 2: Suite2p motion correction (isolated process) ---
    def correct():
        preserve_suite2p(suite_path)
        run_isolated_process(tiff_file=tiff_path, mouse=mouse, trial=run, save_path=folder)
    runner.run("step2_suite2p", correct, lambda: registered_files(reg_path),
        lambda: check_registration(reg_path, tiff_path), dependencies=[tiff_path],
        parameters=suite2p_runtime_signature())

    # --- Step 3: Merge corrected TIFFs ---
    runner.run("step3_merge", lambda: merge_tiffs(suite_path, merged_path),
        [merged_path], lambda: check_tiff(merged_path, tiff_shape(tiff_path)),
        dependencies=registered_files(reg_path))

    # --- Step 4: Affine atlas registration ---
    atlas_root = Path(__file__).resolve().parents[1] / "atlas"
    runner.run("step4_affine",
        lambda: Affine_Atlas(processed_path=folder, mouse=mouse, trial=run),
        [affine_path, preview_path],
        lambda: check_alignment(affine_path, preview_path, (tiff_shape(merged_path)[0], *spatial_shape)),
        dependencies=[merged_path, *[atlas_root / name for name in
            ("Bregma.pkl", "Raw_Area_ID_bin4.npy", "Brain_Area_Masks_bin4.pkl")]])

    # --- Step 5: dRR ---
    runner.run("step5_drr",
        lambda: drr_narrow(processed_path=folder, mouse=mouse, trial=run, week=week, mask_root=ctx.mask_root),
        [filtered_path], lambda: check_filtered(filtered_path, affine_path, spatial_shape),
        dependencies=[affine_path, chamber_path],
        parameters={"algorithm": "drr-gsr-filter-v1", "input_fps": 30, "output_fps": 5})

    return dict(input_path=raw_path)


def run(ctx):
    report = BatchReport(ctx, "01_preprocess")
    total = len(ctx.weeks) * len(ctx.mice) * len(ctx.runs)
    with tqdm(total=total, desc="Processing", leave=False) as pbar:
        for week in ctx.weeks:
            for mouse in ctx.mice:
                for run in ctx.runs:
                    pbar.set_postfix({"Week": week, "Mouse": mouse, "Run": run})
                    ctx.logger.info(f">>> Processing Week: {week}, Mouse: {mouse}, Run: {run}")
                    report.run(lambda: process_recording(ctx, week, mouse, run),
                        week=week, mouse=mouse, run=run,
                        input_root=os.path.join(ctx.raw_root, ctx.dataset))
                    pbar.update(1)
    report.finish()
    if not report.persisted:
        raise OSError("Batch report could not be saved")
    if any(row["status"] == "failed" for row in report.rows):
        raise RuntimeError(f"Stage contains failed records; see {report.path}")
