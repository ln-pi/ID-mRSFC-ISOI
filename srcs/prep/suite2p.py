import os
import subprocess
import sys

from srcs.utils.logger import setup_logger


logger = setup_logger()


def run_isolated_process(tiff_file, mouse, trial, save_path):
    script_content = f'''
import os
import numpy as np
import tifffile
from suite2p import run_s2p, default_settings


with tifffile.TiffFile(r"{tiff_file}") as tif:
    sample_frame = tif.pages[0].asarray()
    nframes = len(tif.pages)

if not np.isfinite(sample_frame).all() or sample_frame.sum() == 0:
    raise ValueError("First frame is non-finite or zero; Suite2p registration would be unstable.")

settings = default_settings()
settings.update({{
    "Ly": sample_frame.shape[0],
    "Lx": sample_frame.shape[1],
    "fs": 30.0,
    "diameter": 10,
    "torch_device": "cpu",
}})
settings["run"].update({{
    "do_registration": 1,
    "do_regmetrics": False,
    "do_detection": False,
    "do_deconvolution": False,
    "multiplane_parallel": False,
}})
settings["io"].update({{
    "combined": False,
    "save_mat": False,
    "save_NWB": False,
    "delete_bin": False,
    "move_bin": False,
}})
settings["registration"].update({{
    "two_step_registration": True,
    "nonrigid": False,
    "block_size": [128, 128],
    "nimg_init": min(500, nframes),
    "1Preg": True,
    "maxregshift": 0.1,
    "maxregshiftNR": 5,
    "smooth_sigma": 1.15,
    "subpixel": 10,
    "phasecorr": True,
    "snr_thresh": 1.2,
    "reg_tif": True,
    "reg_tif_chan2": False,
    "norm_frames": False,
}})
settings["extraction"].update({{
    "batch_size": 500,
    "neuropil_extract": False,
    "inner_neuropil_radius": 0,
}})

db = {{
    "data_path": [os.path.dirname(r"{tiff_file}")],
    "save_path0": r"{os.path.join(save_path, mouse, trial)}",
    "save_folder": "suite2p",
    "input_format": "tif",
    "nplanes": 1,
    "nchannels": 1,
    "functional_chan": 1,
    "keep_movie_raw": False,
}}

run_s2p(db, settings)

plane0 = os.path.join(db["save_path0"], db["save_folder"], "plane0")
reg_tif_dir = os.path.join(plane0, "reg_tif")
if not os.path.isdir(reg_tif_dir) or not any(name.lower().endswith((".tif", ".tiff")) for name in os.listdir(reg_tif_dir)):
    db_path = os.path.join(plane0, "db.npy")
    if not os.path.exists(db_path):
        raise FileNotFoundError(f"Suite2p did not create {{db_path}}")
    plane_db = np.load(db_path, allow_pickle=True).item()
    reg_file = plane_db.get("reg_file", os.path.join(plane0, "data.bin"))
    Ly = int(plane_db.get("Ly", sample_frame.shape[0]))
    Lx = int(plane_db.get("Lx", sample_frame.shape[1]))
    n_binary_frames = os.path.getsize(reg_file) // (np.dtype("int16").itemsize * Ly * Lx)
    if n_binary_frames == 0:
        raise ValueError(f"Suite2p binary is empty: {{reg_file}}")

    os.makedirs(reg_tif_dir, exist_ok=True)
    movie = np.memmap(reg_file, dtype="int16", mode="r", shape=(n_binary_frames, Ly, Lx))
    chunk = int(settings["extraction"]["batch_size"])
    for start in range(0, n_binary_frames, chunk):
        stop = min(start + chunk, n_binary_frames)
        out_file = os.path.join(reg_tif_dir, f"file{{start // chunk:03d}}.tif")
        tifffile.imwrite(out_file, np.asarray(movie[start:stop]), photometric="minisblack")
    del movie

print(f"finished: {mouse}_{trial}")
'''

    temp_script = f"temp_suite2p_{mouse}_{trial}.py"
    with open(temp_script, "w", encoding="utf-8") as f:
        f.write(script_content)

    try:
        logger.info(f"Running in isolated process: {mouse}_{trial}")
        result = subprocess.run(
            [sys.executable, temp_script],
            capture_output=True,
            text=True,
            timeout=24 * 3600,
        )

        if result.returncode == 0:
            logger.info(f"Process completed successfully: {mouse}_{trial}")
        else:
            logger.error(f"Process failed with return code {result.returncode}: {mouse}_{trial}")
            logger.error(f"Error output: {result.stderr[-4000:]}")
            raise RuntimeError(f"Suite2p failed for {mouse}_{trial}; see log for traceback.")
    except subprocess.TimeoutExpired:
        logger.exc(f"Process timed out: {mouse}_{trial}")
        raise

    finally:
        if os.path.exists(temp_script):
            os.remove(temp_script)
