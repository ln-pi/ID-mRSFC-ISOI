"""Run the project's real Suite2p wrapper on synthetic frames (no lab data)."""
import json
import os
from pathlib import Path
import sys
import tempfile
from importlib.metadata import version

# This executable smoke test is intentionally excluded from unittest discovery.
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("MKL_NUM_THREADS", "2")

import numpy as np
from scipy.ndimage import gaussian_filter
import tifffile

from srcs.prep.suite2p import run_isolated_process
from srcs.prep.preprocess import merge_tiffs
from srcs.prep.resume import check_registration, tiff_shape


def main():
    rng = np.random.default_rng(20260924)
    template = gaussian_filter(rng.normal(size=(64, 64)), 1.2)
    template = 5000 + 1500 * template / template.std()
    shifts = rng.integers(-2, 3, size=(240, 2))
    frames = np.stack([
        np.clip(np.roll(template, tuple(shift), axis=(0, 1)) + rng.normal(0, 20, template.shape), 1, 30000)
        for shift in shifts
    ]).astype(np.uint16)
    with tempfile.TemporaryDirectory(prefix="ois_suite2p_smoke_") as tmp:
        folder = Path(tmp) / "Run01" / "A01" / "00w"
        folder.mkdir(parents=True)
        original = folder / "A01_00w.tif"
        tifffile.imwrite(original, frames, photometric="minisblack")
        # Must be excluded by file_list even though it shares the input directory.
        tifffile.imwrite(folder / "A01_00w_merged.tif", frames[:13], photometric="minisblack")
        run_isolated_process(str(original), "A01", "00w", str(folder))
        reg = folder / "suite2p/plane0/reg_tif"
        assert check_registration(reg, original), "Registration output is incomplete"
        merged = folder / "A01_00w_merged.tif"
        merge_tiffs(str(folder / "suite2p"), str(merged))
        assert tiff_shape(merged) == frames.shape
        registered = tifffile.imread(merged)
        assert np.isfinite(registered).all() and registered.std() > 0
        outputs = np.load(folder / "suite2p/plane0/reg_outputs.npy", allow_pickle=True).item()
        assert len(outputs["yoff"]) == len(frames)
        assert len(outputs["xoff"]) == len(frames)
        settings = np.load(folder / "suite2p/plane0/settings.npy", allow_pickle=True).item()
        db = np.load(folder / "suite2p/plane0/db.npy", allow_pickle=True).item()
        assert db["nframes"] == 240, "Merged TIFF was accidentally added to input"
        summary = {
            "python": sys.version, "suite2p": version("suite2p"),
            "torch": version("torch"), "input_shape": list(frames.shape),
            "registered_shape": list(registered.shape),
            "registration_ran": True, "tiff_input_isolated": True,
            "settings": {key: settings[key] for key in ("run", "registration", "io")},
            "keep_movie_raw": db["keep_movie_raw"],
        }
        report = ROOT / "reports/environment_validation/suite2p_smoke.json"
        report.parent.mkdir(parents=True, exist_ok=True)
        report.write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")
        print(json.dumps({k: v for k, v in summary.items() if k != "settings"}, indent=2))
        print("PASS: real CPU registration, TIFF export, merge and resume validation")


if __name__ == "__main__":
    main()
