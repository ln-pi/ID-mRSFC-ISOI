"""Signal normalization and the acquisition chamber mask."""
import glob
import os
import numpy as np
from srcs.utils.logger import setup_logger

logger = setup_logger()


def normalize_data(data, mask=None, axis=-1, eps=1e-8, mask_fill=np.nan):
    """
    Z-score normalize data along a given axis.

    For OIS stacks arranged as H x W x T, the default normalizes each pixel's
    time series. If mask is provided, pixels outside mask are filled with
    mask_fill, which defaults to NaN.
    """
    data = np.asarray(data, dtype=np.float64)

    if mask is None:
        mean = np.nanmean(data, axis=axis, keepdims=True)
        std = np.nanstd(data, axis=axis, keepdims=True)
        normalized = (data - mean) / (std + eps)
        return np.nan_to_num(normalized, nan=0.0, posinf=0.0, neginf=0.0)

    mask = np.asarray(mask, dtype=bool)
    if mask.shape != data.shape[:mask.ndim]:
        logger.raise_(ValueError(
            f"Mask shape {mask.shape} is incompatible with data shape {data.shape}"
        ))

    normalized = np.full(data.shape, mask_fill, dtype=np.float64)
    valid_data = data[mask]
    mean = np.nanmean(valid_data, axis=axis, keepdims=True)
    std = np.nanstd(valid_data, axis=axis, keepdims=True)
    valid_normalized = (valid_data - mean) / (std + eps)
    normalized[mask] = np.nan_to_num(valid_normalized, nan=0.0, posinf=0.0, neginf=0.0)

    return normalized

def load_chamber_mask(mask_root, expected_shape=None):
    """Load the shared chamber mask, allowing filename case differences."""
    from PIL import Image
    files = [
        p for p in glob.glob(os.path.join(mask_root, "*"))
        if os.path.basename(p).casefold() == "chamber_mask.png" and os.path.isfile(p)
    ]
    if len(files) != 1:
        raise FileNotFoundError(f"Expected one Chamber_Mask.png in {mask_root}; found {len(files)}.")
    with Image.open(files[0]) as image:
        mask = np.asarray(image.convert("L")) > 0
        if "A" in image.getbands():
            mask &= np.asarray(image.getchannel("A")) > 0
    if expected_shape is not None and mask.shape != tuple(expected_shape):
        raise ValueError(f"Chamber mask shape {mask.shape} does not match data shape {tuple(expected_shape)}.")
    if not mask.any():
        raise ValueError("Chamber mask contains no valid pixels.")
    return mask
