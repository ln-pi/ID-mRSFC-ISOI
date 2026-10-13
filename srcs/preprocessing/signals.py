import glob, os
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

def should_exclude(name, exclude_list):
    if exclude_list is None:
        return False

    for exclude in exclude_list:
        if isinstance(exclude, dict):
            for side, area in exclude.items():
                if area in name:
                    if side == "LR":
                        return True
                    elif side == "L" and name.endswith("_L"):
                        return True
                    elif side == "R" and name.endswith("_R"):
                        return True
        elif isinstance(exclude, str):
            if "_LR" in exclude:
                area = exclude.replace("_LR", "")
                if area in name:
                    return True
            elif "_L" in exclude:
                area = exclude.replace("_L", "")
                if area in name and name.endswith("_L"):
                    return True
            elif "_R" in exclude:
                area = exclude.replace("_R", "")
                if area in name and name.endswith("_R"):
                    return True
            else:
                if exclude in name:
                    return True

    return False

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

def extract_roi(matrix, mask_root, atlas_exclude=None, fmt="png", mode="interleaved"):

    if matrix.ndim != 3:
        logger.raise_(ValueError("Input matrix must be 3-D (Time, H, W)"))

    if mask_root is None:
        logger.raise_(IOError("FATAL: no mask_root provided!"))

    mask_files = [
        p for p in sorted(glob.glob(os.path.join(mask_root, f"*.{fmt}")))
        if os.path.splitext(os.path.basename(p))[0].casefold() != "chamber_mask"
    ]
    if not mask_files:
        raise FileNotFoundError(f"No ROI masks found in {mask_root} (excluding Chamber_Mask).")

    pairs = []

    for file in mask_files:
        name = os.path.basename(file).split(".")[0]

        if fmt == "png":
            from PIL import Image
            arr = np.array(Image.open(file))
        elif fmt == "npy":
            arr = np.load(file)
        else:
            logger.raise_(ValueError(f"Unsupported format: {fmt}"))

        if arr.shape != matrix.shape[1:]:
            raise ValueError(f"ROI mask {file} shape {arr.shape} does not match data shape {matrix.shape[1:]}.")
        mask_L = arr > 0
        mask_R = np.fliplr(mask_L)

        pairs.append((name, mask_L, mask_R))

    if mode == "interleaved":
        items = [(f"{n}_L", mL) for n, mL, _ in pairs] + [(f"{n}_R", mR) for n, _, mR in pairs]
        items = [x for pair in zip(items[:len(pairs)], items[len(pairs):]) for x in pair]
    elif mode == "symmetric":
        items = [(f"{n}_L", mL) for n, mL, _ in pairs] + [(f"{n}_R", mR) for n, _, mR in reversed(pairs)]
    else:
        logger.raise_(ValueError(f"Unknown mode: {mode}"))

    if atlas_exclude is not None:
        items = [(n, m) for n, m in items if not should_exclude(n, atlas_exclude)]

    roi_ts = []
    names = []

    for name, mask in items:
        time_series = matrix[:, mask]

        valid_voxel = ~np.isnan(time_series).all(axis=0)
        time_series = time_series[:, valid_voxel]

        if time_series.shape[1] == 0:
            roi_mean = np.full(matrix.shape[0], np.nan)
        else:
            roi_mean = np.nanmean(time_series, axis=1)

        roi_ts.append(roi_mean)
        names.append(name)

    roi_time_series = np.stack(roi_ts, axis=1)

    return roi_time_series, names
