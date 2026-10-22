import numpy as np

from scipy import signal

from srcs.utils.logger import setup_logger


logger = setup_logger()


def Filter_2D(matrix, mask, fs):
    sos = signal.butter(N=4, Wn=[0.009 / (fs / 2), (0.08 / (fs / 2))], btype="bandpass", output="sos")

    T, H, W = matrix.shape
    mask = mask.astype(bool)
    data_2D = matrix.reshape(T, -1)
    mask_1D = mask.reshape(-1)
    data_valid = data_2D[:, mask_1D]
    if np.isnan(data_valid).any() or np.isinf(data_valid).any():
        logger.raise_(ValueError("NaN or Inf found inside mask before filtering!"))
    data_filt = signal.sosfiltfilt(sos, data_valid, axis=0)
    filt_2D = np.full((T, H * W), np.nan, dtype=np.float64)
    filt_2D[:, mask_1D] = data_filt
    filt_data = filt_2D.reshape(T, H, W)
    return filt_data
