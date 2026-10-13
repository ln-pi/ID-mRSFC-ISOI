"""Small, atomic result writers for interrupted analysis runs."""
import json
import os
import pickle
import tempfile
import time
import logging

import numpy as np


def _atomic(path, writer, binary=False):
    path = os.path.abspath(os.fspath(path))
    folder = os.path.dirname(path)
    os.makedirs(folder, exist_ok=True)
    handle = tempfile.NamedTemporaryFile(mode="wb" if binary else "w",
        dir=folder, delete=False, **({} if binary else {"encoding": "utf-8"}))
    try:
        with handle:
            writer(handle)
        for attempt in range(6):
            try:
                os.replace(handle.name, path)
                break
            except PermissionError:
                if attempt == 5:
                    raise
                time.sleep(0.05 * (2 ** attempt))
    finally:
        if os.path.exists(handle.name):
            try:
                os.unlink(handle.name)
            except OSError:
                logging.getLogger("project_logger").warning(
                    "Could not remove temporary output file: %s", handle.name)


def _json_value(value):
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    return value


def save_json(path, value):
    _atomic(path, lambda stream: json.dump(_json_value(value), stream,
        ensure_ascii=False, indent=2, allow_nan=False))


def save_pickle(path, value):
    _atomic(path, lambda stream: pickle.dump(value, stream, protocol=pickle.HIGHEST_PROTOCOL), binary=True)


def save_table(path, frame, index=True):
    _atomic(path, lambda stream: frame.to_csv(stream, index=index))

import pandas as pd
from scipy.io import savemat
from srcs.utils.logger import setup_logger
logger=setup_logger()

def mkdir(path, log=True):
    """
    Generate folder at input path if not exist, and do nothing if existed.

    Params:
    - path(str): Path to be generated.
    - log(bool): Whether to log.

    Returns:
    - (bool): Whether made a dir.
    """

    isExists = os.path.exists(path)
    if isExists:
        if log is True:
            logger.info(f"Folder {path} already exists.")
        return False
    else:
        os.makedirs(path)
        return True

def save_as_npy(expath, ai_signals, channel_names, graphs_all):
    """ Save data as numpy format. """
    np.save(os.path.join(expath, "ai_series.npy"), ai_signals)
    for i, c_name in tqdm(enumerate(channel_names), total=len(channel_names), desc="Saving", leave=False):
        c_graphs = graphs_all[c_name]
        avr_graph = c_graphs.mean(0).astype("u2")
        plt.imsave(os.path.join(expath, f"{c_name}.png"), avr_graph, cmap="gray", vmin=avr_graph.min(), vmax=avr_graph.mean() + 5 * avr_graph.std(), format="png")
        np.save(os.path.join(expath, f"{c_name}.npy"), c_graphs)

def save_as_mat(expath, ai_signals, channel_names, graphs_all):
    """ Save data as matlab format. """
    savemat(os.path.join(expath, "ai_series.mat"), {"ai_signals": ai_signals})
    for i, c_name in tqdm(enumerate(channel_names), total=len(channel_names), desc="Saving", leave=False):
        c_graphs = graphs_all[c_name]
        avr_graph = c_graphs.mean(0).astype("u2")
        plt.imsave(os.path.join(expath,f'{c_name}.png'), avr_graph, cmap='gray', vmin=avr_graph.min(), vmax=avr_graph.mean() + 5 * avr_graph.std(), format='png')
        savemat(os.path.join(expath,f'{c_name}.mat'),{f'{c_name}': c_graphs})

def Save_Variable(save_folder, name, variable, extend_name=".pkl"):
    """
    Save a variable as binary data.

    Params:
    - save_folder(str): Save path, which ending with folder name.
    - name(str): File name.
    - variable(Any): Data to save.
    - extend_name(str): Extended name of saved file. Default is ".pkl".
    """

    if not os.path.exists(save_folder):
        os.makedirs(save_folder)
    real_save_path = save_folder + r"\\" + name + extend_name
    fw = open(real_save_path, "wb")
    pickle.dump(variable, fw)
    fw.close()
    return True

def save_tensor(tensor, runs, subs, path):
    np.save(f"{path}\\corr_tensor.npy", tensor)
    np.save(f"{path}\\runs.npy", np.array(runs))
    np.save(f"{path}\\subs.npy", np.array(subs))

def Load_Variable(save_path, file_name=False):
    if file_name is False:
        real_file_path = save_path
    else:
        real_file_path = save_path + r"\\" + file_name
    if os.path.exists(real_file_path):
        pickle_off = open(real_file_path, "rb")
        loaded_file = pd.read_pickle(pickle_off)
        pickle_off.close()
    else:
        loaded_file = False

    return loaded_file
