import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
from scipy.io import savemat
from tqdm import tqdm

from srcs.utils.logger import setup_logger


logger = setup_logger()


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


def one_key_save(inpath, expath):
    """Alternatively, if HD drop problem happens, use this function."""
    from srcs.io.readers import One_Key_Graph_Reader, One_Key_AI_Reader
    raw_r_dics = One_Key_Graph_Reader(inpath, nan=False)
    raw_ai_signals = One_Key_AI_Reader(inpath)

    # NOTE: returns are packed dicts
    logger.info(f"R value dics are in format {type(raw_r_dics)}, keys are: {raw_r_dics.keys()}")
    logger.info(f"AI value dics are in format {type(raw_ai_signals)}, in shape {raw_ai_signals.shape}")


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
