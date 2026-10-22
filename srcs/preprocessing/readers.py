import numpy as np
import pandas as pd
import struct

from tqdm import tqdm

from srcs.utils.logger import setup_logger


logger = setup_logger()


def Info_Reader(txt_path):
    # read in txt
    info_dict = {}
    with open(txt_path, "r") as info_file:
        info_data = info_file.read().split("\n")
    # get capture frequency
    capture_freq = list(filter(lambda x: "FrameRateHz" in x, info_data))[0]
    capture_freq = float(capture_freq.split(":")[-1])
    # get img channel numbers
    img_channels = list(filter(lambda x: 'Illumination' in x, info_data))
    img_channels = list(map(lambda x: x.split(':')[-1].strip(), img_channels))
    ExposureSpeckleMsec = list(filter(lambda x: 'ExposureSpeckleMsec' in x, info_data))[0]
    ExposureSpeckleMsec = float(ExposureSpeckleMsec.split(':')[-1])

    # save element into dict
    info_dict["Imaging_Freq"] = capture_freq
    info_dict["Channel_Names"] = img_channels
    info_dict["ExposureSpeckleMsec"] = ExposureSpeckleMsec

    return info_dict


def Graph_Reader_All(img_file_lists, channel_names = ["Red"], log=True, nan=False):
    """
    Read all graph bin together, returning img stacks of each channel.

    Params:
    - img_file_lists(list): List of all img bin file.
    - channel_names(list): List of all channels' names in capturing.
    - log(bool): Whether to log.
    - nan(bool): Whether keep nan, or fill nan with nearby values/frames.

    Returns:
    - header_info(np.ndarray): Headfile info, in format [version,N_channels,Capture_Rate,Max_Capture_Rate,Max_Capture_Seconds]
    - all_graphs_channel(dict): Dictionary of all img files.
    """

    channel_num = len(channel_names)
    for i, c_img_name in tqdm(enumerate(img_file_lists), total=len(img_file_lists), desc="Extracting", leave=False):
        head_info, c_frame = Graph_Bin_Reader_Core(c_img_name)
        if i == 0:
            all_graphs = c_frame
        else:
            all_graphs = np.concatenate((all_graphs, c_frame), axis=0)

    if log is True:
        logger.info(f"Graph Resolution: {head_info[1]} x {head_info[2]}.")
        logger.info(f"Graph Number: {all_graphs.shape[0]}, in {channel_num} channels.")

    # seperate channel data
    all_channel_data = [[] for _ in range(channel_num)]
    for i in range(all_graphs.shape[0]):
        all_channel_data[i % channel_num].append(all_graphs[i, :, :])

    # save different channel into a dict
    del all_graphs
    all_graphs_channel = {}
    for i, c_channel in tqdm(enumerate(channel_names), total=len(channel_names), desc="Packing", leave=False):
        c_channel_data = np.array(all_channel_data[i])

        if nan is False:
            zero_ids = np.where(c_channel_data.sum(axis=(1, 2)) == 0)[0]
            if 0 in zero_ids:
                logger.info("First frame is zero, fill global average.")
                c_channel_data[0] = c_channel_data.mean(0)
                zero_ids = zero_ids[zero_ids != 0]
            for i, c_id in enumerate(zero_ids):
                c_channel_data[c_id] = c_channel_data[c_id - 1]

        all_graphs_channel[c_channel] = c_channel_data

    return head_info, all_graphs_channel

def Graph_Bin_Reader_Core(filename):
    """
    Read single graph bin file in OIS200 system.

    Params:
    - filename(str): Name of single bin file.

    Returns:
    - header_info(np.ndarray): Headfile info, in format [version,X dim,Y dim,img size in byte,max frame num per bin.]
    - graph_matrix(np.ndarray): Data of all graph array.
    """

    with open(filename, "rb") as file:
        header_bytes = file.read(20)
        data_bytes = file.read()

    # read header
    header = struct.unpack("5i", header_bytes)
    header_info = np.array(header)

    # unpack data, which have 3 int64 headers per frame
    total_framenum = len(data_bytes) // header_info[3]
    x_width = header_info[1]
    y_width = header_info[2]

    all_frames = np.zeros((total_framenum, x_width, y_width), dtype="u2")
    all_heads = np.zeros((total_framenum, 3), dtype="i4")

    for i in range(total_framenum):
        c_graph_bytes = data_bytes[i * header_info[3]:(i + 1) * header_info[3]]
        c_graph_head = np.array(struct.unpack(f"3q", c_graph_bytes[:24]))
        img_data = np.array(struct.unpack(f"{x_width * y_width}H", c_graph_bytes[24:]))
        all_frames[i, :, :] = img_data.reshape((x_width, y_width))
        all_heads[i, :] = c_graph_head

    file.close() # NOTE: with open will auto close, this is just for safety.
    del data_bytes

    num_adj = int(len(all_frames) + all_heads[:, 1].sum())
    all_frames_adj = np.zeros((num_adj, x_width, y_width), dtype="u2")
    miss_num = all_heads[:, 1]
    new_index = 0
    for i in range(len(miss_num)):
        # insert NaN images based on miss_num
        for _ in range(miss_num[i]):
            all_frames_adj[new_index] = np.full((x_width, y_width), np.nan)
            new_index += 1
        # insert the original frame
        all_frames_adj[new_index] = all_frames[i]
        new_index += 1

    del all_frames, all_heads
    return header_info, all_frames_adj


def Load_Variable(save_path):
    with open(save_path, "rb") as stream:
        return pd.read_pickle(stream)
