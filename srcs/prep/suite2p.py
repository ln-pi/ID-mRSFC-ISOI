import os
import subprocess

from srcs.utils.logger import setup_logger
from srcs.utils.utils import parse_args
from srcs.io.readers import load_config


args = parse_args()
config = load_config(args.config)
if config is None:
    raise ValueError("Config is empty (None). Exiting.")
logger = setup_logger(config["paths"]["log"])


def run_isolated_process(tiff_file, mouse, trial, save_path):

    script_content = f'''
import tifffile
import os
from suite2p import run_s2p, default_ops

# 读取图像信息
with tifffile.TiffFile(r"{tiff_file}") as tif:
    sample_frame = tif.pages[0].asarray()

# 配置Suite2p参数
ops = default_ops.default_ops()
ops.update({{
    # === 图像参数 ===
    "Ly": sample_frame.shape[0],      # 图像高度
    "Lx": sample_frame.shape[1],      # 图像宽度
    "fs": 30.0,                       # 帧率（Hz），根据实验调整
    
    # === 运动矫正核心 ===
    "do_registration": 1,
    "two_step_registration": True,
    "nonrigid": True,
    "block_size": [64, 64],           # OIS适合的块大小
    
    # === 模板构建 ===
    "nimg_init": 50,                  # 少量帧建立模板（避免血流信号污染）
    "1Preg": True,                     # 使用更稳健的模板更新
    
    # === 配准参数 ===
    "maxregshift": 0.02,              # 允许更大位移（红光对比度低）
    "smooth_sigma": 0.5,              # 更高平滑（提高信噪比）
    "subpixel": 20,                   # 高子像素精度
    "phasecorr": False,                # 相位相关法对红光更好
    
    # === 非刚性矫正 ===
    "snr_thresh": 1,                # 降低阈值（红光对比度低）
    
    # === 输出控制 ===
    "reg_tif": True,                  # 保存矫正TIFF
    "keep_movie_raw": False,          # 不保存原始二进制（节省空间）
    "save_mat": False,                # 不保存.mat文件
    "save_NWB": False,                # 不保存NWB
    
    # === 必须禁用的功能 ===
    "roidetect": 0,
    "diameter": 10,
    "spikedetect": 0,
    "denoise": 0,
    "smoothness_time": 0,
    "norm_frames": False,             # 不标准化帧（会破坏ΔR/R）
    "neuropil_extract": False,
    "inner_neuropil_radius": 0,

    # === 内存优化 ===
    "batch_size": 3000,                # 根据内存调整（OIS通常帧数少）
    "num_workers": 0,
    "use_gpu": False                # 单线程（避免内存溢出）
}})


# 运行
db = {{
    "data_path": [os.path.dirname(r"{tiff_file}")],
    "save_path0": r"{os.path.join(save_path, mouse, trial)}",
}}

run_s2p(db, ops)
print(f"finished: {mouse}_{trial}")
'''
    
    temp_script = f"temp_suite2p_{mouse}_{trial}.py"
    with open(temp_script, "w", encoding="utf-8") as f:
        f.write(script_content)
    
    try:
        logger.info(f"Running in isolated process: {mouse}_{trial}")
        result = subprocess.run(
            ["python", temp_script],
            capture_output=True,
            text=True,
            timeout=24 * 3600 # Timeout after ONE hour
        )

        if result.returncode == 0:
            logger.info(f"Process completed successfully: {mouse}_{trial}")
        else:
            logger.error(f"Process failed with return code {result.returncode}: {mouse}_{trial}")
            logger.error(f"Error output: {result.stderr[:500]}")
    except subprocess.TimeoutExpired:
        logger.exc(f"Process timed out: {mouse}_{trial}")
    
    finally:
        if os.path.exists(temp_script):
            os.remove(temp_script)
