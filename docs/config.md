# 配置与数组契约

配置使用 YAML；相对路径以 main 根目录解析，绝对路径保持绝对路径语义。`configs/default.yaml` 是可用的 A 系列示例，不再是全 null 的占位模板。Aseries/Nseries 保留当前 10 鼠配置名单，实际分析鼠数由审计决定。

## 配置字段

| 字段 | 含义 |
| --- | --- |
| paths.raw / processed | 每次采集的原始/处理后数据根目录 |
| paths.mask | ROI PNG 与 Chamber_Mask.png；默认 ../data/mask |
| paths.output / log | 独立结果及日志位置 |
| variables.dataset | 路径附加目录 ID，可为空 |
| variables.Mouse_ids / Week_ids / Run_ids | 显式、有序、唯一的采集维度标签 |
| variables.Atlas | ROI 排除集合列表；每个分区独立 edge store |
| preprocess.keepna / resume | 旧读取行为参数；逐步骤校验续跑 |
| correlation.resume | 校验 FC/图缓存后复用 |
| analysis.split_axis | combined 的 week 或 run 拆分轴 |
| analysis.selection_mode | combined 对另一轴采用 single 或 cumulative |
| analysis.statistics_file | combined 聚合结果文件名，默认 statistics.pkl |
| ID.completeness_mode | strict 全配置完整队列，或 cover 最大完整鼠×连续周窗口；只影响 combined |
| ID.group_mode / session_seed | sequential/random、不重叠拆分及种子 |
| ID.shuffle / n_iter | combined 身份重标记；包含并列的上尾 (b+1)/(B+1) |
| Edge.seed | 保留的完整配置身份字段；固定 ID 配对使用 session_seed，不另做随机 edge 拆分 |
| Edge.threshold / gap / cross_only / save_mat | edge 图参数 |
| pairwise.phases | 阶段名称→周标签，仅用于完整数组之后的选择 |
| pairwise.phase_comparisons | 两个不重叠阶段的有向比较 |
| pairwise.week_positions | 可选周数值位置，必须覆盖全部 Week_ids 且严格递增；默认解析 00w 等数字标签 |
| pairwise.window_interpretation | N 疾病时间窗或 A 对齐观察窗的解释记录 |
| descriptive.reference_week | 自身变化参照；默认 03w |
| visualization | FC、combined 与统计图参数；pairwise 报告统一 r 色阶 [-1,1] |

完整 parsed config 仍参与输出哈希，绘图参数也会改变结果目录。没有把源预处理身份、分析身份、绘图身份静默拆为新缓存体系。根入口接受 `--configs` 顺序批处理，默认 Aseries 和 Nseries；`--import-fc` 要求单一 `-c` 并从 edges 开始。`--report-section` 是运行选择参数，不修改配置身份，可单独重绘 pairwise/combined/statistics；`--report-root` 指定批次状态记录目录。

已移除 Session_ids/Trial_ids、数值 group_mode 12/13/14、非空 repeat_id、supplement/DD/ND，以及重复固定配对的 Edge.repeat。旧数字脚本入口改为根 run.py/run.ps1；不对旧 YAML 进行隐式转换。

## FC / edge store

`correlation/<week>_correlation_matrices_<atlas>_excluded.pkl`：行是配置鼠，列是 Run，每个有效单元格为带一致 ROI 行列标签的 Pearson FC 方阵。原始表保持导入字节，以便来源复核和 FC 往返。

`results/edges/<atlas>.npz`，schema=`fc_edges_v1`：

| 键 | 形状/口径 |
| --- | --- |
| F_r | N×W×R×E，原始 Pearson FC 唯一边 |
| F_z | 同形，atanh(clip(F_r,-0.999999,0.999999)) |
| valid_record | N×W×R；表示 FC 有效，不等于模式相关可定义 |
| invalid_reason | N×W×R；缺失、审计失败、ROI 顺序/FC 无效的原因 |
| edge_index | 2×E，numpy.triu_indices(K,k=1) |
| metadata | JSON Unicode scalar；维度顺序、ID、ROI、schema、参数、源 FC/audit SHA256、源配置、Git commit/dirty 状态、代码逐文件哈希、包版本 |

全常量边可以是有效 FC，但无法定义模式 Pearson；其原因在派生 S 中标记。不根据阶段、Run 或完整 14 周预先筛掉动物。

每个 NPZ 原子写入，无 pickle 对象。JSON sidecar 包含数组文件 SHA256；后续步骤同时核对文件、源 FC/audit、配置和代码哈希。源或计算代码变化后显式重建 edges，不能静默复用旧数组。sidecar 的非原子双文件更新时间差最多导致明确校验失败，不会把不一致文件视为有效缓存。

## 每 Run 的 S

`results/pairwise/<atlas>/<run_id>.npz`，schema=`pairwise_similarity_v1`：

| 键 | 形状/口径 |
| --- | --- |
| S_r | N×N×W×W；Pearson_over_edges(F_z[i,u,run],F_z[j,v,run]) |
| S_z | 同形；atanh 仅用 nextafter(1,0) 处理数值边界，无 ±4 展示截断 |
| valid_similarity | 同形；模式 Pearson 可定义，包含有效自比较 |
| statistical_mask | 同形；另排除 i=j 且 u=v 的同一采集自比较 |
| invalid_reason | 同形 uint8；0=定义有效，1=端点无效，2=自比较排除 |
| pattern_reason | N×W；端点缺失/审计失败/常量边原因 |
| metadata | Run ID、分区、鼠/周顺序、尺度、算法、edge store 校验值和代码来源 |

无效位置 S_r/S_z=NaN；同一采集自比较 S_r=1、S_z=NaN。不同采集恰好相关 ±1 时使用 nextafter 边界得到有限 z，并明确记录此规则。周对面板行=i@u、列=j@v，反向面板为转置。

```python
from srcs.analysis.connectome import EdgeStore, load_arrays, SIM_SCHEMA
store = EdgeStore.load("/absolute/output/results/edges/AUDv.npz")
meta, values = load_arrays("/absolute/output/results/pairwise/AUDv/Run01.npz", SIM_SCHEMA)
panel = values["S_r"][:, :, 4, 7]       # 04w → 07w (按保存标签核对位置)
trajectory = values["S_r"][0, 0, :, :]  # 第一只鼠的跨周自身矩阵
```

常规阶段入口会额外校验来源和缓存，直接数组加载只做 schema/维度读取，适合明确知道来源的检查。

combined 原生接口为 `EdgeStore.average_halves(...)`，返回两个 mouse×edge 的 F_z 均值数组，再传给 `srcs.analysis.identification.identify(first, second)`。Phi/DP 通过 `srcs.analysis.edge_contribution.edge_contributions(...)` 直接消费 mouse×week×run×edge 的 F_r。旧的嵌套 DataFrame 半组接口及 Identify/Edgewise_Contribution 类已移除。

## 选择与描述结果

pairwise 的 subjects_by_pair.csv 保留每周对、鼠、Run、窗口、cohort、候选 ID、周间隔、尺度、相似度、自身优势、百分位及识别结果。subjects_by_window.csv 先按周对计算，再逐鼠取均值；同时保留有效/请求周对数和候选数范围。所有窗口/无结果原因及完整来源在 selection.json。

两个策略并行：available_pair 取本周对两端共同有效鼠；fixed_window 取整个所选周窗口均有效鼠。少于两只候选时不伪造排名，记入选择审计。baseline 内部零周对，保留 no_cross_week_pairs 状态。默认阶段内周对计数：onset/early/late 各 3，middle 6，baseline 0；early→middle 为 3×4 个有向周对。

本版没有预设更宽松的“阶段最低有效周数”、单 Run 替代或匹配间隔推断。CSV 保存实际标签对应的周数间隔，供后续相邻周/匹配间隔筛选。它不是实际采集日差；本次输入没有新增采集日期或失访原因资料。

`results/descriptive/topology_<atlas>.npz` 保存逐鼠逐周逐 Run 的相对 03w ΔFC、W×R×E 的鼠间样本方差（ddof=1）、有效鼠数和边标签。CSV 保留同日两个 Run 的模式 r/RMS/平均偏移、同组同周同 Run 的 leave-one-out RMS，以及相对自身参照变化。没有拟合 ICC、随机效应或组间推断。
