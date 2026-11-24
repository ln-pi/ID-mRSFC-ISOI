# 统一流程

本机可直接使用已安装的环境启动，无需先激活 Conda：

```powershell
Set-Location F:\Projects\ID-mRSFC-ISOI\P3_HighDimension
$python = "F:\CondaEnvs\ID-mRSFC-ISOI\python.exe"
& $python run.py --help
& $python run.py -c configs/SeedFC_Aseries_noGSR.yaml --stage seedfc --dry-run
```

上面的 SeedFC 检查命令已在本机验证。源文件含 22 个 seed；当前 FourNetworks 分区通过 exclude 保留指定的 14 个 ROI，noGSR 两组存在有效记录，可运行后续分析。原生 FC 入口还需先在配置中设置实际可用的处理后信号和 mask 路径。

数据来源、跨系列比较和重绘使用同一入口。

```powershell
python run.py -c configs/SeedFC_Aseries_noGSR.yaml --stage seedfc --dry-run
python run.py --compare <A系列输出目录> <N系列输出目录>
python run.py --from <旧输出目录> --stage report
```

输入及插件配置见 [插件契约](plugin_contract.md)。重绘旧输出时，新报告保存在当前工作树默认输出根目录的 `_redraw/`。

工作根目录为 `F:/Projects/ID-mRSFC-ISOI/P3_HighDimension`。入口不依赖 VS Code debug 或手工设置 PYTHONPATH。

```powershell
conda activate ID-mRSFC-ISOI
Set-Location F:\Projects\ID-mRSFC-ISOI\P3_HighDimension
python run.py --dry-run
.\run.ps1 --dry-run
```

`run.ps1` 默认调用已激活环境的 python；也可通过 `ISOI_PYTHON` 指定解释器绝对路径。两种入口传递相同参数和退出码。相对配置及数据路径都以 P3 工作树为基准；参考 P2 的协调器按配置、阶段启动独立 worker。任意阶段失败返回非零码并停止后续工作；执行计划、阶段退出码和整体状态写入 首份配置输出根目录下的 `_batches/<UTC时间_标识>/batch_summary.json`，可用 `--report-root` 改位置。

## 阶段

| 阶段 | 输入 | 输出 |
| --- | --- | --- |
| preprocess | 每个 week/mouse/run 的 info.txt、img*.bin | Suite2p 校正、人工 atlas 配准、FiltDRR |
| fc | 分别计算的每次采集信号、ROI mask | 带 ROI 标签的周 FC 表、逐记录审计 |
| edges | FC 表及审计 | 完整 N×W×Run×E 边数组 |
| pairwise | edge store | 每分区、每 Run 的 N×N×W×W；阶段选择与逐鼠指标 |
| descriptive | edge store | 同日模式相关/RMS/偏移、03w 变化、同组 LOO 偏离、ΔFC/逐边方差 |
| stability | descriptive 与 edge store | 逐边方差分量与稳定性 |
| longitudinal | descriptive 与 edge store | 整鼠轨迹和网络随机斜率 |
| combined | edge store；strict/cover 选择规则 | P1 汇总 ID、shuffle、Phi/DP |
| statistics | 保存的 combined 数值 | 逐鼠自身匹配优势/百分位、重复测量统计与表 |
| report | 保存的 pairwise、combined、statistics 数值 | 阶段图、ID/edge 图、统计图及输入校验索引 |

默认 Aseries→Nseries，从 fc 连续执行。预处理含人工操作，必须显式指定。`--stage` 单阶段，`--start` 连续执行该阶段及后续阶段。只需要新基础数组时，可依次单独执行 edges、pairwise、descriptive、report，不必计算 combined。

```powershell
python run.py -c configs/Nseries.yaml --stage preprocess
python run.py -c configs/Nseries.yaml --stage fc
python run.py -c configs/Nseries.yaml --start edges
python run.py -c configs/Nseries.yaml --stage combined
python run.py -c configs/Nseries.yaml --stage report
python run.py -c configs/Nseries.yaml --stage report --report-section combined
python run.py -c configs/Nseries.yaml --stage report --report-section statistics
```

`report` 默认绘制所有已有数值结果；`--report-section` 可选 `pairwise/descriptive/stability/longitudinal/combined/statistics/all`，已启用的功能插件还可注册报告分节。combined 和 statistics 阶段不再直接绘图，统计数值统一写入 `results/statistics/`，原 `report/tables/` 中初版遗留文件不作为新流程输入。独立重绘从保存的结果读取。阶段相似性展示窗口改变阶段成员时，还会读取经校验的 edge store。FC 阶段仍可生成采集诊断图，生成逻辑位于 visualization，不在 FC 计算核心中。

幅度与拓扑分析：

```powershell
python run.py --stage edges
python run.py --stage descriptive
python run.py --stage report --report-section descriptive
```

默认对每个 Run 及所有 Run 的 Fisher-z FC 均值分别计算欧式/RMS/1−r 距离、LOO ΔFC、群体离散度、逐边样本方差及 ROI/网络方向和贡献。`results/descriptive/` 保存两套队列的全部数值、选择规则和描述性离群示例；`report/descriptive/` 保存 PNG/PDF 和索引。参数和公式见 [配置与数据契约](config.md)。

## SeedFC 批次补充报告

```powershell
conda activate ID-mRSFC-ISOI
python run.py --configs configs/SeedFC_Aseries_noGSR.yaml configs/SeedFC_Nseries_noGSR.yaml --start seedfc
```

所有输入均启用 SeedFC 和 `plugins.function.seedfc_supplement.enabled` 时，完整报告批次在主流程和适用的 compare 完成后，按 atlas 生成补充报告。单系列完整批次也适用。 `--dry-run` 在 `batch_summary.json` 记录补充步骤和输入目录。单阶段和局部报告只执行指定工作。补充失败时，汇总记录失败原因和已生成的产物，可用如下入口从原先输出重试。

```powershell
python -m plugins.function.seedfc_supplement --outputs "<A输出目录>" "<N输出目录>" --atlas FourNetworks --config configs/SeedFC_Aseries_noGSR.yaml
python -m plugins.function.seedfc_supplement --outputs "<A输出目录>" "<N输出目录>" --atlas FourNetworks --config configs/SeedFC_Aseries_noGSR.yaml --redraw
python run.py --from "<已有输出目录>" --stage report --report-section longitudinal --visualization-config "<绘图配置.yaml>"
python run.py --compare "<A输出目录>" "<N输出目录>" --visualization-config "<绘图配置.yaml>"
```

补充分析的默认位置为输出根目录的 `_supplement/numeric/` 和 `_supplement/reports/`。手动补算可用 `--destination` 指定其他根目录。数值清单核验输入 edge store、数值选项及产品散列。修改颜色、DPI 或展示起点后，`--redraw` 要求命中这份已核验的数值缓存，并生成独立的报告。 `--from` 将主流程重绘图写入默认输出根目录的 `_redraw/`。主流程已有报告时，不同选项或输入的新图写入 `rendered_reports/`。跨系列 compare 按已保存输入及展示设置定位结果，直接读取两组共同具备的 pairwise、descriptive 和 longitudinal 产物。

## 显式导入已有 FC

```powershell
python run.py -c configs/Aseries.yaml --stage edges --import-fc "F:\path\to\saved_output" --dry-run
python run.py -c configs/Aseries.yaml --start edges --import-fc "F:\path\to\saved_output"
```

接受单次输出根目录或其 correlation/Correlation 子目录。只读取可信本地 PKL；导入保持原字节，先核验全部表、范围、ROI 顺序及已有目标，再写入 P3 工作树的独立输出。缺失单元格保留为缺失。格式损坏或源审计不一致直接失败，不会为了导入而改变鼠名单。没有源审计时从单元格重建，并注明来源资料缺失。

新路径使用小写 correlation。仅有历史大写目录时允许读取；在区分大小写的文件系统发现两个不同目录则报错。源配置、原始 FC SHA256 和审计被记录。导入路径本身不能证明 GSR/Suite2p 等处理等价，也不表示本次决定了正式主分析的预处理来源。

## 数据和分析口径

基础层不全局删鼠，不拆阶段，不平均 Run。不同分区分开保存；E 是固定 ROI 顺序下严格上三角的唯一边。F_r 为原始 FC，F_z 为边级有限截断后的 atanh。

pairwise 对每个 Run 的 F_z 边向量求 Pearson，r 与相似度 z 分开保存。S[i,j,u,v]=S[j,i,v,u]；指定周对面板通常不对称。同一份 FC 与自身的 r=1 保留以便检查，但 statistical_mask=False、S_z=NaN。缺失和常量模式不能产生零相似度。

阶段从完整 S 选择。阶段内只取 u<v；baseline 内无跨周结果；阶段间保留有方向的矩形周块。输出每周对的候选鼠、周间隔、逐鼠指标及其窗口平均。先算每个周对的排名再汇总，候选集合变化后重算。

本版同时输出 available_pair 和 fixed_window 两种队列。前者只要求本周对两端有效；后者要求所选窗口的所有周有效。没有引入尚未商定的最低有效周数。窗口平均是描述性结果，保留 n_valid_pairs/n_requested_pairs 和候选数范围；周对不作为独立动物。

combined 从 F_z 平均后重新计算相关，保留有效 P1 14 周、7+7 拆分与旧相似度 z 的 ±4 截断。Phi/DP 为兼容 P1 继续平均原始 FC 边；配对与 ID 共用，不重复 50 次相同计算。奇数长度变为不重叠的 floor/ceil 两组，记录组大小，不复用一周制造相同样本量。

## 时间含义和研究范围

00w 是注射第一周；03w 是研究指定的注射后参照，不能称健康基线。N 使用起病/早/中/晚时间窗，A 使用相同观察时间窗。时间窗标签不代表每鼠临床进度相同。

当前实现提供基础数组、幅度/拓扑、稳定性、纵向轨迹与网络随机斜率，以及兼容 combined。stability 对完整平衡队列估计逐边方差分量和条件稳定 ICC；longitudinal 提供动物 block bootstrap 轨迹区间、固定/可用队列对照及网络随机斜率诊断。残差模式的跨周识别、临床表型预测、非线性病程模型和信息性失访校正仍未实现。预处理沿用迁移模块，注册依赖 Suite2p 1.1.0；不因下游升级自动重做原始采集配准或 FC。

稳定性和纵向演化在 descriptive 后运行：

```powershell
python run.py --stage stability
python run.py --stage longitudinal
python run.py --stage report --report-section stability
python run.py --stage report --report-section longitudinal
```

首次升级先重建 edges/descriptive。也可使用 `--start descriptive` 连续运行后续全部阶段；这包含原 combined/statistics，不是仅运行新增两个模块。新阶段的方差模型、固定队列限制和输出详见 [数据契约](config.md)。

旧二进制读取器对缺帧的 unsigned 整数表示限制仍在；不能将 keepna 当成已验证的原始缺失保留方案。预处理实际固定采用 30Hz 输入、0.009–0.08Hz 带通、5Hz 输出；删除了未被实现读取的 bin/filter/robust 配置项，避免造成可配置的假象。

## 验证

```powershell
python -m unittest discover -s validation -v
```

幅度与拓扑检查包含独立公式恒等式、缺失与常量模式、ROI/网络贡献守恒及独立重绘。详细参数见 [数据契约](config.md)。


## 阶段矩阵与 A/N 合并展示

pairwise.summary_phases 指定 baseline、early、middle、late，周名单沿用 pairwise.phases。Run 独立保存，并额外计算所有有向 Run 对。results/pairwise/<atlas>/phase_similarity.npz 保存每个 Run 对、阶段对、队列和平均方法的 N×N 矩阵，以及每格有效周对数。correlate_mean_fc 先逐鼠在阶段内平均 Fisher-z FC 边，再做边模式 Pearson；mean_correlations 先计算各周对的模式 Pearson，再对 r 作算术平均。后者排除同鼠、同 Run、同周的同一采集自比较；baseline 同 Run 的自身对角线缺失，异鼠比较保留。前者同 Run 同阶段的自身对角线屏蔽。跨周自身相关仍保留，各矩阵格的有效周对数分别保存。

available_phase 要求各端点阶段在对应 Run 的全部周均有效；fixed_window 进一步要求四阶段窗口全部周在参与的 Run 中有效。缺失与常数模式保留为 NaN，计数随矩阵保存。

longitudinal 的 phase_subjects.csv 和 phase_dispersion.csv 汇总已算好的每周指标。阶段 LOO RMS 是每周 RMS 的均值，不先平均偏移向量。有符号偏离均值单独保存，允许方向抵消。每指标保留有效周数；阶段图是描述性均值，不平均逐周置信区间。

先运行各系列计算和 report，再使用 python run.py --compare <A输出目录> <N输出目录>。比较输出的 report/series_overlays 保存统一坐标中的每周与四阶段 trajectories、方向热图和跨 Run 阶段矩阵。A/N 使用固定颜色，Run 使用固定点形，模式相关色条在 [-1,1] 的理论范围内由配置决定，方向 RMS 两侧共用非负色条。各组 LOO 模板分别计算。
