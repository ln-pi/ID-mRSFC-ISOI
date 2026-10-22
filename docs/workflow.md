# 统一流程

工作根目录为 `F:/Projects/ID-mRSFC-ISOI/main`。入口不依赖 VS Code debug 或手工设置 PYTHONPATH。

```powershell
conda activate ID-mRSFC-ISOI
Set-Location F:\Projects\ID-mRSFC-ISOI\main
python run.py --dry-run
.\run.ps1 --dry-run
```

`run.ps1` 默认调用已激活环境的 python；也可通过 `ISOI_PYTHON` 指定解释器绝对路径。两种入口传递相同参数和退出码。相对配置及数据路径都以 main 为基准；参考 P2 的协调器按配置、阶段启动独立 worker。任意阶段失败返回非零码并停止后续工作；执行计划、阶段退出码和整体状态写入 `reports/automatic/<UTC时间_标识>/batch_summary.json`，可用 `--report-root` 改位置。

## 阶段

| 阶段 | 输入 | 输出 |
| --- | --- | --- |
| preprocess | 每个 week/mouse/run 的 info.txt、img*.bin | Suite2p 校正、人工 atlas 配准、FiltDRR |
| fc | 分别计算的每次采集信号、ROI mask | 带 ROI 标签的周 FC 表、逐记录审计 |
| edges | FC 表及审计 | 完整 N×W×Run×E 边数组 |
| pairwise | edge store | 每分区、每 Run 的 N×N×W×W；阶段选择与逐鼠指标 |
| descriptive | edge store | 同日模式相关/RMS/偏移、03w 变化、同组 LOO 偏离、ΔFC/逐边方差 |
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

`report` 默认绘制所有已有数值结果；`--report-section` 可选 `pairwise/combined/statistics/all`。combined 和 statistics 阶段不再直接绘图，统计数值统一写入 `results/statistics/`，原 `report/tables/` 中初版遗留文件不作为新流程输入。独立重绘无需原始 FC 表和 edge cache，也不会改写 `results/`。FC 阶段仍可生成采集诊断图，生成逻辑位于 visualization，不在 FC 计算核心中。

## 显式导入已有 FC

```powershell
python run.py -c configs/Aseries.yaml --stage edges --import-fc "F:\path\to\saved_output" --dry-run
python run.py -c configs/Aseries.yaml --start edges --import-fc "F:\path\to\saved_output"
```

接受单次输出根目录或其 correlation/Correlation 子目录。只读取可信本地 PKL；导入保持原字节，先核验全部表、范围、ROI 顺序及已有目标，再写入 main 的独立输出。缺失单元格保留为缺失。格式损坏或源审计不一致直接失败，不会为了导入而改变鼠名单。没有源审计时从单元格重建，并注明来源资料缺失。

新路径使用小写 correlation。仅有历史大写目录时允许读取；在区分大小写的文件系统发现两个不同目录则报错。源配置、原始 FC SHA256 和审计被记录。导入路径本身不能证明 GSR/Suite2p 等处理等价，也不表示本次决定了正式主分析的预处理来源。

## 数据和分析口径

基础层不全局删鼠，不拆阶段，不平均 Run。不同分区分开保存；E 是固定 ROI 顺序下严格上三角的唯一边。F_r 为原始 FC，F_z 为边级有限截断后的 atanh。

pairwise 对每个 Run 的 F_z 边向量求 Pearson，r 与相似度 z 分开保存。S[i,j,u,v]=S[j,i,v,u]；指定周对面板通常不对称。同一份 FC 与自身的 r=1 保留以便检查，但 statistical_mask=False、S_z=NaN。缺失和常量模式不能产生零相似度。

阶段从完整 S 选择。阶段内只取 u<v；baseline 内无跨周结果；阶段间保留有方向的矩形周块。输出每周对的候选鼠、周间隔、逐鼠指标及其窗口平均。先算每个周对的排名再汇总，候选集合变化后重算。

本版同时输出 available_pair 和 fixed_window 两种队列。前者只要求本周对两端有效；后者要求所选窗口的所有周有效。没有引入尚未商定的最低有效周数。窗口平均是描述性结果，保留 n_valid_pairs/n_requested_pairs 和候选数范围；周对不作为独立动物。

combined 从 F_z 平均后重新计算相关，保留有效 P1 14 周、7+7 拆分与旧相似度 z 的 ±4 截断。Phi/DP 为兼容 P1 继续平均原始 FC 边；配对与 ID 共用，不重复 50 次相同计算。奇数长度变为不重叠的 floor/ceil 两组，记录组大小，不复用一周制造相同样本量。

## 时间含义和研究范围

00w 是注射第一周；03w 是研究指定的注射后参照，不能称健康基线。N 使用起病/早/中/晚时间窗，A 使用相同观察时间窗。时间窗标签不代表每鼠临床进度相同。

当前实现提供基础数组、描述分析和兼容 combined。复杂 ICC、方差随机效应、正式动物层重采样、临床表型预测及残差 FC 分析没有因框架统一而自动完成。预处理沿用迁移模块，注册依赖 Suite2p 1.1.0；本次只做预处理模块导入检查，没有执行原始采集的人工配准或重做 FC。

旧二进制读取器对缺帧的 unsigned 整数表示限制仍在；不能将 keepna 当成已验证的原始缺失保留方案。预处理实际固定采用 30Hz 输入、0.009–0.08Hz 带通、5Hz 输出；删除了未被实现读取的 bin/filter/robust 配置项，避免造成可配置的假象。

## 验证

```powershell
python -m unittest discover -s validation -v
python validation/check_existing.py -c configs/Aseries.yaml --source "F:\path\to\saved_P1_A_output"
```

细节见 [现行架构与重构验收](refactor_p2.md)、[初版验收记录](unification_validation.md) 和 [数据契约](config.md)。交接文档保持为原设计记录。
