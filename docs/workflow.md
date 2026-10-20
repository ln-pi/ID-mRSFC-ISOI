# P1 执行流程与 P2 数据互通

## 入口及阶段

根目录 `run.py` 负责参数解析和配置顺序，每个配置使用独立 Python 子进程；`srcs/utils/pipeline.py` 负责阶段执行和状态记录。计算实现在 `srcs/workflow/`，模块导入不启动计算。原 scripts/01–04 只是兼容入口。

| 阶段 | 输入 | 主要输出 |
| --- | --- | --- |
| preprocess | info.txt、img*.bin、mask | processed 下 FiltDRR、预处理恢复记录 |
| fc | FiltDRR、mask | Correlation 周表、metadata/batch_reports |
| identification | FC 周表及审计 | results/combined、report/combined、metadata/cohort_selection.json |
| statistics | 完整识别矩阵及 statistics.pkl | report/tables、逐鼠自身匹配优势/百分位图 |

默认运行 Aseries 然后 Nseries，从 fc 到 statistics。显式 `--start preprocess` 才执行 Suite2p、交互配准及原 P1 的 GSR 预处理。阶段失败停止整个批次，阶段记录写为 failed/interrupted；缺失记录不等于执行异常，仍交给队列筛选。阶段记录是状态记录，不用于自动跳过计算；重新启动 identification 会重算识别。

```powershell
python run.py --configs configs/Aseries.yaml configs/Nseries.yaml --start identification
python run.py -c configs/Aseries.yaml --stage statistics
```

全部相对路径相对于项目根目录；配置哈希仍取原始 YAML 解析值，不因启动目录改变。日志追加到每份配置的 paths.log。`--dry-run` 只检查输入并显示阶段及输出位置，不运行科学计算、不写分析结果（日志仍可能创建）。

## FC 文件契约

文件名：`Correlation/<week>_correlation_matrices_<atlas>_excluded.pkl`。

内容为 pandas DataFrame，索引是 Mouse_ids，列是 Run_ids，每个有效单元格为带相同脑区行列标签的方阵。导入要求周表完整、鼠/Run 顺序与配置一致；矩阵必须有限、对称、对角线为 1、元素在 [-1,1]，且至少 3 个唯一 ROI。同一 atlas 的标签/顺序跨周跨鼠跨 Run 一致。None/NaN 表示缺失采集，损坏矩阵不作为普通缺失静默排除。

导入完整源输出根目录时，若有配置和 FC BatchReport，则检查范围一致；若只有 Correlation，重建每个单元格的状态。源报告中的失败状态继续保留，不因矩阵存在而提升为成功。源 PKL 原样复制，不修改其中 P2 的 config_tag；P1 的配置与来源绑定写入 metadata/fc_import_provenance.json，并在后续 identification 中重新核对文件 SHA256。

```powershell
python run.py -c configs/Aseries.yaml --start identification --import-fc "F:\path\to\source\Correlation"
```

首次导入不需要 raw、processed、ROI mask，也不需要重新运行 fc。后续直接使用 `--start identification`。导入目标中不同内容的同名文件会被拒绝；可通过新的配置输出路径隔离。导入后的输出禁止从 fc 重算，以免把另一套预处理混入同一结果。P1/P2 的完整配置哈希可不同，FC 通过采集维度和矩阵契约匹配。

只复制 FC 周表；FC_Matrix 下逐采集缓存和旧 FC 图不是后续计算所需数据。预处理等价性不在格式校验中推断：P1 自己预处理仍使用 Suite2p/GSR，P2 来源的 FC 使用其自身处理结果。来源配置可用时完整保存在 metadata 中。

## 输出分工

所有运行输出仍位于 `paths.output/<analysis_tag>/<dataset>`。顶层与 P2 一致：Correlation、results、report、metadata。P1 数值结果位于 results/combined；P2 的 results/longitudinal 不被当作 P1 识别结果读取。

- `results/combined/<condition>/numeric`：识别矩阵、逐鼠值、shuffle。
- `results/combined/<condition>/edges`：Phi、DP 数值。
- `results/combined/statistics.pkl`：04 的统计输入。
- `report/combined`：03 的识别和 edge 图。
- `report/tables`：ANOVA、posthoc、自身匹配优势/百分位等表格。
- `metadata`：配置、cohort、导入来源、BatchReport 和阶段状态。

旧输出不会自动迁移或删除。可以从旧 Correlation 显式导入后重新执行后续分析。旧 scripts/05 已移除，DD/ND 仍仅作为可选 legacy 函数。

## 针对性验证

`python -m unittest discover -s validation -v` 覆盖字节不变、缺失队列筛选、重复导入、维度不匹配、矩阵损坏、目标篡改、审计矛盾、目录映射及批次失败停止。

2026-10-20 使用本机已有 P2 输出 A=5a5c0248fb、N=ef021f99d6，分别导入 42 份周表（14 周 × 3 atlas），从项目外启动 identification→statistics。集成验证使用隔离配置和较低绘图分辨率、5 次 shuffle，只验证执行链及数值衔接，不作为正式重采样结果。验证输出位于 data/output/workflow_validation。
