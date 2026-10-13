# ID-mRSFC-ISOI

小鼠静息态连接的纵向 ISOI 分析流程：同周可靠性 → 每周代表连接 → 跨周个体识别 → 统一图表报告 → 表型接口。

统一入口为项目根目录 `run.py`，原 scripts/ 编号入口和 PowerShell 包装器已合并。使用包含项目依赖的 Python 环境运行：

```powershell
cd F:\Projects\ID-mRSFC-ISOI_refined
# 默认顺序运行 Aseries、Nseries，从现有 FiltDRR 开始
python run.py
# 两组已有数值结果，仅补绘全部报告图表
python run.py --stage report
# 单组、单阶段或接续运行
python run.py -c configs/Nseries.yaml --stage longitudinal
python run.py -c configs/Nseries.yaml --start reliability
```

入口自动切换到项目根目录，支持从其他工作目录用绝对路径启动；无需安装包即可运行。可选安装 `python -m pip install --no-deps --no-build-isolation -e .` 后，isoi 使用同一套参数；已经安装旧入口的环境需要重新安装才能更新命令。

所有分析图表由 **report** 阶段生成，包括 FC 热图、识别相似度矩阵、配对图、置换直方图、Phi/DP 环形图与矩阵、箱线图、均值误差曲线、全周对概览、DD/ND 及结构空间、可靠性和逐鼠轨迹。report/figure_index.json 记录图表和校验值；缺图或损坏时再次执行 report 可修复。单独补图要求当前配置的 reliability、longitudinal 数值结果及其阶段清单完整，不重算置换。

- [运行说明](docs/workflow.md)
- [批量自动流程](docs/automatic_workflow.md)
- [图表报告说明](docs/visualization_report.md)
- [配置说明](docs/config.md)

目录：run.py 统一入口，configs/ 配置，srcs/ 功能模块，tests/ 测试，docs/ 文档。原始预处理仍包含人工配准，需显式指定 --start preprocess。历史算法与静态快照保留在 srcs/legacy/。

输出目录统一为 `Correlation/`（保持原样）、`results/`（计算结果）、`report/`（图表与汇总表）、`metadata/`（配置、审计与清单）。删除 Correlation 以外的旧结果后，请用 `python run.py --start fc` 重建 FC 审计并运行后续阶段。

Git 仅纳入源代码、配置、项目使用说明和 `data/mask/`；本地测试、工作记录、校验报告及其他数据不纳入。

PowerShell 实时显示各阶段进度及完成状态。Aseries、Nseries 的日志分别按配置写入 `logs/Aseries.log`、`logs/Nseries.log`，每组所有阶段顺序追加；结束时显示日志路径并保留 warning 汇总。`reports/automatic/` 只保存批量执行摘要，不再保存阶段日志。
