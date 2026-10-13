# 统一运行流程

从包含 numpy、pandas、scipy、matplotlib 等项目依赖的环境运行根目录 run.py。脚本自动将相对配置与数据路径定位到项目根目录。

```powershell
cd F:\Projects\ID-mRSFC-ISOI_refined
python run.py --help
python run.py                         # Aseries + Nseries，默认从 fc 开始
python run.py --dry-run               # 输入预检查与执行计划
python run.py -c configs/Aseries.yaml --start reliability
python run.py --configs configs/Aseries.yaml configs/Nseries.yaml --stage report
python run.py -c configs/Nseries.yaml --stage fc
python run.py -c configs/Nseries.yaml --stage reliability
python run.py -c configs/Nseries.yaml --stage weekly
python run.py -c configs/Nseries.yaml --stage longitudinal
python run.py -c configs/Nseries.yaml --stage report
python run.py -c configs/Nseries.yaml --stage phenotype
python run.py -c configs/Nseries.yaml --start preprocess
python run.py -c configs/Nseries.yaml --import-fc "源输出目录" --start reliability
python -m unittest discover -s tests
```

--stage 只执行指定阶段；--start 从指定阶段接续到末尾，两者互斥。-c/--config 指定单配置；--configs 指定多个配置，两者互斥。未指定配置时使用 Aseries、Nseries。默认不执行含人工配准的 preprocess。

可编辑安装后的 isoi 与 run.py 参数相同；项目路径迁移或入口更新后需要重新安装才能更新环境中的命令。直接运行 run.py 不依赖可编辑安装或 Debugger 的 PYTHONPATH。

各计算阶段保存数值结果；**FC、识别、统计图表统一在 report 阶段生成**。report 读取当前配置的 FC，并校验已保存 reliability、longitudinal 的数值文件与来源清单。缺少数值结果时先执行 --start reliability；仅补图使用 --stage report。表型仍需真实日期和观测数据。

结果按完整配置哈希隔离。report 允许读取路径迁移前生成的清单，但配置、输入文件名与内容、数值输出内容必须一致。修改图表配置仍会改变配置标签，不能直接给旧目录改名。

详情见 [图表报告说明](visualization_report.md) 与 [自动流程](automatic_workflow.md)。


## 统一输出布局

每个配置输出根目录保持不变，内部结构如下：

```text
Correlation/                    原 FC 数值、缓存与图表，位置不变
results/reliability/            同周可靠性
results/weekly/                 每周代表
results/longitudinal/           跨周识别结果、队列与逐鼠表
results/phenotype/              表型输入状态与对齐结果
report/figures/                 identification、statistics、reliability 等图表
report/tables/                  轨迹、结构指标及对照偏离汇总表
report/summary.json             报告摘要
report/figure_index.json        图表完整性索引
metadata/manifests/             阶段完成清单
metadata/batch_reports/         逐条执行记录
metadata/analysis_config.json   本次配置与代码标识
metadata/record_audit.csv        FC 输入审计
```

旧目录不自动迁移或删除。清除 Correlation 以外的旧结果后，从 `python run.py --start fc` 开始重跑，重建 FC 执行记录后顺序运行后续阶段；已有有效 FC 缓存遵循原有校验规则复用。保留 FC 执行记录时，也可从 `--start reliability` 开始。

`run.py --help` 列出各阶段、前置输入和启动示例。`--stage` 只执行一个阶段，`--start` 从指定阶段执行到主流程末尾。

## 版本控制范围

仅跟踪项目使用说明：config.md、workflow.md、automatic_workflow.md、visualization_report.md。交接、迁移工作记录、校验结果等仅本地保留，不跟踪。tests/ 不跟踪；data/ 仅跟踪 mask/，包括新增掩膜。
