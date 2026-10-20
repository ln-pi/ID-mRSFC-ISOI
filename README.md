# ID-mRSFC-ISOI — P1 Combined

ISOI 小鼠皮层功能连接的跨周汇总个体识别分析。数据结构始终为 week × mouse × run；默认在每个 Run 内按相同随机种子划分周次，比较两半的连接模式。Run01/Run02 是同日重复采集，00w 是观测起点。

## 运行

在已安装 `environment.yaml` 依赖的环境中，通过 PowerShell 执行：

```powershell
Set-Location F:\Projects\ID-mRSFC-ISOI\P1_Combined
python run.py --dry-run
python run.py
python run.py --configs configs/Aseries.yaml configs/Nseries.yaml --start identification
python run.py -c configs/Aseries.yaml --stage statistics
```

默认按 Aseries、Nseries 顺序，从 `fc` 开始执行到 `statistics`。`--stage` 只运行一个阶段；`--start` 连续运行指定阶段及其后续阶段。阶段顺序：`preprocess → fc → identification → statistics`。原始数据预处理含人工配准，必须显式选择 `--start preprocess`。配置路径、配置内相对数据路径均以 P1 项目根目录为基准，因此也可从其他目录通过绝对路径启动 `run.py`。

先检查配置中的 raw/processed/mask 路径。A/N 配置的鼠、周、Run、Atlas 名单与当前 P2 配置对齐，P1 的识别参数仍独立保留。发生阶段错误时返回非零退出码，并停止后续阶段和配置；缺失采集记录写入审计，由 strict/cover 队列规则处理。

## 复用 P2 的 Correlation

指定 **P2 某一次运行的输出根目录**，不是整个 `data/output`：

```powershell
python run.py -c configs/Aseries.yaml --start identification --import-fc "F:\path\to\P2_A_output" --dry-run
python run.py -c configs/Aseries.yaml --start identification --import-fc "F:\path\to\P2_A_output"
python run.py -c configs/Nseries.yaml --start identification --import-fc "F:\path\to\P2_N_output"
```

`--import-fc` 也接受 `.../Correlation` 文件夹，包括已经手动复制到 P1 目标目录中的文件夹。首次必须显式导入以创建当前配置的审计；此后直接 `--start identification`。导入仅复制分析所需的周汇总 PKL，保留原始字节、标签和矩阵；无需 raw、FiltDRR 或 FC_Matrix 缓存。不会覆盖内容不同的目标文件。

有源配置/BatchReport 时同时校验采集范围和审计；只有 Correlation 时按缺失单元格重建审计，来源信息标为不可用。矩阵格式损坏或标签不一致直接报错。格式共用不代表 Suite2p/GSR 等预处理等价，源配置保存在导入来源记录中。只导入可信的本地 PKL。

## 输出

与 P2 使用相同顶层结构，仍保留原来的完整配置哈希命名规则：

```text
data/output/<analysis_tag>/<dataset>/
├── Correlation/                  # 两个项目共用的 FC 周汇总表
├── results/combined/             # P1 识别、shuffle、edge 数值及 statistics.pkl
├── report/                       # 图形与统计汇总表（tables/）
│   └── combined/                 # 识别与 edge 图
└── metadata/                     # 配置、队列选择、来源及阶段记录
    ├── batch_reports/
    └── manifests/
```

P2 的纵向结果仍叫 `results/longitudinal`；P1 使用 `results/combined` 表明分析含义不同。两个配置哈希无需相同，导入会将 FC 放入当前 P1 配置对应的输出路径。

数字脚本保留为兼容入口，实现在 `srcs/workflow/`，可以通过 `run(ctx)` 单独调用。DD/ND legacy 功能仍只在 `srcs/utils/legacy_structure.py` 中，不加入正式流程。

## 验证与文档

```powershell
python -m unittest discover -s validation -v
```

[流程与迁移说明](docs/workflow.md) · [配置说明](docs/config.md) · [方法记录](docs/optimization.md) · [引用](CITATION.cff) · [许可证](LICENSE)
