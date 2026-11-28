# P3 批次输出布局

每次独立执行产生一个新批次。默认父目录为配置中的 `paths.output`（常用 `data/output/`）；`--output-root` 可指定其他父目录。批次名为 `Batch_<UTC时间>_<数据组合>_<任务>`，同秒同名时顺延 `_02`、`_03`。系列名仍由来源、系列、预处理分支、周次范围、鼠数、Run 数、短数据身份散列和可选 `--tag` 组成。

```text
data/output/
└── Batch_20261125-054306Z_SeedFC-AN-noGSR_Core-Geometry/
    ├── batch_summary.json
    ├── logs/
    │   ├── batch.log
    │   ├── Aseries-noGSR.log
    │   └── Nseries-noGSR.log
    ├── SeedFC-Aseries-noGSR_00w-13w_10mice_2run_<sha>_Geometry/
    │   ├── correlation/
    │   ├── results/
    │   ├── report/
    │   └── metadata/
    ├── SeedFC-Nseries-noGSR_00w-13w_14mice_2run_<sha>_Geometry/
    │   └── 同上
    └── compare/
        ├── results/
        │   ├── core/
        │   ├── geometry/
        │   └── seedfc_supplement/FourNetworks/
        ├── report/
        │   ├── core/
        │   ├── geometry/
        │   └── seedfc_supplement/FourNetworks/
        └── metadata/
```

目录按实际生成的产品创建。显式 `--compare` 只比较指定的两个已有结果：不同系列同一处理分支为独立比较，同一系列不同处理分支（如 GSR/noGSR）为配对比较。两种结果都直接写入 `compare/`，不为比较对象另建一层目录。GSR/noGSR 需要显式指定，不由 `--configs` 自动触发。SeedFC 双系列补充放入 `compare/` 便于查阅，但不是组间推断；单系列补充放入该系列自己的 `results/seedfc_supplement/<atlas>/` 和 `report/seedfc_supplement/<atlas>/`。没有对比时不创建 `compare/`。清单保留完整校验散列；公开目录不使用长散列命名。

```powershell
# A/N 新批次；可先附加 --dry-run，查看计划（不落盘）。
python run.py --configs configs/SeedFC_Aseries_noGSR.yaml configs/SeedFC_Nseries_noGSR.yaml --start seedfc --function geometry --tag Geometry

# 在明确指定的同一批次运行后续阶段；配置、功能开关和 tag 须与原系列一致。
python run.py -c configs/SeedFC_Aseries_noGSR.yaml --stage geometry --function geometry --tag Geometry --resume-batch "data/output/<批次名>"

# 两条命令均创建新批次，仅从原系列读入，不修改历史批次。
python run.py --compare "data/output/<批次名>/<A系列目录>" "data/output/<批次名>/<N系列目录>" --function geometry
python run.py --from "data/output/<批次名>/<系列目录>" --stage report --report-section geometry --function geometry
```

`--dry-run` 不创建批次摘要或日志；流程只校验所选计划的首阶段。一次正式执行的阶段、退出码和失败原因见批次根的 `batch_summary.json`；续跑前的摘要归档到该批次 `metadata/previous_runs/`。旧配置中的 `paths.log` 不再控制日志位置。若旧脚本使用 `--report-root`，请改用 `--output-root`；兼容参数现在指向整个新批次的父目录，不只移动摘要。

新批次的内部产物清单使用相对于各清单目录的路径；采集目录及独立重绘、比较所引用的外部系列仍记录来源。历史系列按原路径可作为 `--from`、`--compare` 或 `--import-fc` 的输入。补充报告的 `--redraw --numeric-source` 仍需核验数值来源是否与所选系列和选项一致；核验失败时不自动补算。新流程不自动搬移、删除或改写旧输出；旧清单含绝对路径，不宜直接移动。预处理继续使用 `paths.processed` 保存处理后信号，不复制这些信号到每个批次。
