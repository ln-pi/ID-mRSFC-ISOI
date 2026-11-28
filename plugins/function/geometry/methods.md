# FC 高维几何分析

每个有效采集在统一的 ROI 顺序上取 FC 严格上三角并做有限截断的 Fisher-z 变换。几何阶段直接读取边数组中的 F_z。每对 Run 分别计算，不平均 Run，也不混合周次。available_week 使用当周两个 Run 都有效的鼠，fixed_window 要求整个窗口内均有效。

- 幅度距离是逐边差平方均值的平方根。拓扑距离是 1 减边向量 Pearson 相关，常数向量的拓扑距离记为缺失。
- 同周、同 Run 的分散度是逐边样本方差的均值，分母使用 n−1。跨 Run 距离矩阵的同鼠对角线保留真实重复测量距离。
- 每只鼠在两个查询方向分别计算 q 和分离余量。q 是自身比异鼠更近的比例，平局计半分。余量为最近异鼠距离减自身距离。保存完整候选、固定种子的匹配候选和逐鼠留一结果。
- 交叉验证平方差异按两个 Run 的异鼠差向量内积除以边数计算，保存恒等权重、Ledoit–Wolf 噪声权重和留出被比较两鼠重估权重三种视图。负值保留。
- inference 在每周先按鼠平均查询方向，再平均动物。留一敏感性范围取第 5 和 95 百分位。趋势按整鼠轨迹重采样，结果保存于 trends.csv。
- 幅度、形状、可重复差异与识别是不同的指标。PCA 在每分区、队列和 Run 对的全部可用点上一次拟合，用作描述。数值坐标与解释方差随 points.csv 保存，重绘读取已存坐标。

`results/geometry/` 保存边基准、距离、跨周位移、分散度、逐鼠候选比较、交叉验证差异、投影坐标和输入可校验的产品索引；`results/inference/` 保存逐周主终点和趋势。独立 compare 先核对分区、ROI、边顺序及采集轴。A/N 分组抽取整鼠轨迹，GSR/noGSR 检查两侧有效成员并同步抽取同鼠。比较区间标为 sensitivity_range。


## 启用与历史产物

在原生 FC 或 SeedFC 配置中设置 `plugins.function.geometry.enabled: true`，也可通过 `--function geometry` 为当前命令启用。几何数值阶段在 edges 后执行，`inference` 使用已发布的几何结果。A/N 合并展示与 Geometry 比较共享批次的 `compare/`：数值分别在 `results/core/` 和 `results/geometry/`，图分别在 `report/core/` 和 `report/geometry/`。

本地单独安装功能包可在已激活的 Python 环境运行 `python -m pip install <P3工作树>/plugins/function/geometry`。功能包通过 `plugins.function.geometry` 提供清单和实现，依赖 scikit-learn。原始结果保存在原有目录中。

```powershell
python run.py --configs configs/SeedFC_Aseries_noGSR.yaml configs/SeedFC_Nseries_noGSR.yaml --start seedfc --function geometry
python run.py --from "<已有输出目录>" --stage report --report-section geometry --function geometry
python run.py --compare "<A输出目录>" "<N输出目录>" --function geometry
```

历史重绘读取并核验已保存的几何产品。几何比较的每个序列需包含 `results/geometry/index.json` 和 `results/inference/index.json`。
