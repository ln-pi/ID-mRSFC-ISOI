# FC 高维几何分析

每个有效采集在统一的 ROI 顺序上取 FC 严格上三角并做有限截断的 Fisher-z 变换。几何阶段直接读取边数组中的 F_z。每对 Run 分别计算，不平均 Run，也不混合周次。available_week 使用当周两个 Run 都有效的鼠，fixed_window 要求整个窗口内均有效。

- 幅度距离是逐边差平方均值的平方根。拓扑距离是 1 减边向量 Pearson 相关，常数向量的拓扑距离记为缺失。
- 同周、同 Run 的分散度是逐边样本方差的均值，分母使用 n−1。跨 Run 距离矩阵的同鼠对角线保留真实重复测量距离。
- 每只鼠在两个查询方向分别计算 q 和分离余量。q 是自身比异鼠更近的比例，平局计半分。余量为最近异鼠距离减自身距离。保存完整候选、固定种子的匹配候选和逐鼠留一结果。
- 交叉验证平方差异按两个 Run 的异鼠差向量内积除以边数计算，保存恒等权重、Ledoit–Wolf 噪声权重和留出被比较两鼠重估权重三种视图。负值保留。
- inference 在每周先按鼠平均查询方向，再平均动物。留一敏感性范围取第 5 和 95 百分位。趋势按整鼠轨迹重采样，结果保存于 trends.csv。
- 幅度、形状、可重复差异与识别是不同的指标。PCA 在每分区、队列和 Run 对的全部可用点上一次拟合，用作描述。数值坐标与解释方差随 points.csv 保存，重绘读取已存坐标。

`results/geometry/` 保存边基准、距离、跨周位移、分散度、逐鼠候选比较、交叉验证差异、投影坐标和输入可校验的产品索引；`results/inference/` 保存逐周主终点和趋势。独立 compare 先核对分区、ROI、边顺序及采集轴。A/N 分组抽取整鼠轨迹，GSR/noGSR 检查两侧有效成员并同步抽取同鼠。比较区间标为 sensitivity_range。


## 四阶段汇总与合并展示

pairwise.summary_phases 指定 baseline、early、middle、late，周名单沿用 pairwise.phases。Run 独立保存，并额外计算所有有向 Run 对。results/pairwise/<atlas>/phase_similarity.npz 保存每个 Run 对、阶段对、队列和平均方法的 N×N 矩阵，以及每格有效周对数。correlate_mean_fc 先逐鼠在阶段内平均 Fisher-z FC 边，再做边模式 Pearson；mean_correlations 先计算各周对的模式 Pearson，再对 r 作算术平均。后者排除同鼠、同 Run、同周的同一采集自比较；baseline 同 Run 的自身对角线缺失，异鼠比较保留。前者同 Run 同阶段的自身对角线屏蔽。跨周自身相关仍保留，各矩阵格的有效周对数分别保存。

available_phase 要求各端点阶段在对应 Run 的全部周均有效；fixed_window 进一步要求四阶段窗口全部周在参与的 Run 中有效。缺失与常数模式保留为 NaN，计数随矩阵保存。

longitudinal 的 phase_subjects.csv 和 phase_dispersion.csv 汇总已算好的每周指标。阶段 LOO RMS 是每周 RMS 的均值，不先平均偏移向量。有符号偏离均值单独保存，允许方向抵消。每指标保留有效周数；阶段图是描述性均值，不平均逐周置信区间。

先运行各系列计算和 report，再使用 python run.py --compare <A输出目录> <N输出目录>。比较输出的 report/series_overlays 保存统一坐标中的每周与四阶段 trajectories、方向热图和跨 Run 阶段矩阵。A/N 使用固定颜色，Run 使用固定点形，模式相关矩阵统一 [-1,1]，方向 RMS 两侧共用非负 cbar。各组 LOO 模板分别计算。

批量运行两份 A/N 配置并完成全部 report 时，如果两组的几何、推断和描述产物齐全、分支相同，协调器会自动执行合并比较，并将 comparison 路径写入 batch_summary.json。也可用 --compare 独立重画。合并图位于各比较目录的 report/series_overlays，避免多次比较相互覆盖。
