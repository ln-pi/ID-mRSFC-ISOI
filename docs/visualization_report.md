# report 图表与恢复

所有主分析可视化集中在 report 阶段，入口为 `python run.py --stage report`。该命令默认处理 Aseries 和 Nseries；加 -c configs/Aseries.yaml 可仅处理一组。

默认使用4个绘图进程，可通过 `--plot-workers 1` 改为串行或指定其他正整数。进程数是执行选项，不改变配置标签或分析结果；show=true 时自动串行显示。

| 输出位置 | 图表 |
| --- | --- |
| Correlation/Figures/ | 每条有效记录的 FC 热图，兼容原有图表缓存 |
| report/figures/identification/<模式>/<周对>/<atlas>/ | similarity、paired、permutation、phi/dp circular、phi/dp matrix（系统间选中边比例） |
| report/figures/identification/<模式>/<周对>/ | relative_effect、group_individual 箱线图 |
| report/figures/statistics/<atlas>/ | 起始周/相邻周/窗口的自身与他鼠均值误差曲线；all_pairs 准确率与优势矩阵 |
| report/figures/statistics/structure/ | DD、ND 箱线图和结构空间散点图 |
| report/tables/structure_metrics.csv | 基于已保存相似度矩阵的 DD/ND 描述指标 |
| report/figures/reliability/<atlas>/ | 同周 Run 识别准确率概览 |
| report/figures/statistics/trajectories.png | 全部有效观测的逐鼠变化轨迹 |
| report/figure_index.json | 图片相对路径、类型、绘图签名和文件 SHA256 |

目录按比较模式和周对组织，例如 `observation_start/00w_to_01w/SM_VIS_AUD/`；windows 使用 `windows/<atlas>/`。不追加哈希后缀，不自动改写名称；非法路径或重名会明确报错。起始周、相邻周、全部周对分别保留图表；其重复周对不表示独立证据。Edge.save_mat: false 可关闭 Edge 矩阵图；ID.shuffle: false 的结果不生成置换图。

report 校验 03/05 清单的配置标签、FC 输入内容、数值输出哈希，然后使用已经保存的相似度、逐鼠指标、置换分布与 Phi/DP。绘图代码变化不要求重算识别或置换。数值来源缺失或不匹配时明确报错，需重新执行相关计算阶段。

report 阶段清单纳入全部图片。缺图或文件内容损坏会使 report 重新进入执行；逐图索引复用内容与签名均匹配的文件，恢复其他图表。中断时索引保留已完成项；未写入索引的最后少量图片可能重新生成。visualization.defaults.save/show 决定保存与显示；均为 false 时索引标记 disabled，report 仅输出表格与摘要。

图表沿用 id_*、edge_*、stat_* 样式。相似度颜色默认随数据范围变化，可设置 id_similarity_matrix.vmin/vmax 固定范围，避免旧固定范围裁剪新结果。绘图不自动执行 Wilcoxon 或配对 t 检验，也不自动标星。置换图显示已有 Monte Carlo p 值，未添加新的多重比较校正。

均值误差曲线沿用配置中的 SEM/SD 和相对效应梯度选项。梯度按有序比较的位置计算，并非每天的变化率；单一窗口没有梯度，改画相对效应。all_pairs 使用周×周下三角概览（对角格留空），不将无顺序含义的周对串成时间曲线。箱线图中的 Other/Self 对应他鼠与自身相似度，不是 A/N 生物学组间比较。Phi/DP 仅表示探索性识别贡献。

Edge 矩阵图保留旧圆图配套矩阵的定义：按 ROI 名称前缀聚合系统，用严格大于阈值的有符号边计数，除以两倍选中边数，显示下三角。它不是原始 ROI×ROI 的 Phi/DP 数值热图。

DD/ND 和结构空间图按唯一观察周/窗口对去重，沿用旧算法的最小值平移后对角占比与标准化距对角线距离定义。ND 依赖鼠排列顺序，因此仅作描述；重复周对或共享鼠不作为独立统计样本。图表使用 supplement 中的 atlas、组名和排除配置，均不自动检验。原历史均分算法仍保留在 --stage legacy。

周对统计矩阵仅显示下三角，行表示第二周、列表示第一周；对角格保留边界但不着色。DP/Phi 边比例矩阵保留下三角及其对角值。两类图都使用旧版阶梯形边界和单元格格线，上三角无色块与格线，不保留完整方形外框。识别相似度矩阵仍完整显示。
