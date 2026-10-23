# 配置与数组契约

配置使用 YAML；相对路径以 main 根目录解析，绝对路径保持绝对路径语义。`configs/default.yaml` 是可用的 A 系列示例，不再是全 null 的占位模板。Aseries/Nseries 保留当前 10 鼠配置名单，实际分析鼠数由审计决定。

## 配置字段

| 字段 | 含义 |
| --- | --- |
| paths.raw / processed | 每次采集的原始/处理后数据根目录 |
| paths.mask | ROI PNG 与 Chamber_Mask.png；默认 ../data/mask |
| paths.output / log | 独立结果及日志位置 |
| variables.dataset | 路径附加目录 ID，可为空 |
| variables.Mouse_ids / Week_ids / Run_ids | 显式、有序、唯一的采集维度标签 |
| variables.Atlas | ROI 排除集合列表；每个分区独立 edge store |
| preprocess.keepna / resume | 旧读取行为参数；逐步骤校验续跑 |
| correlation.resume | 校验 FC/图缓存后复用 |
| analysis.split_axis | combined 的 week 或 run 拆分轴 |
| analysis.selection_mode | combined 对另一轴采用 single 或 cumulative |
| analysis.statistics_file | combined 聚合结果文件名，默认 statistics.pkl |
| ID.completeness_mode | strict 全配置完整队列，或 cover 最大完整鼠×连续周窗口；只影响 combined |
| ID.group_mode / session_seed | sequential/random、不重叠拆分及种子 |
| ID.shuffle / n_iter | combined 身份重标记；包含并列的上尾 (b+1)/(B+1) |
| Edge.seed | 保留的完整配置身份字段；固定 ID 配对使用 session_seed，不另做随机 edge 拆分 |
| Edge.threshold / gap / cross_only / save_mat | edge 图参数 |
| pairwise.phases | 阶段名称→周标签，仅用于完整数组之后的选择 |
| pairwise.phase_comparisons | 两个不重叠阶段的有向比较 |
| pairwise.week_positions | 可选周数值位置，必须覆盖全部 Week_ids 且严格递增；默认解析 00w 等数字标签 |
| pairwise.window_interpretation | N 疾病时间窗或 A 对齐观察窗的解释记录 |
| descriptive.reference_week | 自身变化参照；默认 03w |
| visualization | FC、combined 与统计图参数；pairwise 报告统一 r 色阶 [-1,1] |

完整 parsed config 仍参与输出哈希，绘图参数也会改变结果目录。没有把源预处理身份、分析身份、绘图身份静默拆为新缓存体系。根入口接受 `--configs` 顺序批处理，默认 Aseries 和 Nseries；`--import-fc` 要求单一 `-c` 并从 edges 开始。`--report-section` 是运行选择参数，不修改配置身份，可单独重绘 pairwise/combined/statistics；`--report-root` 指定批次状态记录目录。

已移除 Session_ids/Trial_ids、数值 group_mode 12/13/14、非空 repeat_id、supplement/DD/ND，以及重复固定配对的 Edge.repeat。旧数字脚本入口改为根 run.py/run.ps1；不对旧 YAML 进行隐式转换。

## FC / edge store

`correlation/<week>_correlation_matrices_<atlas>_excluded.pkl`：行是配置鼠，列是 Run，每个有效单元格为带一致 ROI 行列标签的 Pearson FC 方阵。原始表保持导入字节，以便来源复核和 FC 往返。

`results/edges/<atlas>.npz`，schema=`fc_edges_v1`：

| 键 | 形状/口径 |
| --- | --- |
| F_r | N×W×R×E，原始 Pearson FC 唯一边 |
| F_z | 同形，atanh(clip(F_r,-0.999999,0.999999)) |
| valid_record | N×W×R；表示 FC 有效，不等于模式相关可定义 |
| invalid_reason | N×W×R；缺失、审计失败、ROI 顺序/FC 无效的原因 |
| edge_index | 2×E，numpy.triu_indices(K,k=1) |
| metadata | JSON Unicode scalar；维度顺序、ID、ROI、schema、参数、源 FC/audit SHA256、源配置、Git commit/dirty 状态、代码逐文件哈希、包版本 |

全常量边可以是有效 FC，但无法定义模式 Pearson；其原因在派生 S 中标记。不根据阶段、Run 或完整 14 周预先筛掉动物。

每个 NPZ 原子写入，无 pickle 对象。JSON sidecar 包含数组文件 SHA256；后续步骤同时核对文件、源 FC/audit、配置和代码哈希。源或计算代码变化后显式重建 edges，不能静默复用旧数组。sidecar 的非原子双文件更新时间差最多导致明确校验失败，不会把不一致文件视为有效缓存。

## 每 Run 的 S

`results/pairwise/<atlas>/<run_id>.npz`，schema=`pairwise_similarity_v1`：

| 键 | 形状/口径 |
| --- | --- |
| S_r | N×N×W×W；Pearson_over_edges(F_z[i,u,run],F_z[j,v,run]) |
| S_z | 同形；atanh 仅用 nextafter(1,0) 处理数值边界，无 ±4 展示截断 |
| valid_similarity | 同形；模式 Pearson 可定义，包含有效自比较 |
| statistical_mask | 同形；另排除 i=j 且 u=v 的同一采集自比较 |
| invalid_reason | 同形 uint8；0=定义有效，1=端点无效，2=自比较排除 |
| pattern_reason | N×W；端点缺失/审计失败/常量边原因 |
| metadata | Run ID、分区、鼠/周顺序、尺度、算法、edge store 校验值和代码来源 |

无效位置 S_r/S_z=NaN；同一采集自比较 S_r=1、S_z=NaN。不同采集恰好相关 ±1 时使用 nextafter 边界得到有限 z，并明确记录此规则。周对面板行=i@u、列=j@v，反向面板为转置。

```python
from srcs.analysis.connectome import EdgeStore, load_arrays, SIM_SCHEMA
store = EdgeStore.load("/absolute/output/results/edges/AUDv.npz")
meta, values = load_arrays("/absolute/output/results/pairwise/AUDv/Run01.npz", SIM_SCHEMA)
panel = values["S_r"][:, :, 4, 7]       # 04w → 07w (按保存标签核对位置)
trajectory = values["S_r"][0, 0, :, :]  # 第一只鼠的跨周自身矩阵
```

常规阶段入口会额外校验来源和缓存，直接数组加载只做 schema/维度读取，适合明确知道来源的检查。

combined 原生接口为 `EdgeStore.average_halves(...)`，返回两个 mouse×edge 的 F_z 均值数组，再传给 `srcs.analysis.identification.identify(first, second)`。Phi/DP 通过 `srcs.analysis.edge_contribution.edge_contributions(...)` 直接消费 mouse×week×run×edge 的 F_r。旧的嵌套 DataFrame 半组接口及 Identify/Edgewise_Contribution 类已移除。

## 选择与描述结果

pairwise 的 subjects_by_pair.csv 保留每周对、鼠、Run、窗口、cohort、候选 ID、周间隔、尺度、相似度、自身优势、百分位及识别结果。subjects_by_window.csv 先按周对计算，再逐鼠取均值；同时保留有效/请求周对数和候选数范围。所有窗口/无结果原因及完整来源在 selection.json。

两个策略并行：available_pair 取本周对两端共同有效鼠；fixed_window 取整个所选周窗口均有效鼠。少于两只候选时不伪造排名，记入选择审计。baseline 内部零周对，保留 no_cross_week_pairs 状态。默认阶段内周对计数：onset/early/late 各 3，middle 6，baseline 0；early→middle 为 3×4 个有向周对。

本版没有预设更宽松的“阶段最低有效周数”、单 Run 替代或匹配间隔推断。CSV 保存实际标签对应的周数间隔，供后续相邻周/匹配间隔筛选。它不是实际采集日差；本次输入没有新增采集日期或失访原因资料。

## 幅度与拓扑（descriptive v2）

不改配置即可执行 `python run.py --stage descriptive`；数值结果通过 `--stage report --report-section descriptive` 独立绘图。代码更新后，先执行 `--stage edges` 重建代码身份已变化的缓存。

可选配置（省略时采用下列默认策略；weeks 默认为全部配置周）：

```yaml
descriptive:
  reference_week: 03w
  run_views: [per_run, mean_fc]
  # weeks: [03w, 04w, 05w]
  # network_mapping:
  #   PF-PL_L: association
  #   ... 每个实际使用的 ROI 都必须显式映射
```

`per_run` 保留各 Run；`mean_fc` 在 Fisher-z FC 边空间平均所有配置 Run，再计算指标，要求所有 Run 均有效。不得将平均 FC 的距离当作各 Run 距离的平均。`weeks` 按基础配置时间顺序选择分析窗口，不改变 edge store；自身参照周可以位于所选窗口之外。

两套队列：`available_week` 使用本周/本视图有效鼠；`fixed_window` 使用选定所有周均有效的鼠，按各视图分别确定。均保留完整鼠轴，以掩码/NaN 表达不参与位置。LOO 参照名单逐鼠记录，不混合 A/N 或其他不同配置。

定义（f 为 Fisher-z FC 唯一边向量）：

| 指标 | 公式/含义 |
| --- | --- |
| Euclidean | sqrt(sum_e((f_i−f_j)^2)) |
| RMS | Euclidean / sqrt(E) |
| correlation distance | 1−Pearson_over_edges(f_i,f_j)，常量模式不可定义 |
| residual_group_loo | f_i−mean(f_j,j≠i)，同队列、同周、同视图 |
| individual deviation | RMS(residual_group_loo) |
| group dispersion | mean_i(mean_e((f_i−全组均值)^2))，分母 N |
| edge_variance | var_i(f,ddof=1)，分母 N−1 |
| delta_from_reference | f_i,t−f_i,reference，独立于群体 LOO ΔFC |

N<2 时离散度、样本方差、LOO 偏离均标记未定义，不用单鼠制造零群体差异；单鼠自身欧式/RMS 距离仍为 0。缺失记录不填零，常量模式仍可计算幅度距离。距离矩阵主对角线不参与唯一成对距离汇总。

输出目录 `results/descriptive/v2/`，schema=`descriptive_magnitude_topology_v2`：

- `<atlas>/<cohort>.npz`：`distance_euclidean/rms/correlation[W,V,N,N]`；`residual_group_loo/delta_from_reference[N,W,V,E]`；`edge_variance/group_mean[W,V,E]`；`valid_record[N,W,V]`、`n_valid_mice[W,V]`、`edge_index`。V 的 ID、Run 汇总规则和网络映射保存在 metadata。
- `distances.csv`：每周/视图/队列的唯一非自身鼠对及三种距离，常量模式原因。
- `subjects.csv`：LOO RMS/欧式/平均有符号偏离、正负边比例、偏离度降序排名、参照鼠名单及自身参照周变化。
- `dispersion.csv`：群体离散度、平均逐边方差、有效鼠数、唯一鼠对数和状态。
- `regional.csv`：逐鼠 ΔFC 和群体逐边方差的 ROI、network_pair、hemisphere_pair 汇总。
- `edges_<atlas>.csv`：边到两端 ROI、网络、半球的映射。
- `examples.csv`：每周/视图/队列中 LOO RMS 最大的描述性示例；保存全部精确并列名单，绘图按配置鼠顺序破同分。不实施离群显著性检验或自动剔鼠。
- `reliability.csv`：所有同日 Run 对的原始 FC 重复性。
- `selection.json`、`methods.json`、`index.json`：实际队列、公式、来源、数据哈希和数值产品索引。

ROI 贡献将每条边的平方偏离能量（或边方差）各分一半给两端，再除以全体唯一边总能量；各 ROI 贡献之和为 1。网络对/半球对按无序对归属唯一边，各层贡献也分别为 1。另报边数和每边均值，避免边多的模块被误判为每边差异大。总能量为 0 时贡献比例为 NaN。

ΔFC 正值表示相对参照 FC 更强，负值表示更弱；不解释为兴奋/抑制。有符号均值与正/负方向 RMS、正负边比例同时保存，避免方向抵消。方差是无方向量，不生成“负方差”。默认按 atlas 已定义的 AUD/LAT/MED/PF/SM/VIS 解剖分组；这不是独立估计的功能网络。自定义 `network_mapping` 必须覆盖所用 ROI；未知 atlas 标签无映射时只提供可确定的 ROI 汇总。

`report/descriptive/` 生成每视图的两队列离散度/样本数、逐鼠偏离轨迹、首个有效周的三种距离矩阵与全部小鼠 ΔFC 图、全周方差图、ROI/网络方差和示例方向图。所有周的全部鼠与两队列数值保存在上述数组/表中；图表使用保存产品，核查哈希后独立重绘。图中矩阵对角线不代表自连接，保持遮罩。

旧 `results/descriptive/topology_<atlas>.npz` 等 v1 文件不作为 v2 输入，不会自动删除。v2 的 `residual_group_loo` 才是群体 ΔFC；旧 `delta_from_reference` 的自身参照含义不变。`descriptive` 阶段不估计置信区间、ICC、随机效应或组间推断；稳定性与纵向模型在下面的独立阶段计算。

## 稳定性与纵向演化

`stability` 和 `longitudinal` 读取校验后的 descriptive v2 与 edge store；共同采用 `descriptive.weeks` 选择的窗口。时间位置沿用 `pairwise.week_positions`，默认从 `00w` 等数字周标签提取，非数字标签必须显式配置。位置是数值时间单位，不自动等于实际采集日期。

```yaml
stability:
  n_bootstrap: 200
  seed: 42
longitudinal:
  n_bootstrap: 200
  seed: 42
  model_network_pairs: within  # within=网络内均值；all=另加网络间均值
  maxiter: 300
```

这些都是省略时采用的默认值，不要求修改现有配置。bootstrap 次数可为 0（不生成区间）或至少 20，种子为非负整数。默认 200 次用于可重复的初步区间，正式研究可增加次数检查 Monte Carlo 精度。

### Run 重复性与 ΔFC 跨周相似

每个同日 Run 对使用两次 Run 均有效的共同鼠名单，在各 Run 内重新计算 LOO ΔFC，避免把两次不同的模板组成误当成 Run 差异。保存原始 FC 与 ΔFC 的模式 r、RMS、有符号差和 Run 均值—差值数据；这些是同日重复性，不能代替跨周稳定性。

`residual_similarity[i,view,u,v] = Pearson_edges(ΔFC[i,u,view], ΔFC[i,v,view])`。只保留可定义模式，缺失/常量模式为 NaN；同一记录主对角线始终排除。输出 dynamic `available_week` 和固定成员 `fixed_window` 两套结果。动态模板下相似度也可能受模板成员变化影响，不能与固定成员结果混称。

### Edgewise ICC 与稳定个体成分

逐边模型定义：

```text
y[i,t,r,e] = fixed_week[t,e] + fixed_Run[r,e]
             + u[i,e] + v[i,t,e] + error[i,t,r,e]
stable_ICC[e] = sigma_mouse² / (sigma_mouse² + sigma_mouse_week² + sigma_error²)
same_week_ICC[e] = (sigma_mouse² + sigma_mouse_week²) / total
mean_Run_stable_ICC[e] = sigma_mouse² / (sigma_mouse² + sigma_mouse_week² + sigma_error²/R)
```

这是调整分类周次、Run 固定效应后的单 Run 方差比例，不冒用未调整的传统 ICC 类型名称。`u` 是时间不变的鼠随机截距；`v` 是鼠×周共同项，也可能包含同日成像状态；Run 残差不被直接命名为纯技术噪声。

估计器使用所选周与所有 Run 均完整的平衡固定队列，至少 3 鼠、2 周、2 Run。不会为了凑齐平衡数据插补。无法支持的设计明确写 `unsupported`，方差/ICC 为 NaN；动态可用观测仍用于其他分析。实际保留/排除名单和缺失原因记录在 cohorts.csv。

在平衡设计下，扣除固定周次和 Run 后存在三个正交限制对比子空间：

| 子空间 | 自由度 | 特征方差 |
| --- | --- | --- |
| cell 内 Run 残差，扣除 Run 固定主效应 | (N×W−1)(R−1) | λ0=σerror² |
| mouse×week 交互 | (N−1)(W−1) | λ1=σerror²+Rσmouse_week² |
| mouse | N−1 | λ2=σerror²+Rσmouse_week²+R×Wσmouse² |

相应平方和分别来自 `y−cell_mean−Run_mean+grand_mean`、`cell_mean−mouse_mean−week_mean+grand_mean` 及 `mouse_mean−grand_mean`（后两者乘以 R 和 R×W）。REML 限制似然为 `0.5*sum_k[df_k*log(λk)+SS_k/λk]`，约束 `0≤λ0≤λ1≤λ2`。若无约束均方违反顺序，则按自由度加权合并相邻子空间，再反解非负方差分量。这是包含零边界的约束 REML，不是把负矩估计逐项截成零。独立检查将其与完整观测协方差矩阵的直接 REML 数值优化比较。

同时保存绝对方差、ICC、零边界状态、各子空间均方/自由度、群体周均值、Run 偏移、`stable_mouse_blup[N_complete,E]`、`mouse_week_blup[N_complete,W,E]` 与 `conditional_residual[N_complete,W,R,E]`。稳定成分是条件预测/收缩估计，不能等同于独立观测到的纯神经身份特征。

方差与 stable ICC 的 95% 区间使用完整动物轨迹 block bootstrap，每次重算模型；同一动物的周和 Run 不拆散。这是 percentile 区间，不替代零边界参数的正式假设检验。网络汇总分别命名 `mean_edge_stable_icc` 和 `ratio_of_summed_variances`，后者为该网络各边稳定方差之和除以总方差之和；两者不混用。

### 纵向轨迹、收敛/分化与随机斜率

逐周群体离散度沿用描述阶段公式；对整只动物的时间轨迹重采样，每次重新计算本周群体中心和离散度，保存点估计及 percentile 95% 区间。逐鼠 LOO RMS、有符号偏离及排名保持（共同鼠之间的 Spearman r）同时保存，不对共享周对生成伪独立 p 值。

对离散度随数值时间的整体 OLS 斜率使用相同 bootstrap。95% 区间完全大于/小于 0 时，分别标记 increasing/decreasing_dispersion，否则为 inconclusive。该标签描述整体线性趋势，不证明全程单调或疾病导致分化/收敛；必须结合可用/固定队列、Run 重复性与独立 QC。尚未补齐的临床/QC 变量不会被软件假装调整。

网络随机斜率模型以每次采集的**原始 Fisher-z FC 网络均值**为响应，不把共享 LOO 模板的偏离分数作为独立响应：

```text
network_mean[i,t,r] = beta0 + beta_time*time + fixed_Run[r]
                      + u0[i] + u1[i]*time + v[i,t] + error[i,t,r]
```

默认拟合 all_edges 与各解剖/自定义网络内部均值；all 可加入网络间均值。使用 [Statsmodels MixedLM](https://www.statsmodels.org/stable/generated/statsmodels.regression.mixed_linear_model.MixedLM.html) 的 REML，鼠随机截距和斜率允许相关，另含鼠×周随机截距。使用全部有效 Run 观测，至少 4 鼠、每鼠至少 3 周且存在足够同周 Run 重复。时间在拟合时中心化/标准化以改善计算，输出斜率和协方差还原到原时间单位。

保存固定斜率、随机截距/斜率方差及协方差、鼠×周与残差方差、条件逐鼠斜率、拟合值/残差、优化器尝试、收敛、边界、Hessian/奇异警告。未收敛时不生成个体 BLUP 斜率；不支持/全优化器失败有明确状态。线性模型不自动证明合理，需检查残差随时间的结构和方差不齐；必要时另行预设分类时间/分段模型。死亡相关缺失不会因使用全部有效观测或 REML 自动得到校正。

### 输出与运行

```powershell
python run.py --stage edges
python run.py --stage descriptive
python run.py --stage stability
python run.py --stage longitudinal
python run.py --stage report --report-section stability
python run.py --stage report --report-section longitudinal
```

| 位置 | 内容 |
| --- | --- |
| results/stability/repeatability.csv | 原始/共同模板 ΔFC 的 Run 重复性与均值—差值 |
| results/stability/residual_pairs.csv | 逐鼠、视图、队列的唯一跨周 ΔFC 相似 |
| results/stability/edge_icc.csv / regional_icc.csv | 逐边方差/ICC/区间、区域汇总 |
| results/stability/cohorts.csv | 平衡模型的保留/排除原因 |
| results/stability/&lt;atlas&gt;/variance_components.npz | 方差、ICC、稳定/鼠×周 BLUP、残差、设计标签 |
| results/stability/&lt;atlas&gt;/residual_&lt;cohort&gt;.npz | N×V×W×W 残差相似矩阵及有效掩码 |
| results/longitudinal/dispersion.csv / trends.csv | 逐周区间、整体斜率区间和方向标签 |
| results/longitudinal/subjects.csv / rank_pairs.csv | 逐鼠轨迹、描述性斜率、共同鼠排名保持 |
| results/longitudinal/random_slopes.csv / mouse_slopes.csv | 网络模型、条件逐鼠斜率 |
| results/longitudinal/models.json / model_diagnostics.csv | 结构化模型诊断与逐记录拟合值/残差 |
| results/longitudinal/&lt;atlas&gt;/trajectories_&lt;cohort&gt;.npz | W×V 离散度、2×W×V 区间、N×W×V 偏离度、群体 FC/相邻周变化 |
| report/stability / report/longitudinal | PNG/PDF、索引；无需源 FC/edge cache 的独立重绘 |

两阶段各有 methods/index.json，保存 schema、数值输入 SHA256、完整规则、代码和 Statsmodels 等依赖版本。它们与幅度/拓扑 v2 分开保存，不改写旧指标的含义。
