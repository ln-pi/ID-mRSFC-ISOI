# 新配置说明

原配置的路径、鼠/周/Run 标识、atlas 排除规则和绘图属性继续保留。没有 None 分类；AUDv 的排除沿用三个 atlas 配置。新主分析属性如下：

预处理实际带通参数统一到 `preprocess.signal`：input_fps=30、output_fps=5、highpass=0.009、lowpass=0.08，默认与旧算法的硬编码完全相同。旧 `preprocess.filter/bin/robust` 曾未参与该 dRR 函数，迁移配置删除这些不生效属性，避免误以为修改它们会改变滤波。空间采样与配准仍沿用 bin4 图谱；若改用 bin1/bin2，需要同步匹配图谱、数据与掩膜尺寸，不在这次数值保持的迁移中自动切换。

```yaml
ID:
  completeness_mode: strict  # 或 cover
  session_seed: 42          # 置换的随机种子，纳入完整配置哈希
  shuffle: true
  n_iter: 1000
longitudinal:
  comparisons: [observation_start, adjacent, all_pairs]
  early_weeks: []
  late_weeks: []
  control_representatives: null
study:
  group: unspecified        # 明确设为 model/control 等，不能只根据鼠名猜组别
  acquisition_dates: {}     # 例如 {A01: {00w: '2026-01-01'}}
  missing_reasons: {}       # 例如 {'A01/01w/Run01': technical_failure}
phenotype:
  observations: null
  events: null
  tolerance_days: 7
```

`observation_start`：固定完整队列起始周与后续周；`adjacent`：所选连续周中的相邻周；`all_pairs`：全部周对。可加入 `windows`，并指定互不重叠的 `early_weeks` 和 `late_weeks`；窗口不在所选队列时标记不可计算，不能静默删周后换问题。

`control_representatives` 可显式指向对照组的 `Weekly/representatives.pkl`，必须来自可信项目输出；只使用同周、完整的对照观测，ROI 顺序必须一致。对照鼠自身进行分析时排除自己，输出有效对照数量。未配置时不生成对照偏离指标。

表型模板在 `configs/phenotype_template.csv`，字段 mouse/date/measure/value；事件模板在 `configs/events_template.csv`，字段 mouse/entry_date/event_date/censored，censored 为 0 或 1，删失个体 event_date 表示最后观察日。表型按 mouse、测量名取采集日之前最近、容差范围内的测量，避免使用未来数据。真实采集日缺失时不进行伪日期匹配。事件接口仅验证并输出数据，不自动拟合模型。

当前主流程不随机分周，历史分半入口和实现已移除。配置文件保持原样：`analysis.split_axis`、`analysis.selection_mode`、`ID.group_mode/repeat_id`、`Edge.seed/repeat` 不再控制任何历史分析。为保持既有结果位置，输出标签仍使用原配置哈希和命名规则（包括 group_mode、repeat_id、Edge.seed），这些名称不表示主流程执行随机分半或 Edge 重复采样。`ID.session_seed` 继续控制置换随机数；Phi/DP 使用固定周对。`analysis.statistics_file` 及 supplement 中当前报告使用的绘图设置保留。

本次代码清理不迁移已有数据和结果目录。阶段清单包含代码哈希，因此再次执行数值阶段时可能重新计算；report 仍按已保存数值的配置、输入及输出哈希校验来源。
