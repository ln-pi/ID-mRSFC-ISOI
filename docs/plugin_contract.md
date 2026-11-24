# SeedFC 来源契约

noGSR 的两组入口分别为 `configs/SeedFC_Aseries_noGSR.yaml`（A25、A26、A27、A28、A29、A30、A31、A32、A33、A35）和 `configs/SeedFC_Nseries_noGSR.yaml`（B20、B22、B23、B24、B25、B30、B31、N12、N13、N15、N17、N20、N27、N29）。两份配置都读取 organized_original_ids，并分别写入系列日志。两组均直接使用原编号，按指定顺序保留 10 只对照鼠和 14 只疾病鼠。

```powershell
python run.py --configs configs/SeedFC_Aseries_noGSR.yaml configs/SeedFC_Nseries_noGSR.yaml --stage seedfc --dry-run
python run.py --configs configs/SeedFC_Aseries_noGSR.yaml configs/SeedFC_Nseries_noGSR.yaml --stage seedfc
```

第一条只检查，第二条生成两组原生 FC 周表和缺失审计。配置现采用单一 FourNetworks 分区，通过 partitions.exclude 排除 AUD、TEA、FRP、SSS 双侧节点，实际分析保留 14 个 ROI、91 条唯一边。源文件保持 22 个 seed。两个 noGSR 配置可用 --start seedfc 连续运行。

`configs/SeedFC_Aseries_noGSR.yaml` 从 `../data/SeedFC_001_008/organized_original_ids` 读取整理好的采集。源目录下按 branch、series、week、run、mouse 定位 `FC_matrices.npz`。读取 `r` 与 `seed_names`，在写周表前核对各采集的标签、顺序、矩阵和可用的身份字段。输入已包含的 `z` 不会再次变换成原始相关。

每个分区显式列出排除的 seed 名。同一分区需要至少三个保留的 seed。未知排除标签或跨记录标签不一致将终止转换。缺少文件或分区内存在非有限连接的记录保留在审计中，对应周表单元格为空。完整有限矩阵还须对称、相关系数在允许范围内、主对角线为 1。

FourNetworks 保留 Motor（MOP、MOS）、Somatosensory（SSP）、Visual（VIS）、DefaultMode（RSP、ACA、PL）各 ROI 的左右半球。排除发生在每个分区的 FC 有效性判定之前。按原编号清单核验，noGSR Aseries 有 277 个有效采集、3 个缺失文件，Nseries 有 392 个有效采集。该网络映射是预先指定的分组，图中的 ROI 与网络对按上述顺序排列。

`--stage seedfc --dry-run` 只核验源文件。`--start seedfc` 先生成带标签的原生周表、逐记录审计和来源摘要，再进入 edges。启用该来源的配置与原生 preprocess/fc 互斥。


## 插件目录与配置

`plugins/source/<id>/plugin.yaml` 定义 FC 来源。单次运行可选择一个来源，空配置走原生 FC。来源声明输入验证入口和参与数据身份的数值源码。 SeedFC 启用位置为 `plugins.source.seedfc.enabled`，旧版 `plugins.seedfc` 配置继续可读。

`plugins/function/<id>/plugin.yaml` 定义分析或绘图能力。每个功能单独启用，可声明数值阶段、依赖阶段、产物、报告分节和两系列比较入口，或登记批次收尾任务。功能阶段使用同一份已校验的 edge store 和已发布的结果。功能开关属于阶段配置，不进入 FC 或 edge store 数据身份。

```yaml
plugins:
  function:
    seedfc_supplement:
      enabled: true
```

这段配置合并到已有 SeedFC 配置中的 `plugins` 下，原有 `plugins.source.seedfc.options` 提供采集路径、系列和分区。

主流程从已启用插件的清单生成阶段列表。功能插件可随项目提供，也可以作为 Python 包安装到同一环境，包内须提供对应命名空间和 `plugin.yaml`。显式启用后，启动时校验功能包和专用依赖。普通运行、历史重绘和比较都可通过 `--function <id>` 为当前命令启用已安装功能。
