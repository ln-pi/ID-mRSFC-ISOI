# 多配置自动流程

统一使用项目根目录 run.py。已移除 scripts/run_series.py、scripts/run_series.ps1 及编号脚本。

```powershell
python F:\Projects\ID-mRSFC-ISOI_refined\run.py
python F:\Projects\ID-mRSFC-ISOI_refined\run.py --dry-run
python F:\Projects\ID-mRSFC-ISOI_refined\run.py --start reliability
python F:\Projects\ID-mRSFC-ISOI_refined\run.py --stage report
python F:\Projects\ID-mRSFC-ISOI_refined\run.py --configs configs/Aseries.yaml configs/Nseries.yaml --stage report
```

使用当前 Python 环境运行，不再由 PowerShell 自动选择 Conda 环境。默认依次处理 Aseries、Nseries，每阶段独立子进程，入口会切换到项目根目录。一组失败后停止该组下游，继续另一组；没有表型数据时 phenotype 跳过。

每次启动在 reports/automatic/<时间戳>_<标识>/batch_summary.json 保存执行摘要，不在 reports/ 下保存阶段日志。可用 --report-root 指定其他摘要目录。

日志使用现有 logger 组件及各配置的 paths.log：Aseries 的所有阶段依次追加到 logs/Aseries.log，Nseries 则追加到 logs/Nseries.log。每次启动写入运行分隔标题，保留历史内容；批量流程不按阶段拆分或轮转日志。子进程的日志、普通输出与异常统一由主进程写入该组日志，不重复创建日志文件。

PowerShell 中显示 `<阶段> ... [进度条]`，按当前工作环节显示实际完成数量；初始化、校验和保存时没有确定总数则显示活动指示与耗时。阶段结束显示 complete、completed_with_warnings、skipped、failed 等状态。人工配准输入提示会直接显示，输入期间暂停进度刷新。全部结束或中断时显示项目状态及各组日志路径。

缺失输入等警告仍汇总为 completed_with_warnings；退出码 0 表示完成（可能有警告），1 表示失败或阻塞，130 表示中断。

--dry-run 保存计划并核对输入。对于 report，还核对已有数值结果及其来源，不计算识别或生成图片。首次分析从 fc 开始；已有完整结果仅需 --stage report。绘图中断后再次执行同一命令，可复用逐图校验通过的输出。

FC 输入继续采用 <mouse>_<run>_<week>_FiltDRR.npy，week/mouse/run 三层目录可交换顺序。缺失输入记录到 FC 报告，完整队列仍按 strict/cover 选择；轨迹保留不完整鼠的有效观测。

