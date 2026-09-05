# ALBA 第九步 v3 — 从 v2 继续的闭环数值验证

主报告：`REPORT_v3.md`。主模型：`stage9b_port/alba9b_full_loop.slx`。
本包的完整参考模型已实际运行于 MATLAB/Simulink R2025b，不是将 Python 输出冒充 Simulink。
15 mL/s 定位为 **Upper-feasible / near-limit case under nominal model assumptions**。
所有结论仅限条件模型；未验证硬件，也没有进行第十步。

## 打开与重建

在 MATLAB 将当前目录设为本工程根目录：

```matlab
addpath(pwd,fullfile(pwd,'stage9b_port'));
[~,S_vec]=alba9_scenario(300,2,'rise_s',3,'terminal_manager',1,'stop_remaining_mL',0.05);
V0_m3=0.0005;
open_system(fullfile('stage9b_port','alba9b_full_loop.slx'));
% 运行单例，完整 MAT 存到 results/manual_run（不要重用交付结果标签）：
run_v3_case(300,2,'manual_run','rise_s',3,'terminal_manager',1,'stop_remaining_mL',0.05);
```

完整一键重建/运行：先复制整个工程到新目录，安装 Python 依赖，再在 MATLAB 调用：

```text
python -m pip install --target .python_deps -r requirements-v3.txt
```

```matlab
rebuild_v3('D:\conda\python.exe') % 改为本机 Python 路径；MATLAB/Simulink 必须可用
```

依赖验证环境：Python 3.13.11、NumPy 2.5.2、SciPy 1.18.1、Matplotlib 3.11.1。
Python 第三方二进制库未打入 zip；全部项目依赖源码、冻结参数、参考模型源码均包含。
本次实际运行命令与版本见 `logs`、每例 `*.run.log`；一键入口用于复现，不表示又额外执行过一遍所有历史试验。

## 文件索引

- `results/D_diagnostic_baseline`：保留 Kp=50、Ki=600、Kaw=10，2 s 启动，修诊断后的三种控制器。
- `results/D_convergence_tight`：完整模型收紧容差/步长复跑。
- `results/A_startup_*`：只调整启动时间。
- `results/B_terminal_*`：只调整终止管理；固定步频管理变体仍非 PI。
- `results/B_guard_localised`：相同 0.25 mL 失败案例的零交叉精确定位。
- `results/C_combined_300`：先组合验证，再冻结选择。
- `results/FINAL_verified`：同一控制器的 200/300/400 mL，含尾流结束精确定位；`FINAL_verified_tight` 为最终方案收敛检查。
- `results/FINAL_selected`：保留定位修正前记录；其中 200 mL 的仿真结束时间对照未通过 2 ms 判据，见报告。不是最终验收目录。
- `validation/python_trajectories`：独立 Python 原始 MAT、同时间对齐误差 MAT、逐阶段及事件时间对照 JSON。
- `validation/ROM_full_*`：未删除启动段的 ROM—完整模型偏差，不作 ROM 全轨迹准确性的证明。
- `validation/frozen_parameter_diff.*`、`provenance_v2`：原精度参数和 ROM 系数审计。
- `poster`：最终方案流量/累计收集量 PNG、SVG、PDF。
- `frozen_reference`、`frozen_rom`：冻结参考源码/已辨识 ROM，未重新辨识。
- `baseline_v2_reports` 和根目录旧 `REPORT_A/REPORT_B`、`rebuild_all.m`、`alba9_closed_loop.slx`：v2 历史资料/ROM 入口，**不是 v3 主结论/入口**。

完整原始数据以每例 MAT 为准，包括连续日志、Simulink SimulationOutput、状态与最终运行点；CSV 为附加的抽样视图。
失败案例没有补写尾流或目标体积。`debug_*` 是保留的诊断开发过程，不纳入有效比较。
关闭命令不等于断电；实际速度、力与转矩不会被人工清零。
