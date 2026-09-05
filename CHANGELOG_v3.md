# 本轮修改清单

原 v2 目录、用户压缩包、实验数据、固件均未覆盖。未做 SINDy/Fixed-ID；未增加 phi*x_dot；未开展第十步。

## D：只修程序与诊断

- `alba9_controller.m`：按总原始指令限幅前后记录饱和，保留原 PI 增益及限制。终止管理默认关闭。
- `alba9_scenario.m`：追加显式停止/诊断和可选终止参数，原前 19 项顺序不变。
- `stage9b_port/alba9b_plant_block.m`：只扩展有单位守卫、未截断重叠量和压力输出；原 RHS 不变。
- `stage9b_port/alba9b_log_index.m`：46 个日志字段，分别记录真实速度、加速度、守卫、估计体积等。
- `alba9v3_extend_model.m`、`stage9b_port/build_alba9b_model.m`：反馈估计体积来自独立 ∫Q；机械停止驻留计时；守卫和尾流结束零交叉定位；保留单次关阀/Q 复位；保存完整运行状态。
- 最后连线审计同时把参考调度及关阀阈值改接独立 ∫Q，不读取真实体积；理想模型下等价，仍实际复跑三组与严格容差结果，并自动断言连接来源。此前结果保留在 `debug_pre_wiring_audit`。
- `alba9v3_diagnose.m`、`v3_classify_status.m`：分开四类事件、显式完成分类、分阶段误差和有单位守卫，闭阀后运动/压力与独立守恒检查。
- `v3_unit_tests.m`、`v3_verify_points.m`：15 项单元测试、6117 个逐点状态及边界邻域测试，包含原 RHS/Jacobian/守卫。
- `run_v3_case.m`、`v3_prepare.m`：真实 Simulink 编译/运行，记录环境、许可证、完整 MAT、状态、日志与警告。

## A：启动优化单独保存

- 2 s 基线不变，对纯 PI/前馈＋PI 分别增加 2.5 s、3 s 参考。
- 增益、物理对象、驱动限制不变；同时检查自然间隙及重新接触，未人为消除起步延迟。

## B：终止优化单独保存

- 固定步频原基线不动；另设清楚命名的剩余体积终止管理变体，不加入 PI。
- 比较距离关阀阈值 0.05 / 0.25 mL 时停止收线请求，仍受原开阀斜率和闭阀制动限制。
- 两项 0.25 mL PI 失败保留；守卫精确定位后依然失败，原粗定位结果另保留。

## C / FINAL：组合与独立验证

- 300 mL 组合通过后冻结：纯 PI，3 s，0.05 mL；同一参数用于 200/400 mL。
- `FINAL_selected` 保留首次最终试验，200 mL 的仿真完成时刻差 2.924 ms 超过 2 ms 判据。
- 修正尾流结束零交叉，不放宽判据。重跑结果为 `FINAL_verified`；不要混用两个目录。
- `v3_verify.py` 独立调用冻结 Python 对象，进行整段/逐阶段与事件时间对照；原始 Python 与 Simulink MAT 分开保存。
- `v3_post.py`、`v3_deliver.py` 提供科学绘图、时间对齐误差、报告、原始文件清单和 zip 校验；`rebuild_v3.m` 是完整重建入口。

## 保留的开发问题和环境说明

- 沙箱内 MATLAB 启动曾出现 File system inconsistency；在用户授权环境运行后可用，真实版本 R2025b。
- 早期停止计时与复位放在同一个 MATLAB Function 引入人工代数环提示，现已拆开；原失败运行输出保留在 `results/debug_attempt_01`。
- 初次日志矩阵转置处理出错，现按时间维统一；未将该开发输出算成有效验证。
- 一次批处理使用不存在的参数名 `rise`，在场景初始化时失败（正确参数名 `rise_s`），当时组合仿真尚未开始；随后用明确的 `run_v3_final_validation` 入口运行。
- 本轮第三方 Python 库装在工程 `.python_deps`，打包仅提供精确 requirements，不打入第三方二进制。项目参考模型与配置源码全部包含。
- 各有效 MAT 均包含 `simulation_metadata`、`simulation_output`、`last_warning`/`last_warning_id`（无警告时为空），每例独立日志与显式状态。

`source_changes_v3.csv` 给出代码文件与 v2 的哈希差异；`MANIFEST_v3.json` 给出交付文件最终哈希。旧 `manifest.json` 属于 v2 历史资料，不作 v3 文件校验依据。
