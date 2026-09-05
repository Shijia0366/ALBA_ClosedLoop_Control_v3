# rom3_rhs.m 整数类型修复（2026-09-02）

`ALBA_Step8_SINDy_v1.zip` 里原版 `rom3_rhs.m` 在 MATLAB R2025a 下无法运行。
原因是 `make_results.py` 用 `savemat` 导出时，`stiffness_powers` 和 `steps_per_rev`
落成了 int64，而 MATLAB 对整数类型的运算规则与 NumPy 不同。

两处改动（只改数据类型，模型系数、结构、物理参数全部未动）：

| 行 | 原版 | 修复后 | 后果 |
|---|---|---|---|
| 11 | `(F/100).^powers(:)` | `(F/100).^double(powers(:))` | 原版直接报错：double 底数不能用 int64 指数 |
| 13 | `2*pi*r/steps_per_rev` | `2*pi*r/double(steps_per_rev)` | 原版**不报错但静默出错**：double/int64 返回 int64，ku=4.9e-6 会被四舍五入成 0，电机输入被整条抹掉 |

第二处更危险，因为它不会报错。

## 修复后的验证结果（MATLAB R2025a，ode15s）

五条保留测试轨迹全部跑完，与 Python（Radau）预测的差异：

- 最大流量差 2.4e-6 mL/s
- 最大排量差 2.4e-6 mL
- 对参考模型的流量 RMSE：0.000911 / 0.000770 / 0.001094 / 0.001112 / 0.000789 mL/s
  —— 与 `results/test_metrics.json` 里 Python 的数值逐位一致

这是同一冻结模型的跨语言数值复现，不是新的实验验证。

## 根因建议

若要重新导出 MAT，应在 `make_results.py` 里把
`stiffness_powers=np.array(model.powers)` 改成 `np.array(model.powers,dtype=float)`，
`steps_per_rev` 同理转 float，从源头避免整数类型跨语言传播。
