"""Build an evidence-based report, comparison figures, inventory and portable zip."""
from pathlib import Path
import sys,json,csv,hashlib,zipfile,datetime,math,os
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'.python_deps'))
import numpy as np
from scipy.io import loadmat
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2),encoding='utf-8')
def diag(tag,target=300,mode=2):return read(ROOT/'results'/tag/f'case_{target}mL_mode{mode}.json')['diagnostics']
def fmt(v,n=5):return '—' if v is None or not np.isfinite(v) else f'{v:.{n}f}'
def table(headers,rows):return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+['| '+' | '.join(map(str,r))+' |' for r in rows])+'\n'
def csvwrite(name,rows):
    if not rows:return
    with (ROOT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def python_figures():
    for tag in ['D_diagnostic_baseline','FINAL_verified','B_guard_localised']:
        for p in (ROOT/'validation/python_trajectories'/tag).glob('*_aligned_errors.mat'):
            d=loadmat(p,simplify_cells=True);t=np.asarray(d['t']).ravel();err=d['error'];y=d['simulink_y'];same=np.asarray(d['same_branch'],bool)
            closed=y[:,9]>.5;tc=t[closed][0] if np.any(closed) else None
            for scope in ['full','startup','closure']:
                if scope=='closure' and tc is None:continue
                dest=p.with_name(p.stem+'_'+scope+'.png')
                if dest.exists() and dest.stat().st_mtime>=p.stat().st_mtime:continue
                fig,axs=plt.subplots(4,2,figsize=(11,9),sharex=True)
                cols=[1,12,5,25,24,4,23,39];units=['Q (mL/s)','Collected (mL)','Command (STEP/s)','Actual speed (rad/s)','Torque (N m)','Contact force (N)','Phase (rad)','Overlap (m)']
                for ax,col,unit in zip(axs.ravel(),cols,units):
                    for c in [False,True]:
                        mask=same & (closed==c);ax.plot(t[mask],err[mask,col],lw=.7)
                    ax.axhline(0,color='gray',lw=.5);ax.grid(alpha=.2);ax.set_ylabel('Error: '+unit)
                    if scope=='startup':ax.set_xlim(0,3)
                    elif scope=='closure':ax.set_xlim(max(0,tc-2),min(t[-1],tc+2))
                axs[-1,0].set_xlabel('Time (s)');axs[-1,1].set_xlabel('Time (s)')
                fig.suptitle(f'Python minus Simulink | {tag} | {p.stem}\nSame-time comparison; no interpolation across valve reset',fontsize=10)
                fig.tight_layout(rect=(0,0,1,.95));fig.savefig(dest,dpi=160);plt.close(fig)

def report():
    env=read(ROOT/'logs/environment.json');pt=read(ROOT/'validation/point_verification.json');ut=read(ROOT/'validation/unit_tests.json')
    audit=read(ROOT/'validation/frozen_parameter_diff.json');selection=read(ROOT/'validation/selected_controller.json')
    checks=[]
    for tag,modes,targets in [('D_diagnostic_baseline',[1,2,3],[300]),('FINAL_verified',[2],[200,300,400]),('B_guard_localised',[2,3],[300])]:
        for target in targets:
            for mode in modes:
                p=ROOT/'validation/python_trajectories'/tag/f'case_{target}mL_mode{mode}_comparison.json';a=read(p)
                checks.append(dict(case=f'{tag}/{target}/mode{mode}',passed=a['passed'],worst_scaled_error=max(v['maximum_scaled_error'] for v in a['phase_results'].values()),max_event_error_s=max(abs(r['difference_s']) for r in a['event_comparison'] if r['difference_s'] is not None and np.isfinite(r['difference_s']))))
    conv=[]
    for label in ['convergence_300_mode1','convergence_300_mode2','convergence_300_mode3','final_verified_convergence_300']:
        a=read(ROOT/'validation'/label/'comparison.json');conv.append(dict(case=label,passed=a['passed'],worst_scaled_error=max(v['maximum_scaled_error'] for v in a['phases'].values())))
    final=[diag('FINAL_verified',t) for t in [200,300,400]]
    balance_keys=['state_integral_vs_bladder_mL','flow_integral_vs_bladder_mL','tail_integral_vs_drained_mL','total_inventory_identity_error_mL','V_est_vs_state_mL']
    balance_checks=[dict(target_mL=d['target_mL'],absolute_tolerance_mL=.002,
        maximum_error_mL=max(abs(d['balance'][k]) for k in balance_keys),
        passed=all(abs(d['balance'][k])<=.002 for k in balance_keys) and d['balance']['minimum_pipe_inventory_mL']>=0) for d in final]
    allpass=bool(pt['passed'] and ut['failed']==0 and audit['all_unchanged'] and all(x['passed'] for x in checks+conv+balance_checks) and all(d['completed'] and d['meets_plateau_5pct'] and d['meets_volume_2mL'] for d in final))
    write(ROOT/'validation/acceptance_v3.json',dict(passed=allpass,scope='Final conditional full-model numerical validation, NOT hardware validation',trajectory_checks=checks,convergence_checks=conv,balance_checks=balance_checks,known_failed_design='0.25 mL early stop with PI/FF+PI',step10_started=False))
    summary=[];guards=[];events=[]
    for p in sorted((ROOT/'results').glob('*/case*.json')):
        if '.status.' in p.name or p.parent.name.startswith('debug'):continue
        a=read(p)
        if 'diagnostics' not in a:continue
        d=a['diagnostics'];r={'tag':p.parent.name,'target_mL':d['target_mL'],'mode':d['mode'],'status':d['status']}
        for k in ['plateau_rmse_mL_s','plateau_err_pct','collected_mL','volume_error_mL','stop_time_s','max_abs_torque_Nm','torque_margin_Nm','torque_utilisation_pct','max_contact_force_N','electrical_phase_boundary_margin_deg','postclose_peak_force_N','postclose_travel_m','postclose_pressure_increase_Pa','command_start_s','command_stop_s','mechanical_stop_s','t_close_s','tail_end_s','upper_saturation_s','lower_saturation_s','slew_limited_s']:
            r[k]=d.get(k)
        r['plateau_overshoot_pct']=max(0,d['plateau']['overshoot_pct']);r['plateau_band_entry_s']=d['plateau']['settled_in_5pct_band_s'];summary.append(r)
        for g in d['guards']:guards.append(dict(tag=p.parent.name,target_mL=d['target_mL'],mode=d['mode'],**g))
        be=d.get('boundary_events',[]);be=[be] if isinstance(be,dict) else be
        for e in be:events.append(dict(tag=p.parent.name,target_mL=d['target_mL'],mode=d['mode'],**e))
    csvwrite('controller_comparison_v3.csv',summary);csvwrite('boundary_margins_v3.csv',guards);csvwrite('boundary_events_v3.csv',events)
    txt=['# ALBA 第九步 v3 验证与优化报告\n',f'生成：{datetime.datetime.now().isoformat(timespec="seconds")}。实际运行环境：MATLAB {env["version"]}，Simulink 许可证测试={env["simulink_license_test"]}。\n',
    '## 结论\n',f'最终条件模型数值验收：**{"通过" if allpass else "存在未通过项，见 acceptance_v3.json"}**。主模型为 `stage9b_port/alba9b_full_loop.slx`，最终原始结果为 `results/FINAL_verified`。所有曲线来自实际 Simulink 仿真；Python 为独立参考校验，分开保存。\n',
    '300 mL 筛选后选择 **纯流量 PI，3 s 平滑启动，距离关阀阈值 0.05 mL 时停止收线指令**。Kp=50、Ki=600、Kaw=10 不变。冻结后验证 200/400 mL，没有逐目标调参。选择纯 PI 是为降低启动峰值转矩；前馈＋PI 的平台 RMSE 可以更小，并非被预设为一定最好。\n',
    '15 mL/s 定位始终为 **Upper-feasible / near-limit case under nominal model assumptions**。模型内通过不等于硬件可行，末端近零体积误差不代表实物控制精度。本轮未开展第十步、未改 ESP32、未重新辨识 SINDy、未增加 phi*x_dot。\n',
    '## 冻结物理对象与单位\n',
    '63 项数值参数/ROM 数组精确比较全部不变；原参数 MAT、两个物理常量源码与 v2 的 SHA-256 相同，见 `validation/frozen_parameter_diff.json`。独立原始 v2 未改动。未收到单独的 Claude 未完成包。\n',
    '- Rh=1.6440083045256217e9 Pa·s/m³，统一名称 **bench effective resistance**。保留原 Re 相关直管修正和避免重复层流损失处理，未新增阀/流量计损失。\n- h=+0.10 m，助流压头 +981 Pa；裸轴半径 0.0025 m、3200 STEP/rev。\n- 条件性平坦转矩包络 0.29658517619749203 N·m，未提高；12 V、约 0.67 A 仍是工况/估计，不是实测运行转矩曲线。\n- Ecoflex、70/66 mm、270 g 自重、摩擦、阻尼、惯性、传力几何、液压惯性均沿用原模型。旧 15 N 不作本阶段电机上限；118 N 仅为 ROM 刚度审查范围，完整模型未截力。\n- t=0 阀全开并启用控制，初始鱼线张紧、无额外预载；自然被动排液/接触间隙均保留，未设置等待。\n',
    '控制器为 v2 的连续时间 PI（Ts=0，**未模拟离散采样/固件调度**），误差 mL/s，输出 STEP/s。n_ff=Q_ref/(ku*Ax(V_est)) 时在接口内部将 Q 转 SI；ku=2πr/3200，Ax=D²/2。V_est=V0−∫Q_actual dt，只用理想 Q 反馈，不输入真实接触力。总原始指令先限幅，实际受限 n 用于抗积分饱和。指令 0–3000 STEP/s，开阀斜率 6000 STEP/s²、闭阀制动 2000 STEP/s²，tau_slew=0.005 s 是指令整形参数，不是物理电机时间常数。\n',
    '## 程序修正（与优化分开）\n',
    '1. 指令停止、实际机械停止、关阀和尾流结束单独记录。机械停止要求 |omega|≤0.001 rad/s、|acc|≤0.001 m/s² 连续 0.1 s，不能仅用 n≈0；事件时间按原始记录给出采样括区。n=0 不等于断电，实际速度、转矩、接触力未清零。\n2. 完成分类区分 normal_completion、boundary_triggered、timeout、solver_failure、manual_interruption；另保留 implementation_failure/unknown。正常完成需要单次关阀、尾流结束、运动稳定且无边界触发。\n3. 饱和按限幅前后指令判断，上下限及六种完成分类均有单元测试。各守卫分别报告原单位，禁止混合单位排序。\n4. 相位指标准确称为“距电相位模型边界的余量”，不是经典稳定性相位裕度。\n5. 保留 v2 关阀锁存和单次 Q 复位；增加守卫零交叉定位，不更改边界数值。\n6. 200 mL 初次最终验证仅仿真结束时间相差 2.924 ms，超过预设 2 ms；保留 `FINAL_selected` 及对应失败 JSON。补上尾流结束面零交叉后重跑 `FINAL_verified`，没有放宽判据或改变物理轨迹/控制参数。\n',
    '7. 最后连线审计将 Ref 与 ValveCmp 的体积输入也从 v2 的 V0−V 明确改为独立 ∫Q。理想模型下两者数值等价，但这样控制侧不再读取真实体积；原实现留在 debug_pre_wiring_audit，最终三组与严格容差结果均复跑。connection_audit.json 同时自动断言这两条连线来源。\n',
    '开发期出现过人工代数环提示、日志矩阵方向错误和一次调用参数名 rise（正确名 rise_s）的入口错误；它们不是有效结果，已修复，早期失败输出保留在 debug 目录/运行说明。最终模型编译与完整运行以当前日志和原始 MAT 为准。\n',
    '## 数值验证\n',f'逐点测试 {pt["points"]} 点通过，单元测试 {ut["passed"]} 项通过。包括负重叠、重新接触邻域、约 116 N/120 N、78.6°/90° 邻域、102 mL/Vmin 邻域及各边界两侧。越界点仅用于方程移植比较，仿真没有越界继续。逐变量 abs_tol+rel_tol*|reference| 见 `point_test_spec.json`、`point_verification.json`。\n',
    '全轨迹比较保留启动、平台、减速、关阀后四段。在同一时间对齐独立 Python 解与 Simulink 原始日志，排除仅由事件时间微差造成的不同阀状态样本，不跨 Q 复位平滑插值。除以下表格外，每例 JSON 给出逐变量绝对误差、相对判据和各事件差。误差/容差比≤1 为通过。\n',
    table(['案例','通过','最大误差/容差','最大事件差 (s)'],[[r['case'],r['passed'],fmt(r['worst_scaled_error'],4),fmt(r['max_event_error_s'],7)] for r in checks]),
    '完整对象使用 ode15s；名义 RelTol=2e-8、MaxStep=0.025 s；收敛复跑 RelTol=2e-10、MaxStep=0.005 s、各状态 AbsTol 再缩小 100 倍。绝对容差按状态单位分别指定，见每例 solver 字段。\n',
    table(['收敛复跑','通过','最大误差/容差'],[[r['case'],r['passed'],fmt(r['worst_scaled_error'],4)] for r in conv]),
    'ROM 没有额外串联执行器，也没有重训。300 mL 启动段同时间最大误差：前馈＋PI 为 Q≈2.2024 mL/s、F≈20.0135 N；纯 PI 为 Q≈2.3498 mL/s、F≈21.2353 N。两个独立峰值接近不能说明整条曲线准确。完整模型在被动排水后出现负重叠和重新接触；消去快电机动态的 ROM 不保留同等间隙记忆，这是启动偏差的重要限制。见 `validation/ROM_full_*` 时间对齐误差 MAT/图。ROM 不单独证明电机不失步。\n',
    '## 300 mL 公平基线\n',
    table(['控制器（2 s）','平台 RMSE mL/s','最差偏差 %','平台超调 %','进入并保持 ±5% 时刻 s','峰值转矩 N·m','包络占比 %','完成 s'],[[diag('D_diagnostic_baseline',mode=m)['mode_name'],fmt(diag('D_diagnostic_baseline',mode=m)['plateau_rmse_mL_s']),fmt(diag('D_diagnostic_baseline',mode=m)['plateau_err_pct'],3),fmt(diag('D_diagnostic_baseline',mode=m)['plateau']['overshoot_pct'],3),fmt(diag('D_diagnostic_baseline',mode=m)['plateau']['settled_in_5pct_band_s'],4),fmt(diag('D_diagnostic_baseline',mode=m)['max_abs_torque_Nm'],8),fmt(diag('D_diagnostic_baseline',mode=m)['torque_utilisation_pct'],3),fmt(diag('D_diagnostic_baseline',mode=m)['stop_time_s'],4)] for m in [1,2,3]]),
    '固定步频保持原始 606.94 STEP/s 和受限启动、统一阈值/尾流；平台 ±5% 失败，如实保留。纯 PI 和前馈＋PI 均满足基线平台指标。表中负的峰值偏差表示全平台未超过目标，其正向超调为 0（CSV 已按非负超调输出）。进入带时刻是在平台窗口内确认并保持的时刻，不是启动过程中第一次瞬时穿越目标带。\n',
    '## 分开的优化试验\n',
    table(['启动时长 s','控制器','峰值转矩 N·m','占比 %','最大径向间隙 µm','平台 RMSE mL/s'],[[rise,diag(tag,mode=m)['mode_name'],fmt(diag(tag,mode=m)['max_abs_torque_Nm'],8),fmt(diag(tag,mode=m)['torque_utilisation_pct'],3),fmt(diag(tag,mode=m)['maximum_radial_gap_um'],3),fmt(diag(tag,mode=m)['plateau_rmse_mL_s'])] for rise,tag in [(2,'D_diagnostic_baseline'),(2.5,'A_startup_2p5'),(3,'A_startup_3')] for m in [2,3]]),
    '更慢启动降低了这组名义模型的转矩峰，但增大被动排水形成的间隙；未以消除约 0.5 s 延迟为优化目标。终止优化仍保持 2 s 启动，先与启动优化分开：\n',
    table(['终止方案','控制器','状态','已收集 mL','关阀后峰值接触力 N','关阀后继续行程 µm'],[[tag,diag(tag,mode=m)['mode_name'],diag(tag,mode=m)['status'],fmt(diag(tag,mode=m)['collected_mL'],6),fmt(diag(tag,mode=m)['postclose_peak_force_N'],4),fmt(None if diag(tag,mode=m)['postclose_travel_m'] is None else diag(tag,mode=m)['postclose_travel_m']*1e6,3)] for tag,modes in [('D_diagnostic_baseline',[1,2,3]),('B_terminal_0p05',[1,2,3]),('B_guard_localised',[2,3])] for m in modes]),
    '固定步频的管理变体按剩余体积降低原固定步频，不加入流量 PI；它仍不能满足平台指标。0.25 mL 提前停止的两种 PI 在约 297.908 mL、Q_actual≈2.10 mL/s 时触发张力边界，阀未关，F 已降至 0。失败并未被补写成 300 mL，也未放宽张力边界。该事件不是 118 N 力上限触发。完整事件值见 `boundary_events_v3.csv`。0.05 mL 变体降低关阀后的继续压紧，随后与 3 s 启动组合。\n',
    '## 冻结最终方案：200 / 300 / 400 mL\n',
    table(['目标 mL','状态','平台 RMSE mL/s','最差偏差 %','体积误差 mL','完成 s','转矩余量 N·m'],[[d['target_mL'],d['status'],fmt(d['plateau_rmse_mL_s']),fmt(d['plateau_err_pct'],3),fmt(d['volume_error_mL'],8),fmt(d['stop_time_s'],5),fmt(d['torque_margin_Nm'],8)] for d in final]),
    table(['目标 mL','指令开始 s','关阀 s','机械停止 s','指令停止 s','尾流结束 s','最大步频 STEP/s','最大实际 rpm'],[[d['target_mL'],fmt(d['command_start_s'],5),fmt(d['t_close_s'],5),fmt(d['mechanical_stop_s'],5),fmt(d['command_stop_s'],5),fmt(d['tail_end_s'],5),fmt(d['max_n_steps_s'],2),fmt(d['max_rotor_rpm'],3)] for d in final]),
    '机械停止可能早于指令 n 完全趋近零，这是电机实际运动与 STEP 指令分开建模的结果。停止事件含数值阈值，不代表瞬间严格静止或断电。\n',
    table(['目标 mL','关阀后峰值力 N','关阀后继续行程 µm','关阀后压力增量 Pa','剩余膀胱量 mL','尾流时长 s'],[[d['target_mL'],fmt(d['postclose_peak_force_N'],4),fmt(d['postclose_travel_m']*1e6,3),fmt(d['postclose_pressure_increase_Pa'],2),fmt(d['min_bladder_mL'],6),fmt(d['tail_duration_s'],5)] for d in final]),
    '400 mL 的关阀后继续压紧仍明显（约 24.4 N、约 4.9 kPa 压力增量），不能因为排量误差很小就说终止问题已完全消失。下一轮硬件/测量验证前需保留这一限制；本轮不开展参数不确定性扫描。\n',
    '### 各守卫与饱和（分别有单位）\n',
    table(['目标 mL','守卫','单位','最小裕度','触发'],[[d['target_mL'],g['name'],g['unit'],f'{g["minimum_margin"]:.9g}',g['violated']] for d in final for g in d['guards']]),
    table(['目标 mL','上限饱和 s','下限饱和 s','斜率限幅 s','最大开阀加速度 STEP/s²','最大闭阀制动 STEP/s²','距电相位模型边界 °'],[[d['target_mL'],fmt(d['upper_saturation_s'],5),fmt(d['lower_saturation_s'],5),fmt(d['slew_limited_s'],5),fmt(d['max_accel_open_steps_s2'],3),fmt(d['max_accel_closed_steps_s2'],3),fmt(d['electrical_phase_boundary_margin_deg'],4)] for d in final]),
    '### 守恒与尾水\n',
    table(['目标 mL','独立 ∫Q 状态−膀胱排量 mL','日志 Q 积分误差 mL','尾流积分误差 mL','总库存恒等误差 mL','最小管路库存 mL'],[[d['target_mL'],f'{d["balance"]["state_integral_vs_bladder_mL"]:.3g}',f'{d["balance"]["flow_integral_vs_bladder_mL"]:.3g}',f'{d["balance"]["tail_integral_vs_drained_mL"]:.3g}',f'{d["balance"]["total_inventory_identity_error_mL"]:.3g}',fmt(d['balance']['minimum_pipe_inventory_mL'],6)] for d in final]),
    '总库存恒等式与独立日志积分检查分开报告，前者本身不是独立数值证明。阀关后 Q_valve=0，尾流从管路库存扣除，不再次扣膀胱。名义 198/298/398 mL 过阀后关阀；保留 v2 锁存 5 µs 延迟，解释约 1.5e-5 mL 的名义超量。\n',
    '参考减速按剩余体积开始（R0=28.5 mL，到 r_hold=1.5 mL 达 3 mL/s），早于关阀阈值；并非把实际流量设为参考。尾流约 2 mL，半余弦 T=2Vtail/q0 是假设；q0 为关阀前 1 ms 跟踪值而非精确左极限。最终 PI 尾流时长约 1.35 s，与“约 1.33 s”接近，但不是实验通用常数。原固定步频关阀流量较高，对应 T≈0.46 s；其管理变体流量更低、T≈2.4 s，体现假设随工况变化，不能解释成实测保证。\n',
    '## 交付与下一步\n',
    '入口和依赖说明见 `README_v3.md`；可一键 `rebuild_v3(pythonExe)` 在新副本重建并实际运行。原始 Simulink MAT 保存所有日志/状态/SimulationOutput，CSV 只是抽样视图。完整对照 CSV、各阶段图、Python 原始 MAT/误差 MAT、版本/运行日志和参数差异均在包内。poster 图在 `poster`，保留条件模型注记。\n',
    '当前证据支持同一冻结控制器在这三个名义工况内通过数值与设计评估，仍受平坦转矩假设、理想流量反馈、连续控制器、阀/尾水近似和完整模型有效范围限制。ROM 启动偏差、0.25 mL 提前停止失败、400 mL 闭阀后压紧均未隐去。建议先审阅这些结果并决定后续测量/硬件测试安排；本轮在第九步交付处停止。\n']
    (ROOT/'REPORT_v3.md').write_text('\n'.join(txt),encoding='utf-8')
    print('V3 ACCEPTANCE',allpass,'cases',len(summary),'trajectory checks',len(checks),flush=True)
    python_figures()
    return allpass

def included():
    excluded={'.python_deps','slprj','__pycache__','.git'}
    files=[]
    for base,dirs,names in os.walk(ROOT):
        dirs[:]=[d for d in dirs if d not in excluded]
        files.extend(Path(base)/n for n in names if Path(n).suffix not in ['.slxc','.pyc'] and n!='MANIFEST_v3.json')
    return sorted(files)

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()

def package():
    files=included();manifest=[]
    for p in files:manifest.append(dict(path=p.relative_to(ROOT).as_posix(),bytes=p.stat().st_size,sha256=sha(p)))
    write(ROOT/'MANIFEST_v3.json',dict(created=datetime.datetime.now().isoformat(),files=manifest,excluded=['Third-party .python_deps','MATLAB slprj/slxc cache','Python bytecode'],primary_model='stage9b_port/alba9b_full_loop.slx',primary_report='REPORT_v3.md'))
    zipfile_path=ROOT.parent/'ALBA_Step9_ClosedLoop_v3.zip'
    with zipfile.ZipFile(zipfile_path,'w',zipfile.ZIP_DEFLATED,compresslevel=4,allowZip64=True) as z:
        for p in files+[ROOT/'MANIFEST_v3.json']:z.write(p,Path(ROOT.name)/p.relative_to(ROOT))
    with zipfile.ZipFile(zipfile_path) as z:
        bad=z.testzip();assert bad is None,bad;count=len(z.namelist())
    digest=sha(zipfile_path)
    write(ROOT.parent/'ALBA_Step9_ClosedLoop_v3_package.json',dict(zip=str(zipfile_path),bytes=zipfile_path.stat().st_size,sha256=digest,entries=count,crc_verified=True))
    print('PACKAGE',str(zipfile_path),zipfile_path.stat().st_size,digest,flush=True)

if __name__=='__main__':
    if sys.argv[1]=='report':sys.exit(0 if report() else 1)
    elif sys.argv[1]=='package':package()
