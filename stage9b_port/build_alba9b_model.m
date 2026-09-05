function mdl = build_alba9b_model(mdl)
%BUILD_ALBA9B_MODEL Stage-9 B: the SAME controller and termination logic on the
%FULL five-state reference plant, ported from the frozen Python source and
%verified pointwise against it (see alba9b_verify_port / port_verification.json).
%
%   WHAT IS DIFFERENT FROM MODEL A
%     The plant is no longer ROM3. The motor acts only through the torque-phase
%     relation, so the ACTUAL cable speed and the flow are both solved:
%         phi    = (Nfull/4)*(2*pi*count/steps_per_rev - x/r)
%         torque = torque_cap*sin(phi)
%         M*a    = torque/r - damping*v - friction*tanh(v/vs) - Fn
%         Ih*dQ  = pe + pg + Fn/A + rho*g*h - loss(Q)
%     count is the integral of the STEP/s command, so the command enters as a
%     commanded POSITION, exactly as the frozen model has it. Nothing sets the
%     actual speed equal to the commanded speed and nothing prescribes Q.
%
%   WHAT IS IDENTICAL TO MODEL A
%     The reference generator, the controller (gains, structure, anti-windup,
%     slew and saturation), the monotone valve latch, the single Q reset, the
%     closing-instant capture and the tail. No retuning in this stage.
%
%   GUARDS ARE REPORTED, NOT ENFORCED. If the phase leaves +-pi/2, the cable
%   tension goes negative, the travel or the sphere limit is reached, the run
%   STOPS there and the failure is kept. Nothing is clipped to make a curve look
%   acceptable.
if nargin < 1, mdl = 'alba9b_full_loop'; end

here = fileparts(mfilename('fullpath'));
slx = fullfile(here,[mdl '.slx']);
if bdIsLoaded(mdl), close_system(mdl,0); end
if isfile(slx), delete(slx); end
load_system('simulink');
new_system(mdl);
add = @(src,name,pos,varargin) add_block(src,[mdl '/' name],'Position',pos,varargin{:});

% ---------------- sources ----------------
add('simulink/Sources/Clock','Clock',[30 470 60 500]);
add('simulink/Sources/Constant','S',[30 540 90 570],'Value','S_vec');
add('simulink/Sources/Constant','V0',[30 300 90 330],'Value','V0_m3');
add('simulink/Sources/Constant','CloseThresh',[30 380 130 410],'Value','S_vec(1)-S_vec(2)');
add('simulink/Sources/Constant','Half',[1040 410 1075 440],'Value','0.5');

% ---------------- plant states (SI, exactly as the frozen model) ----------
add('simulink/Continuous/Integrator','Int_count',[760 20 800 60],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-6');
add('simulink/Continuous/Integrator','Int_x',[760 80 800 120],'InitialCondition','0', ...
    'AbsoluteTolerance','2e-11');
add('simulink/Continuous/Integrator','Int_v',[760 140 800 180],'InitialCondition','0', ...
    'AbsoluteTolerance','2e-10');
add('simulink/Continuous/Integrator','Int_V',[760 200 800 240],'InitialCondition','V0_m3', ...
    'AbsoluteTolerance','2e-13');
add('simulink/Continuous/Integrator','Int_Q',[760 260 800 305],'InitialCondition','0', ...
    'AbsoluteTolerance','2e-14','ExternalReset','rising');
add('simulink/Continuous/Integrator','Int_intQ',[760 320 800 360],'InitialCondition','0', ...
    'AbsoluteTolerance','2e-13');

% ---------------- controller states ----------------
add('simulink/Continuous/Integrator','Int_n',[470 60 510 100],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-4','LimitOutput','on', ...
    'UpperSaturationLimit','S_vec(12)','LowerSaturationLimit','0');
add('simulink/Continuous/Integrator','Int_ni',[470 130 510 170],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-4');

% ---------------- valve latch and tail states ----------------
add('simulink/Continuous/Integrator','Int_latch',[1090 328 1125 362],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-9','LimitOutput','on', ...
    'UpperSaturationLimit','1','LowerSaturationLimit','0');
add('simulink/Continuous/Integrator','Int_tc',[760 520 800 560],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-7');
add('simulink/Continuous/Integrator','Int_q0',[760 580 800 620],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-7');
add('simulink/Continuous/Integrator','Int_tail',[760 700 800 740],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-9');

% ---------------- unit interfaces ----------------
add('simulink/Math Operations/Gain','V_mL',[850 200 890 240],'Gain','1e6');
add('simulink/Math Operations/Gain','Q_mL',[850 265 890 300],'Gain','1e6');
add('simulink/Math Operations/Gain','intQ_mL',[850 320 890 360],'Gain','1e6');
add('simulink/Math Operations/Sum','DischargedSI',[930 290 950 310],'Inputs','+-');
add('simulink/Math Operations/Gain','Discharged',[975 285 1015 315],'Gain','1e6');
add('simulink/Math Operations/Sum','Collected',[880 760 900 780],'Inputs','++');

% ---------------- valve: raw event -> monotone latch -> flag ----------------
add('simulink/Logic and Bit Operations/Relational Operator','ValveCmp',[930 330 970 360], ...
    'Operator','>=','ZeroCross','on');
add('simulink/Signal Attributes/Data Type Conversion','CmpD',[985 330 1025 360], ...
    'OutDataTypeStr','double');
add('simulink/Math Operations/Gain','LatchGain',[1040 335 1075 355],'Gain','S_vec(18)');
add('simulink/Logic and Bit Operations/Relational Operator','ValveCmp2',[1150 330 1190 360], ...
    'Operator','>=','ZeroCross','on');
add('simulink/Signal Attributes/Data Type Conversion','ValveD',[1210 330 1250 360], ...
    'OutDataTypeStr','double');

% ---------------- functions ----------------
add('simulink/User-Defined Functions/MATLAB Function','Ref',  [250 350 370 430]);
add('simulink/User-Defined Functions/MATLAB Function','Ctrl', [250 60 400 220]);
add('simulink/User-Defined Functions/MATLAB Function','Plant',[560 20 700 360]);
add('simulink/User-Defined Functions/MATLAB Function','Latch',[560 520 690 620]);
add('simulink/User-Defined Functions/MATLAB Function','Tail', [560 660 690 790]);

setChart(mdl,'Ref',[ ...
 'function Qref = fcn(t,discharged,valve_closed,S)\n' ...
 '%%#codegen\n' ...
 'Qref = alba9_reference(t,discharged,valve_closed,S);\n']);

setChart(mdl,'Ctrl',[ ...
 'function [n_in,dni,n_raw,e_track,sat_flag,slew_flag] = fcn(Qref,Q,V,n,ni,valve_closed,S)\n' ...
 '%%#codegen\n' ...
 '%% identical to model A: same gains, same structure, no retuning in stage B\n' ...
 '[n_in,dni,n_raw,e_track,~,sat_flag,slew_flag] = ...\n' ...
 '    alba9_controller(Qref,Q,V,n,ni,valve_closed,S);\n']);

setChart(mdl,'Plant',[ ...
 'function [dx,dv,dV,dQ,Fn,phi,torque,omega,tension,acc,gmin,gidx,stop] = ...\n' ...
 '        fcn(x,v,V,Q,count,valve_closed)\n' ...
 '%%#codegen\n' ...
 '[dx,dv,dV,dQ,Fn,phi,torque,omega,tension,acc,gmin,gidx,stop] = ...\n' ...
 '    alba9b_plant_block(x,v,V,Q,count,valve_closed);\n']);

setChart(mdl,'Latch',[ ...
 'function [dtc,dq0] = fcn(valve_closed,Q,q0,S)\n' ...
 '%%#codegen\n' ...
 'isopen = 1 - min(max(valve_closed,0),1);\n' ...
 'dtc = isopen;\n' ...
 'dq0 = isopen*(Q-q0)/S(19);\n']);

setChart(mdl,'Tail',[ ...
 'function [q_tail,done,q_outlet] = fcn(t,valve_closed,q0,tclose,drained,n,Q,S)\n' ...
 '%%#codegen\n' ...
 '[q_tail,done,q_outlet] = alba9_tail(t,valve_closed,q0,tclose,drained,n,Q,S);\n']);

% ---------------- stop conditions ----------------
add('simulink/Sinks/Stop Simulation','StopDone',[760 820 790 850]);
add('simulink/Sinks/Stop Simulation','StopGuard',[760 870 790 900]);

% ---------------- logging ----------------
add('simulink/Signal Routing/Mux','Mux',[1320 20 1325 800],'Inputs','32');
add('simulink/Sinks/To Workspace','ToWs',[1380 400 1450 430], ...
    'VariableName','alba9b_log','SaveFormat','Structure With Time','SampleTime','-1');

% ---------------- wiring ----------------
L = @(a,b) add_line(mdl,a,b,'autorouting','on');

% plant inputs and state derivatives
L('Int_x/1','Plant/1');  L('Int_v/1','Plant/2');  L('Int_V/1','Plant/3');
L('Int_Q/1','Plant/4');  L('Int_count/1','Plant/5'); L('ValveD/1','Plant/6');
L('Plant/1','Int_x/1');  L('Plant/2','Int_v/1');  L('Plant/3','Int_V/1');
L('Plant/4','Int_Q/1');  L('ValveD/1','Int_Q/2');
L('Int_n/1','Int_count/1');           % commanded position = integral of STEP/s
L('Int_Q/1','Int_intQ/1');            % independent flow integral
L('Plant/13','StopGuard/1');          % stop on any domain violation, do not clip

% unit interfaces
L('Int_V/1','V_mL/1');  L('Int_Q/1','Q_mL/1');  L('Int_intQ/1','intQ_mL/1');
L('V0/1','DischargedSI/1');  L('Int_V/1','DischargedSI/2');
L('DischargedSI/1','Discharged/1');

% valve latch
L('Discharged/1','ValveCmp/1');  L('CloseThresh/1','ValveCmp/2');
L('ValveCmp/1','CmpD/1');        L('CmpD/1','LatchGain/1');
L('LatchGain/1','Int_latch/1');
L('Int_latch/1','ValveCmp2/1');  L('Half/1','ValveCmp2/2');
L('ValveCmp2/1','ValveD/1');

% reference and controller
L('Clock/1','Ref/1');  L('Discharged/1','Ref/2');  L('ValveD/1','Ref/3');  L('S/1','Ref/4');
L('Ref/1','Ctrl/1');   L('Q_mL/1','Ctrl/2');   L('V_mL/1','Ctrl/3');
L('Int_n/1','Ctrl/4'); L('Int_ni/1','Ctrl/5'); L('ValveD/1','Ctrl/6'); L('S/1','Ctrl/7');
L('Ctrl/1','Int_n/1'); L('Ctrl/2','Int_ni/1');

% closing-instant latches and tail
L('ValveD/1','Latch/1');  L('Q_mL/1','Latch/2');  L('Int_q0/1','Latch/3');  L('S/1','Latch/4');
L('Latch/1','Int_tc/1');  L('Latch/2','Int_q0/1');
L('Clock/1','Tail/1');    L('ValveD/1','Tail/2');    L('Int_q0/1','Tail/3');
L('Int_tc/1','Tail/4');   L('Int_tail/1','Tail/5');  L('Int_n/1','Tail/6');
L('Q_mL/1','Tail/7');     L('S/1','Tail/8');
L('Tail/1','Int_tail/1'); L('Tail/2','StopDone/1');
L('Discharged/1','Collected/1');  L('Int_tail/1','Collected/2');

% logging, order per alba9b_log_index
L('Ref/1','Mux/1');         L('Q_mL/1','Mux/2');      L('V_mL/1','Mux/3');
L('Discharged/1','Mux/4');  L('Plant/5','Mux/5');     L('Int_n/1','Mux/6');
L('Ctrl/3','Mux/7');        L('Int_x/1','Mux/8');     L('Int_v/1','Mux/9');
L('ValveD/1','Mux/10');     L('Tail/1','Mux/11');     L('Int_tail/1','Mux/12');
L('Collected/1','Mux/13');  L('Ctrl/4','Mux/14');     L('Ctrl/1','Mux/15');
L('Int_ni/1','Mux/16');     L('Int_q0/1','Mux/17');   L('Int_tc/1','Mux/18');
L('Tail/3','Mux/19');       L('Int_latch/1','Mux/20');L('CmpD/1','Mux/21');
L('Ctrl/5','Mux/22');       L('Ctrl/6','Mux/23');     L('Plant/6','Mux/24');
L('Plant/7','Mux/25');      L('Plant/8','Mux/26');    L('Plant/9','Mux/27');
L('Plant/11','Mux/28');     L('Plant/12','Mux/29');   L('intQ_mL/1','Mux/30');
L('Int_count/1','Mux/31');  L('Plant/10','Mux/32');
L('Mux/1','ToWs/1');

% ---------------- solver ----------------
set_param(mdl, ...
    'SolverType','Variable-step','Solver','ode15s', ...
    'RelTol','2e-8','AbsTol','1e-9','MaxStep','0.025', ...
    'ZeroCrossControl','UseLocalSettings','StopTime','120', ...
    'SaveOutput','off','SaveTime','off','SignalLogging','off');

alba9v3_extend_model(mdl);
save_system(mdl,slx);
fprintf('Built %s\n',slx);
end

% -------------------------------------------------------------------------
function setChart(mdl,name,code)
chart = sfroot().find('-isa','Stateflow.EMChart','Path',[mdl '/' name]);
chart.Script = sprintf(code);
end
