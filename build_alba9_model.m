function mdl = build_alba9_model(mdl)
%BUILD_ALBA9_MODEL Stage-9 v2 closed-loop Simulink model.
%   The SCRIPT is the source of truth; the .slx is generated and may be deleted
%   and rebuilt at any time. Physics, controller, reference and tail live in
%   shared .m files, so this model and the ode15s cross-check in
%   alba9_check_ode execute literally the same equations.
%
%   STATES
%     f        contact force / 100 N. Int_f is a SATURATED integrator with lower
%              limit 0: that IS the unilateral-contact condition. Clamping dF/dt
%              inside the equation instead gives ode15s a step discontinuity in f
%              whose numerical Jacobian stalls the state at zero, so the chart
%              passes clamp_contact = false.
%     V, Q     bladder water (mL) and valve flow (mL/s).
%     n, ni    command state and PI integral. Int_n is saturated to the
%              conditional STEP/s ceiling and fed a slew-limited input, so the
%              rate limit and the saturation act on the one signal the
%              anti-windup measures against.
%     latch    v2. Monotone valve latch, see below.
%     tc, q0   closing time and pre-closure flow, frozen at closure.
%     tail     tail volume delivered.
%
%   VALVE, v2. The raw comparator alone was not safe. After closure the bladder
%   volume freezes exactly ON the threshold and the solver's zero-crossing search
%   re-evaluates on the open side, so the raw signal chatters: the v1 logs show
%   12 closed->open transitions at 200 mL and 21 at 300 mL pure PI. Each one
%   re-fired the Q reset and restarted the closing-time and q0 latches.
%   v2 feeds the raw comparator through a MONOTONE latch integrator (input >= 0,
%   output clamped to [0 1]) and takes the valve flag as latch >= 0.5. A monotone
%   state crosses 0.5 exactly once and never returns, so one open->closed and zero
%   closed->open transitions are guaranteed BY CONSTRUCTION, not by tolerance
%   tuning, log filtering or smoothing. The raw comparator is still logged
%   (channel cmp_raw) so the chatter remains visible rather than hidden.
%   The latch gain sets a deterministic closing lag of 0.5/latch_gain seconds,
%   reported in the diagnostics.
if nargin < 1
    mdl = 'alba9_closed_loop';
end

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
add('simulink/Sources/Constant','V0',[30 300 90 330],'Value','V0mL');
add('simulink/Sources/Constant','CloseThresh',[30 380 130 410],'Value','S_vec(1)-S_vec(2)');
add('simulink/Sources/Constant','Half',[1040 410 1075 440],'Value','0.5');

% ---------------- states ----------------
add('simulink/Continuous/Integrator','Int_f',[760 60 800 100],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-9','LimitOutput','on', ...
    'UpperSaturationLimit','inf','LowerSaturationLimit','0');
add('simulink/Continuous/Integrator','Int_V',[760 130 800 170],'InitialCondition','V0mL', ...
    'AbsoluteTolerance','1e-7');
add('simulink/Continuous/Integrator','Int_Q',[760 200 800 250],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-7','ExternalReset','rising');
add('simulink/Continuous/Integrator','Int_n',[470 60 510 100],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-4','LimitOutput','on', ...
    'UpperSaturationLimit','S_vec(12)','LowerSaturationLimit','0');
add('simulink/Continuous/Integrator','Int_ni',[470 130 510 170],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-4');
add('simulink/Continuous/Integrator','Int_latch',[1090 328 1125 362], ...
    'InitialCondition','0','AbsoluteTolerance','1e-9','LimitOutput','on', ...
    'UpperSaturationLimit','1','LowerSaturationLimit','0');
add('simulink/Continuous/Integrator','Int_tc',[760 500 800 540],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-7');
add('simulink/Continuous/Integrator','Int_q0',[760 560 800 600],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-7');
add('simulink/Continuous/Integrator','Int_tail',[760 680 800 720],'InitialCondition','0', ...
    'AbsoluteTolerance','1e-9');

% ---------------- algebra ----------------
add('simulink/Math Operations/Sum','Discharged',[880 290 900 310],'Inputs','+-');
add('simulink/Math Operations/Sum','Collected',[880 740 900 760],'Inputs','++');

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
add('simulink/User-Defined Functions/MATLAB Function','Plant',[560 60 690 250]);
add('simulink/User-Defined Functions/MATLAB Function','Latch',[560 500 690 610]);
add('simulink/User-Defined Functions/MATLAB Function','Tail', [560 640 690 780]);

setChart(mdl,'Ref',[ ...
 'function Qref = fcn(t,discharged,valve_closed,S)\n' ...
 '%%#codegen\n' ...
 'Qref = alba9_reference(t,discharged,valve_closed,S);\n']);

setChart(mdl,'Ctrl',[ ...
 'function [n_in,dni,n_raw,e_track,sat_flag,slew_flag] = fcn(Qref,Q,V,n,ni,valve_closed,S)\n' ...
 '%%#codegen\n' ...
 '[n_in,dni,n_raw,e_track,~,sat_flag,slew_flag] = ...\n' ...
 '    alba9_controller(Qref,Q,V,n,ni,valve_closed,S);\n']);

setChart(mdl,'Plant',[ ...
 'function [df,dV,dQ,Ax,keff,F] = fcn(f,V,Q,n,valve_closed)\n' ...
 '%%#codegen\n' ...
 'P = alba9_const();\n' ...
 '%% false: Int_f is a saturated integrator, it enforces f >= 0 itself\n' ...
 '[df,dV,dQ,Ax,keff,F] = alba9_plant(f,V,Q,n,valve_closed,P,false);\n']);

setChart(mdl,'Latch',[ ...
 'function [dtc,dq0] = fcn(valve_closed,Q,q0,S)\n' ...
 '%%#codegen\n' ...
 '%% Both integrators simply stop when the valve latches closed: Int_tc then\n' ...
 '%% holds the closing time and Int_q0 holds the flow it was tracking. Because\n' ...
 '%% the valve flag is monotone, they freeze exactly once and never restart.\n' ...
 '%% Int_q0 is a first-order tracker with time constant S(19), so the captured\n' ...
 '%% q0 is a LAGGED pre-closure flow, not the exact left limit of Q.\n' ...
 'isopen = 1 - min(max(valve_closed,0),1);\n' ...
 'dtc = isopen;\n' ...
 'dq0 = isopen*(Q-q0)/S(19);\n']);

setChart(mdl,'Tail',[ ...
 'function [q_tail,done,q_outlet] = fcn(t,valve_closed,q0,tclose,drained,n,Q,S)\n' ...
 '%%#codegen\n' ...
 '[q_tail,done,q_outlet] = alba9_tail(t,valve_closed,q0,tclose,drained,n,Q,S);\n']);

% ---------------- stop and logging ----------------
add('simulink/Sinks/Stop Simulation','Stop',[760 800 790 830]);
add('simulink/Signal Routing/Mux','Mux',[1320 60 1325 780],'Inputs','23');
add('simulink/Sinks/To Workspace','ToWs',[1380 400 1450 430], ...
    'VariableName','alba9_log','SaveFormat','Structure With Time','SampleTime','-1');

% ---------------- wiring ----------------
L = @(a,b) add_line(mdl,a,b,'autorouting','on');

% plant
L('Int_f/1','Plant/1');  L('Int_V/1','Plant/2');  L('Int_Q/1','Plant/3');
L('Int_n/1','Plant/4');  L('ValveD/1','Plant/5');
L('Plant/1','Int_f/1');  L('Plant/2','Int_V/1');  L('Plant/3','Int_Q/1');
L('ValveD/1','Int_Q/2');           % single reset: the latched flag rises once

% discharged volume, raw event, monotone latch, valve flag
L('V0/1','Discharged/1');  L('Int_V/1','Discharged/2');
L('Discharged/1','ValveCmp/1');  L('CloseThresh/1','ValveCmp/2');
L('ValveCmp/1','CmpD/1');        L('CmpD/1','LatchGain/1');
L('LatchGain/1','Int_latch/1');
L('Int_latch/1','ValveCmp2/1');  L('Half/1','ValveCmp2/2');
L('ValveCmp2/1','ValveD/1');

% reference
L('Clock/1','Ref/1');  L('Discharged/1','Ref/2');  L('ValveD/1','Ref/3');  L('S/1','Ref/4');

% controller
L('Ref/1','Ctrl/1');   L('Int_Q/1','Ctrl/2');   L('Int_V/1','Ctrl/3');
L('Int_n/1','Ctrl/4'); L('Int_ni/1','Ctrl/5');  L('ValveD/1','Ctrl/6');  L('S/1','Ctrl/7');
L('Ctrl/1','Int_n/1'); L('Ctrl/2','Int_ni/1');

% closing-instant latches
L('ValveD/1','Latch/1');  L('Int_Q/1','Latch/2');  L('Int_q0/1','Latch/3');  L('S/1','Latch/4');
L('Latch/1','Int_tc/1');  L('Latch/2','Int_q0/1');

% outlet tail
L('Clock/1','Tail/1');    L('ValveD/1','Tail/2');    L('Int_q0/1','Tail/3');
L('Int_tc/1','Tail/4');   L('Int_tail/1','Tail/5');  L('Int_n/1','Tail/6');
L('Int_Q/1','Tail/7');    L('S/1','Tail/8');
L('Tail/1','Int_tail/1'); L('Tail/2','Stop/1');
L('Discharged/1','Collected/1');  L('Int_tail/1','Collected/2');

% logging; the channel order is documented by alba9_log_index
L('Ref/1','Mux/1');         L('Int_Q/1','Mux/2');      L('Int_V/1','Mux/3');
L('Discharged/1','Mux/4');  L('Plant/6','Mux/5');      L('Int_n/1','Mux/6');
L('Ctrl/3','Mux/7');        L('Plant/5','Mux/8');      L('Plant/4','Mux/9');
L('ValveD/1','Mux/10');     L('Tail/1','Mux/11');      L('Int_tail/1','Mux/12');
L('Collected/1','Mux/13');  L('Ctrl/4','Mux/14');      L('Ctrl/1','Mux/15');
L('Int_ni/1','Mux/16');     L('Int_q0/1','Mux/17');    L('Int_tc/1','Mux/18');
L('Tail/3','Mux/19');       L('Int_latch/1','Mux/20'); L('CmpD/1','Mux/21');
L('Ctrl/5','Mux/22');       L('Ctrl/6','Mux/23');
L('Mux/1','ToWs/1');

% ---------------- solver ----------------
set_param(mdl, ...
    'SolverType','Variable-step','Solver','ode15s', ...
    'RelTol','2e-8','AbsTol','1e-7','MaxStep','0.025', ...
    'ZeroCrossControl','UseLocalSettings','StopTime','120', ...
    'SaveOutput','off','SaveTime','off','SignalLogging','off');

save_system(mdl,slx);
fprintf('Built %s\n',slx);
end

% -------------------------------------------------------------------------
function setChart(mdl,name,code)
%SETCHART Write the body of a MATLAB Function block.
chart = sfroot().find('-isa','Stateflow.EMChart','Path',[mdl '/' name]);
chart.Script = sprintf(code);
end
