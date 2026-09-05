function alba9v3_extend_model(mdl)
% v3 diagnostic additions to the v2 model; unchanged physical RHS and latch.
ch=sfroot().find('-isa','Stateflow.EMChart','Path',[mdl '/Plant']);
ch.Script=sprintf(['function [dx,dv,dV,dQ,Fn,phi,torque,omega,tension,acc,nviol,gidx,stop,g,gap,pressure] = fcn(x,v,V,Q,count,closed)\n' ...
'%%#codegen\n[dx,dv,dV,dQ,Fn,phi,torque,omega,tension,acc,nviol,gidx,stop,g,gap,pressure]=alba9b_plant_block(x,v,V,Q,count,closed);\nend\n']);
add=@(src,n,pos,varargin)add_block(src,[mdl '/' n],'Position',pos,varargin{:});
L=@(a,b)add_line(mdl,a,b,'autorouting','on');
delete_line(mdl,'V_mL/1','Ctrl/3');
add('simulink/Sources/Constant','V0_est_mL',[25 940 110 970],'Value','V0_m3*1e6');
add('simulink/Math Operations/Sum','V_est_mL',[180 940 210 975],'Inputs','+-');
L('V0_est_mL/1','V_est_mL/1'); L('intQ_mL/1','V_est_mL/2'); L('V_est_mL/1','Ctrl/3');
% The reference scheduler and valve threshold must also use measured-flow
% integration. True V remains available only to the plant and truth logging.
delete_line(mdl,'Discharged/1','Ref/2');L('intQ_mL/1','Ref/2');
delete_line(mdl,'Discharged/1','ValveCmp/1');L('intQ_mL/1','ValveCmp/1');
add('simulink/User-Defined Functions/MATLAB Function','Completion',[550 950 730 1120]);
ch=sfroot().find('-isa','Stateflow.EMChart','Path',[mdl '/Completion']);
ch.Script=sprintf(['function [quiet,notquiet,done,term] = fcn(n,omega,acc,vc,dwell,t,tc,q0,drained,discharged,S)\n' ...
'%%#codegen\nquiet=double(abs(n)<=S(22) && abs(omega)<=S(23) && abs(acc)<=S(24));\n' ...
'notquiet=1-quiet;\nT=inf; if q0>0, T=2*S(3)/q0; end\n' ...
'done=double(vc>0.5 && t>=tc+T && drained>=S(3)-1e-6 && dwell>=S(25));\n' ...
'term=double(S(20)>0.5 && S(1)-S(2)-discharged<=S(21));\nend\n']);
add('simulink/Continuous/Integrator','Int_quiet',[810 960 860 1000], ...
    'InitialCondition','0','AbsoluteTolerance','1e-8','ExternalReset','rising');
L('Int_n/1','Completion/1');L('Plant/8','Completion/2');L('Plant/10','Completion/3');
L('ValveD/1','Completion/4');L('Int_quiet/1','Completion/5');L('Clock/1','Completion/6');
L('Int_tc/1','Completion/7');L('Int_q0/1','Completion/8');L('Int_tail/1','Completion/9');
L('intQ_mL/1','Completion/10');L('S/1','Completion/11');
% Separate the reset predicate from the timer-dependent done calculation.
% A single MATLAB Function made an artificial direct-feedthrough reset loop.
add('simulink/User-Defined Functions/MATLAB Function','QuietTest',[260 960 430 1040]);
ch=sfroot().find('-isa','Stateflow.EMChart','Path',[mdl '/QuietTest']);
ch.Script=sprintf(['function [quiet,notquiet]=fcn(n,w,a,S)\n%%#codegen\n' ...
 'quiet=double(abs(n)<=S(22) && abs(w)<=S(23) && abs(a)<=S(24));\nnotquiet=1-quiet;\nend\n']);
L('Int_n/1','QuietTest/1');L('Plant/8','QuietTest/2');L('Plant/10','QuietTest/3');L('S/1','QuietTest/4');
L('QuietTest/1','Int_quiet/1');L('QuietTest/2','Int_quiet/2');
delete_line(mdl,'Tail/2','StopDone/1'); L('Completion/3','StopDone/1');
% Completion is an event as well: locate the continuous tail-end surface.
% Otherwise a MATLAB Function Boolean waits for the next major solver step.
add('simulink/User-Defined Functions/MATLAB Function','TailEndSurface',[480 1190 650 1270]);
ch=sfroot().find('-isa','Stateflow.EMChart','Path',[mdl '/TailEndSurface']);
ch.Script=sprintf(['function margin=fcn(t,tc,q0,closed,S)\n%%#codegen\n' ...
 'margin=-1; if closed>0.5 && q0>0, margin=t-tc-2*S(3)/q0; end\nend\n']);
add('simulink/Sources/Constant','TailZero',[500 1300 550 1330],'Value','0');
add('simulink/Logic and Bit Operations/Relational Operator','TailEndCrossing',[710 1190 750 1230], ...
 'Operator','>=','ZeroCross','on');
add('simulink/Logic and Bit Operations/Logical Operator','DoneLocated',[820 1200 860 1250],'Operator','AND','Inputs','2');
L('Clock/1','TailEndSurface/1');L('Int_tc/1','TailEndSurface/2');L('Int_q0/1','TailEndSurface/3');
L('ValveD/1','TailEndSurface/4');L('S/1','TailEndSurface/5');
L('TailEndSurface/1','TailEndCrossing/1');L('TailZero/1','TailEndCrossing/2');
delete_line(mdl,'Completion/3','StopDone/1');
L('Completion/3','DoneLocated/1');L('TailEndCrossing/1','DoneLocated/2');L('DoneLocated/1','StopDone/1');
% Locate guard surfaces with Simulink zero crossings, rather than relying on
% a Boolean computed inside MATLAB Function only at the next major step.
delete_line(mdl,'Plant/13','StopGuard/1');
add('simulink/Sources/Constant','GuardZero',[30 1140 90 1170],'Value','0');
add('simulink/Logic and Bit Operations/Relational Operator','GuardCrossing',[170 1110 220 1150], ...
 'Operator','<=','ZeroCross','on');
add('simulink/User-Defined Functions/MATLAB Function','AnyGuard',[290 1100 420 1150]);
ch=sfroot().find('-isa','Stateflow.EMChart','Path',[mdl '/AnyGuard']);
ch.Script=sprintf('function stop=fcn(flags)\n%%#codegen\nstop=any(flags);\nend\n');
L('Plant/14','GuardCrossing/1');L('GuardZero/1','GuardCrossing/2');
L('GuardCrossing/1','AnyGuard/1');L('AnyGuard/1','StopGuard/1');
set_param([mdl '/Mux'],'Inputs','40');
% Mux input 33 is a 7-vector; remaining inputs are scalar: 46 columns total.
L('Plant/14','Mux/33');L('Plant/15','Mux/34');L('Plant/16','Mux/35');
L('V_est_mL/1','Mux/36');L('Int_quiet/1','Mux/37');L('Completion/3','Mux/38');
L('Completion/4','Mux/39');L('Completion/1','Mux/40');
set_param(mdl,'SaveState','on','StateSaveName','xout','SaveFinalState','on', ...
 'FinalStateName','xFinal','SaveOperatingPoint','on','ReturnWorkspaceOutputs','on');
Simulink.Annotation(mdl,['v3: full frozen dynamics; no imposed flow or speed. ' ...
 'Controller uses ideal Q and V0-integral(Q). ' ...
 'bench effective resistance unchanged. Guards retain individual SI units. ' ...
 'Done requires actual motion settled for 0.1 s and tail ended; n=0 is not power-off.']);
end
