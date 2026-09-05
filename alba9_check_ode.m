function out = alba9_check_ode(S,opts)
%ALBA9_CHECK_ODE Independent ode15s realisation of the same closed loop.
%
%   It exists to cross-check the Simulink model, not to replace it. It calls the
%   SAME frozen plant / controller / reference / tail functions, so a discrepancy
%   is a framework bug, not physics.
%
%   It is deliberately NOT a replica of the Simulink implementation:
%     * the closing volume is found by a terminal EVENT, so the closing time and
%       the pre-closure flow are exact to solver tolerance. q0 here is the true
%       LEFT LIMIT of Q at the closing instant, captured before Q is zeroed.
%     * it therefore has none of the Simulink model's 0.5/latch_gain closing lag
%       and none of its tau_q0 tracking lag.
%   Comparing the two is how those two implementation choices get quantified
%   instead of assumed harmless.
%
%   The unilateral contact condition is applied the plain-ODE way here
%   (clamp_contact = true), because a single ODE right-hand side has no
%   saturated integrator to enforce it.
if nargin < 1 || isempty(S), S = alba9_scenario(300,3); end
if nargin < 2, opts = struct(); end
tight = isfield(opts,'tight') && opts.tight;

vec = packS(S);
P = alba9_const();
V0mL = P.V0*1e6;

if tight
    rtol = 1e-10; atol = [1e-11 1e-9 1e-9 1e-6 1e-6 1e-11]; maxstep = 0.0125;
else
    rtol = 2e-8;  atol = [1e-9 1e-7 1e-7 1e-4 1e-4 1e-9];   maxstep = 0.025;
end

% states: [f; Vm; Qm; n; ni; drained]
z0 = [0; V0mL; 0; 0; 0; 0];
o1 = odeset('RelTol',rtol,'AbsTol',atol,'MaxStep',maxstep, ...
            'Events',@(t,z)closeEvent(t,z,vec,V0mL));
sol1 = ode15s(@(t,z)rhs(t,z,vec,P,0,0,0),[0 400],z0,o1);
if isempty(sol1.xe)
    error('alba9_check_ode:noClosure','The valve closing volume was never reached');
end
tc = sol1.xe(1);
ze = sol1.ye(:,1);
q0_exact = ze(3);                       % TRUE left limit of Q, before any reset

z2 = ze; z2(3) = 0;                     % ideal instantaneous closure
T_tail = 2*vec(3)/q0_exact;
o2 = odeset('RelTol',rtol,'AbsTol',atol,'MaxStep',min(maxstep,0.01));
sol2 = ode15s(@(t,z)rhs(t,z,vec,P,1,q0_exact,tc),[tc tc+max(T_tail,2.5)+0.5],z2,o2);

t1 = sol1.x(:); t2 = sol2.x(:);
z1 = sol1.y.';  zz2 = sol2.y.';
t = [t1; t2(2:end)];
z = [z1; zz2(2:end,:)];
closed = [zeros(numel(t1),1); ones(numel(t2)-1,1)];

out = struct();
out.t = t;
out.f = z(:,1);            out.F_N = max(z(:,1)*100,0);
out.V_mL = z(:,2);         out.Q_mL_s = z(:,3);
out.n_steps_s = z(:,4);    out.ni = z(:,5);
out.drained_mL = z(:,6);
out.discharged_mL = V0mL - z(:,2);
out.collected_mL = out.discharged_mL + out.drained_mL;
out.valve_closed = closed;
out.q_outlet_mL_s = out.Q_mL_s;
out.q_outlet_mL_s(closed>0.5) = arrayfun(@(k) ...
    alba9_tail(t(k),1,q0_exact,tc,z(k,6),z(k,4),z(k,3),vec), find(closed>0.5));
out.Qref_mL_s = arrayfun(@(k) ...
    alba9_reference(t(k),out.discharged_mL(k),closed(k),vec),(1:numel(t)).').';
out.Qref_mL_s = out.Qref_mL_s(:);

out.t_close = tc;
out.q0_exact_left_limit = q0_exact;
out.tail_duration_s = T_tail;
out.sol_open = sol1;
out.sol_closed = sol2;
out.scenario = S;
out.tight = tight;
out.final_collected_mL = out.collected_mL(end);
end

% -------------------------------------------------------------------------
function dz = rhs(t,z,vec,P,closed,q0,tc)
f = z(1); Vm = z(2); Qm = z(3); n = z(4); ni = z(5); drained = z(6);
discharged = P.V0*1e6 - Vm;
Qref = alba9_reference(t,discharged,closed,vec);
[n_in,dni] = alba9_controller(Qref,Qm,Vm,n,ni,closed,vec);
[df,dVm,dQm] = alba9_plant(f,Vm,Qm,n,closed,P);      % clamp_contact = true here
if closed > 0.5
    q_tail = alba9_tail(t,1,q0,tc,drained,n,Qm,vec);
else
    q_tail = 0;
end
% the command state is clamped at its limits the way the saturated integrator is
if (n >= vec(12) && n_in > 0) || (n <= 0 && n_in < 0)
    n_in = 0;
end
dz = [df; dVm; dQm; n_in; dni; q_tail];
end

function [value,isterminal,direction] = closeEvent(~,z,vec,V0mL)
value = (vec(1)-vec(2)) - (V0mL - z(2));   % remaining volume before closure
isterminal = 1;
direction = -1;
end

function vec = packS(S)
[~,vec] = alba9_scenario(S.target_mL,S.mode);
names = {'target_mL','tail_comp_mL','tail_vol_mL','q_plateau','q_close', ...
         'rise_s','R0_mL','r_hold_mL','Kp','Ki','Kaw','n_max','a_max', ...
         'a_brake','tau_slew','mode','n_fixed','latch_gain','tau_q0'};
for i = 1:numel(names)
    vec(i) = S.(names{i});
end
end
