function c = alba9_compare(t,y,ref,label)
%ALBA9_COMPARE Full-trajectory comparison against a reference solution.
%
%   Compares the WHOLE trajectory, not only the final volume: valve flow,
%   bladder volume, contact force and motor command, on the open phase and the
%   closed phase separately.
%
%   The two runs close the valve at slightly different instants, so a window
%   around closure is EXCLUDED rather than compared: interpolating across a flow
%   discontinuity would manufacture an error that means nothing. The width of
%   that window is reported, so the exclusion is visible instead of silent.
%
%   ref is either an alba9_check_ode output (with dense sol_open / sol_closed)
%   or another (t,y) Simulink log wrapped by alba9_wrap_log.
%
%   THIS IS NUMERICAL VERIFICATION, NOT EXPERIMENTAL VALIDATION. Agreement means
%   two implementations of the same equations agree.
if nargin < 4, label = 'reference'; end
ix = alba9_log_index();

tc_a = t(find(y(:,ix.valve_closed) > 0.5,1));
tc_b = ref.t_close;
guard = max(abs(tc_a-tc_b)*2, 1e-6);

openMask   = t < min(tc_a,tc_b) - guard;
closedMask = t > max(tc_a,tc_b) + guard;

c = struct('label',label,'t_close_difference_s',abs(tc_a-tc_b), ...
           'excluded_window_s',2*guard, ...
           'open_samples',sum(openMask),'closed_samples',sum(closedMask));

[c.open, c.open_note]     = cmpPhase(t(openMask),  y(openMask,:),  ref.sol_open,  ix);
[c.closed,c.closed_note]  = cmpPhase(t(closedMask),y(closedMask,:),ref.sol_closed,ix);

c.worst_flow_mL_s   = max([c.open.dQ_mL_s,      c.closed.dQ_mL_s]);
c.worst_volume_mL   = max([c.open.dV_mL,        c.closed.dV_mL]);
c.worst_force_N     = max([c.open.dF_N,         c.closed.dF_N]);
c.worst_command_sps = max([c.open.dn_steps_s,   c.closed.dn_steps_s]);
c.final_collected_difference_mL = abs(y(end,ix.collected) - ref.collected_mL(end));
end

% -------------------------------------------------------------------------
function [s,note] = cmpPhase(tq,yq,sol,ix)
note = '';
s = struct('samples',numel(tq),'dQ_mL_s',0,'dV_mL',0,'dF_N',0,'dn_steps_s',0);
if isempty(tq)
    note = 'no samples outside the excluded closure window';
    return;
end
tq = min(max(tq,sol.x(1)),sol.x(end));
z = deval(sol,tq);                       % [f; Vm; Qm; n; ni; drained]
s.dQ_mL_s    = max(abs(yq(:,ix.Q) - z(3,:).'));
s.dV_mL      = max(abs(yq(:,ix.V) - z(2,:).'));
s.dF_N       = max(abs(yq(:,ix.F) - max(z(1,:).'*100,0)));
s.dn_steps_s = max(abs(yq(:,ix.n) - z(4,:).'));
end
