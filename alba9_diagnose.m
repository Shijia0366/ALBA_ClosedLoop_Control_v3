function d = alba9_diagnose(t,y,S)
%ALBA9_DIAGNOSE Stage-9 v2 diagnostics.
%
%   Everything here is measured against DESIGN targets on a reduced surrogate of
%   the conditional reference model. None of it is experimental validation and
%   none of these limits is a medical or regulatory criterion.
%
%   What this CANNOT check: running torque, electrical phase, rotor speed, cable
%   travel and stall margin. ROM3 has no motor states - the fast motor
%   coordinates were eliminated when it was identified. Those checks belong to
%   the replay of the SAME controller on the full five-state reference model.
%
%   v2 CHANGES, all of them fixes to the MEASUREMENT rather than the control:
%     * phases are separated (startup / plateau / deceleration / closed) and each
%       error statistic is TIME-WEIGHTED. v1 used an unweighted sample RMSE,
%       which over-weights whatever the variable-step solver sampled densely.
%     * the open and closed phases are integrated separately, so nothing is ever
%       interpolated or smoothed across the flow discontinuity at closure.
%     * the command-acceleration envelope is checked separately against 6000
%       STEP/s^2 while open and 2000 STEP/s^2 while braking, instead of passing
%       everything against the larger of the two.
%     * the old "mass_balance_error" was the construction identity
%       collected = discharged + tail. It is kept under an honest name and three
%       genuinely independent balances are added.
%     * valve transitions, the actual captured q0, and the saturation and slew
%       states are reported as events, not only as summary maxima.
ix = alba9_log_index();
P  = alba9_const();

d = struct();
d.target_mL = S.target_mL;
d.mode      = S.mode;
d.mode_name = modeName(S.mode);

% ---------------------------------------------------------------- valve events
vc  = y(:,ix.valve_closed) > 0.5;
raw = y(:,ix.cmp_raw)      > 0.5;
d.valve_open_to_closed = sum(diff(double(vc)) > 0);
d.valve_closed_to_open = sum(diff(double(vc)) < 0);
d.raw_comparator_open_to_closed = sum(diff(double(raw)) > 0);
d.raw_comparator_closed_to_open = sum(diff(double(raw)) < 0);
d.valve_single_transition = (d.valve_open_to_closed == 1) && (d.valve_closed_to_open == 0);
if ~any(vc)
    error('alba9_diagnose:neverClosed','The valve never closed in this run');
end
ic = find(vc,1);

% Rows belonging to each side of the discontinuity, de-duplicated WITHIN the
% side. Locating the zero crossing makes the solver emit several rows at the
% same timestamp; keeping the last row on each side preserves both one-sided
% limits and never mixes them.
iOpen   = sideIndex(t,~vc);
iClosed = sideIndex(t, vc);

d.t_close_s        = t(ic);
d.q_close_left_mL_s= y(iOpen(end),ix.Q);      % last OPEN-side sample of Q
d.q0_captured_mL_s = y(end,ix.q0);            % what the model actually used
d.q0_capture_note  = sprintf( ...
    ['q0 is a first-order tracker with time constant %.3g s frozen at closure, ' ...
     'not the exact left limit of Q. Difference from the last open-side sample: ' ...
     '%.3e mL/s.'],S.tau_q0,d.q0_captured_mL_s-d.q_close_left_mL_s);
d.latch_lag_s      = 0.5/S.latch_gain;
d.latch_lag_volume_mL = d.latch_lag_s*d.q_close_left_mL_s;
d.cycle_end_s      = t(end);

% motor stop and tail end
k = find(vc & y(:,ix.n) <= 1e-6,1);
d.motor_stop_s = ternary(isempty(k),NaN,t(max(k,1)));
k = find(y(:,ix.drained) >= S.tail_vol_mL-1e-6,1);
d.tail_end_s = ternary(isempty(k),NaN,t(max(k,1)));
d.tail_duration_s = 2*S.tail_vol_mL/max(d.q0_captured_mL_s,eps);
d.tail_duration_basis = ['T = 2*V_tail/q0. Reproduces the user-reported ~1.33 s only ' ...
    'near q0 = 3 mL/s; at any other closing flow this is an EXTRAPOLATION of the ' ...
    'assumed half-cosine rule, not an observed drain time.'];

% ---------------------------------------------------------------- volumes
d.collected_mL         = y(end,ix.collected);
d.bladder_discharge_mL = y(end,ix.discharged);
d.tail_drained_mL      = y(end,ix.drained);
d.volume_error_mL      = d.collected_mL - S.target_mL;
d.volume_error_is_bookkeeping_only = true;
d.volume_error_note = ['Near zero BY CONSTRUCTION: the valve closes on an ideal volume ' ...
    'event and the tail delivered equals the tail compensated. This is not a ' ...
    'demonstration of volumetric control accuracy. See tail_sensitivity.'];
tails = [1.0 1.5 2.0 2.5 3.0];
d.tail_sensitivity = struct('actual_tail_mL',num2cell(tails), ...
                            'final_volume_error_mL',num2cell(tails-S.tail_comp_mL));

% ---------------------------------------------------------------- phases
Qref = y(:,ix.Qref);
ph_start = iOpen(t(iOpen) <= S.rise_s);
ph_plat  = iOpen(t(iOpen) >  S.rise_s & Qref(iOpen) >= S.q_plateau-1e-6);
ph_dec   = iOpen(Qref(iOpen) <  S.q_plateau-1e-6 & t(iOpen) > S.rise_s);

d.phase_definition = struct( ...
    'startup',     sprintf('valve open and t <= rise_s = %.3g s',S.rise_s), ...
    'plateau',     'valve open, past startup, reference still at the design plateau', ...
    'deceleration','valve open, reference below the plateau', ...
    'closed',      'after the latched valve closure');
d.startup      = phaseStats(t,y,ph_start,ix,S);
d.plateau      = phaseStats(t,y,ph_plat, ix,S);
d.deceleration = phaseStats(t,y,ph_dec,  ix,S);

d.plateau_err_pct   = d.plateau.max_abs_rel_error_pct;
d.plateau_rmse_mL_s = d.plateau.rmse_time_weighted_mL_s;
d.plateau_rmse_unweighted_mL_s = d.plateau.rmse_unweighted_mL_s;
d.plateau_min_mL_s  = d.plateau.min_Q_mL_s;
d.plateau_max_mL_s  = d.plateau.max_Q_mL_s;

% ---------------------------------------------------------------- actuator
d.max_n_steps_s = max(y(:,ix.n));
d.max_accel_open_steps_s2   = maxAbs(y(iOpen,  ix.n_dot));
d.max_accel_closed_steps_s2 = maxAbs(y(iClosed,ix.n_dot));
d.accel_open_ok   = d.max_accel_open_steps_s2   <= S.a_max  + 1e-6;
d.accel_closed_ok = d.max_accel_closed_steps_s2 <= S.a_brake+ 1e-6;
d.n_saturated_time_s = phaseTime(t,iOpen,y(iOpen,ix.sat_flag) ~= 0);
d.slew_limited_time_s = phaseTime(t,iOpen,y(iOpen,ix.slew_flag) > 0.5);
d.final_ni = y(end,ix.ni);
d.max_abs_ni = maxAbs(y(:,ix.ni));

% ---------------------------------------------------------------- ROM domain
d.max_F_N  = max(y(:,ix.F));
d.min_keff = min(y(:,ix.keff));
d.F_within_audited_118N = d.max_F_N <= 118;
d.audited_range_note = ['118 N is the upper end of the range over which the LEARNED ' ...
    'stiffness was checked positive. It is not a motor force limit and not a ' ...
    'saturation in the model; exceeding it means the ROM is being used outside ' ...
    'its audited domain.'];
d.min_V_mL = min(y(:,ix.V));
d.V_above_Vmin = d.min_V_mL > P.Vmin*1e6;

% ---------------------------------------------------------------- balances
tube_mL = P.pipe_area*P.pipe_length*1e6;
intQ    = trapzOn(t,y(:,ix.Q),     iOpen);
intTail = trapzOn(t,y(:,ix.q_tail),iClosed);
d.balance = struct();
d.balance.accounting_identity_error_mL = maxAbs( ...
    y(:,ix.collected) - (y(:,ix.discharged)+y(:,ix.drained)));
d.balance.accounting_identity_note = ['This is the CONSTRUCTION identity ' ...
    'collected = discharged + tail. It cannot fail and proves nothing physical. ' ...
    'Kept only so a regression would show up.'];
d.balance.valve_flow_integral_mL = intQ;
d.balance.bladder_decrease_mL    = y(iOpen(1),ix.V) - y(iOpen(end),ix.V);
d.balance.valve_flow_vs_bladder_mL = intQ - d.balance.bladder_decrease_mL;
d.balance.tail_integral_mL = intTail;
d.balance.tail_vs_drained_mL = intTail - y(end,ix.drained);
d.balance.tube_initial_mL = tube_mL;
d.balance.total_inventory_error_mL = maxAbs( ...
    (y(:,ix.V) + (tube_mL - y(:,ix.drained)) + y(:,ix.collected)) - (P.V0*1e6 + tube_mL));
d.balance.note = ['valve_flow_vs_bladder and tail_vs_drained are genuine numerical ' ...
    'checks: an integral recomputed from the log against a state the solver ' ...
    'integrated. total_inventory follows from those two plus the accounting ' ...
    'identity, and the tail is never deducted from the bladder a second time.'];

% ---------------------------------------------------------------- acceptance
d.meets_plateau_5pct = d.plateau_err_pct <= 5;
d.meets_volume_2mL   = abs(d.volume_error_mL) <= 2;
end

% -------------------------------------------------------------------------
function idx = sideIndex(t,mask)
%SIDEINDEX Rows on one side of the closure, de-duplicated within that side.
idx = find(mask);
if isempty(idx), return; end
keep = [diff(t(idx)) > 0; true];
idx = idx(keep);
end

function s = phaseStats(t,y,idx,ix,S)
s = struct('samples',numel(idx),'duration_s',0,'min_Q_mL_s',NaN,'max_Q_mL_s',NaN, ...
    'rmse_time_weighted_mL_s',NaN,'rmse_unweighted_mL_s',NaN, ...
    'max_abs_error_mL_s',NaN,'max_abs_rel_error_pct',NaN,'t_start_s',NaN,'t_end_s',NaN);
if numel(idx) < 2, return; end
tt = t(idx); e = y(idx,ix.e_track); q = y(idx,ix.Q);
s.t_start_s = tt(1); s.t_end_s = tt(end);
s.duration_s = tt(end)-tt(1);
s.min_Q_mL_s = min(q); s.max_Q_mL_s = max(q);
s.rmse_time_weighted_mL_s = sqrt(trapz(tt,e.^2)/s.duration_s);
s.rmse_unweighted_mL_s = sqrt(mean(e.^2));
s.max_abs_error_mL_s = max(abs(e));
s.max_abs_rel_error_pct = 100*max(abs(q/S.q_plateau-1));
end

function v = maxAbs(x)
if isempty(x), v = NaN; else, v = max(abs(x)); end
end

function v = trapzOn(t,x,idx)
if numel(idx) < 2, v = 0; else, v = trapz(t(idx),x(idx)); end
end

function s = phaseTime(t,idx,mask)
if numel(idx) < 2 || ~any(mask), s = 0; return; end
tt = t(idx);
dt = [diff(tt); 0];
s = sum(dt(mask));
end

function v = ternary(c,a,b)
if c, v = a; else, v = b; end
end

function s = modeName(m)
switch round(m)
    case 1, s = 'fixed step rate';
    case 2, s = 'pure PI';
    otherwise, s = 'geometric FF + PI';
end
end
