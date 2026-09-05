function d = alba9b_diagnose(t,y,S,stopTime)
%ALBA9B_DIAGNOSE Diagnostics for the FULL five-state reference closed loop.
%
%   Reports everything ROM3 could not: electrical phase, torque against the
%   conditional envelope, ACTUAL rotor speed, cable travel and tension, and the
%   model-validity guards.
%
%   A run that stops early because a guard crossed zero is reported as INCOMPLETE
%   with the guard named. Nothing is clipped, smoothed or retried to make it
%   finish. 15 mL/s is treated as a near-limit design case, not as a case that
%   must succeed.
ix = alba9b_log_index();
P  = alba9_const();
B  = alba9b_const();
names = alba9b_guard_names();

d = struct();
d.target_mL = S.target_mL;
d.mode = S.mode;
d.mode_name = modeName(S.mode);
d.plant = 'full five-state reference, ported and pointwise-verified';

vc = y(:,ix.valve_closed) > 0.5;
d.valve_open_to_closed = sum(diff(double(vc)) > 0);
d.valve_closed_to_open = sum(diff(double(vc)) < 0);
d.valve_single_transition = (d.valve_open_to_closed == 1) && (d.valve_closed_to_open == 0);

% ---- did it finish, or did the model leave its declared domain? -------------
[gm,k] = min(y(:,ix.gmin));
d.min_guard_value = gm;
d.min_guard_name  = names{max(round(y(k,ix.gidx)),1)};
d.min_guard_time_s = t(k);
d.guard_violated = gm <= 0;
d.completed = ~d.guard_violated && any(vc) && t(end) < stopTime - 1e-9;
d.stopped_early_at_s = t(end);

% ---- actuator, the whole point of stage B ---------------------------------
d.max_abs_phase_rad = max(abs(y(:,ix.phi)));
d.max_abs_phase_deg = d.max_abs_phase_rad*180/pi;
d.phase_margin_deg  = 90 - d.max_abs_phase_deg;
d.max_abs_torque_Nm = max(abs(y(:,ix.torque)));
d.torque_envelope_Nm = B.torque_cap_Nm;
d.torque_utilisation_pct = 100*d.max_abs_torque_Nm/B.torque_cap_Nm;
d.torque_margin_Nm = B.torque_cap_Nm - d.max_abs_torque_Nm;
d.torque_basis = ['Conditional FLAT envelope, not a measured 12 V running-torque ' ...
    'curve. 0.67 A is an assumed current, not a measured limit. A high utilisation ' ...
    'here is a warning about the assumption as much as about the controller.'];
d.max_rotor_rpm = max(abs(y(:,ix.omega)))*60/(2*pi);
d.max_cable_speed_m_s = max(abs(y(:,ix.v)));
d.max_cable_travel_m = max(y(:,ix.x));
d.travel_limit_m = B.xmax;
d.travel_margin_m = B.xmax - d.max_cable_travel_m;
d.min_cable_tension_N = min(y(:,ix.tension));
d.tension_note = ['A cable cannot push. The frozen guard allows -1e-5 N of ' ...
    'numerical slack; a genuinely negative tension means contact was lost.'];
d.max_contact_force_N = max(y(:,ix.Fn));
d.min_bladder_mL = min(y(:,ix.V));

% ---- flow and volume ------------------------------------------------------
iOpen   = sideIndex(t,~vc);
iClosed = sideIndex(t, vc);
if any(vc)
    ic = find(vc,1);
    d.t_close_s = t(ic);
    d.q_close_left_mL_s = y(iOpen(end),ix.Q);
else
    d.t_close_s = NaN; d.q_close_left_mL_s = NaN;
end
d.q0_captured_mL_s = y(end,ix.q0);
d.collected_mL = y(end,ix.collected);
d.volume_error_mL = d.collected_mL - S.target_mL;
d.volume_error_is_bookkeeping_only = true;
d.tail_duration_s = 2*S.tail_vol_mL/max(d.q0_captured_mL_s,eps);
kk = find(vc & y(:,ix.n) <= 1e-6,1);
d.motor_stop_s = ternary(isempty(kk),NaN,t(max(kk,1)));
kk = find(y(:,ix.drained) >= S.tail_vol_mL-1e-6,1);
d.tail_end_s = ternary(isempty(kk),NaN,t(max(kk,1)));

Qref = y(:,ix.Qref);
ph_start = iOpen(t(iOpen) <= S.rise_s);
ph_plat  = iOpen(t(iOpen) >  S.rise_s & Qref(iOpen) >= S.q_plateau-1e-6);
ph_dec   = iOpen(Qref(iOpen) <  S.q_plateau-1e-6 & t(iOpen) > S.rise_s);
d.startup      = phaseStats(t,y,ph_start,ix,S);
d.plateau      = phaseStats(t,y,ph_plat, ix,S);
d.deceleration = phaseStats(t,y,ph_dec,  ix,S);
d.plateau_err_pct = d.plateau.max_abs_rel_error_pct;
d.plateau_rmse_mL_s = d.plateau.rmse_time_weighted_mL_s;
d.meets_plateau_5pct = d.plateau_err_pct <= 5;

% ---- command --------------------------------------------------------------
d.max_n_steps_s = max(y(:,ix.n));
d.max_accel_open_steps_s2   = maxAbs(y(iOpen,  ix.n_dot));
d.max_accel_closed_steps_s2 = maxAbs(y(iClosed,ix.n_dot));
d.accel_open_ok   = d.max_accel_open_steps_s2   <= S.a_max  + 1e-6;
d.accel_closed_ok = d.max_accel_closed_steps_s2 <= S.a_brake+ 1e-6;

% ---- balances, now with a genuinely independent flow integral -------------
tube_mL = P.pipe_area*P.pipe_length*1e6;
d.balance = struct();
d.balance.state_integral_vs_bladder_mL = maxAbs(y(:,ix.intQ) - y(:,ix.discharged));
d.balance.state_integral_note = ['The flow was integrated a SECOND time as its own ' ...
    'solver state and compared with the bladder decrease. Unlike a post-hoc ' ...
    'trapezoid this is independent of the log sampling.'];
d.balance.trapz_valve_flow_vs_bladder_mL = trapzOn(t,y(:,ix.Q),iOpen) - ...
    (y(iOpen(1),ix.V) - y(iOpen(end),ix.V));
d.balance.tail_vs_drained_mL = trapzOn(t,y(:,ix.q_tail),iClosed) - y(end,ix.drained);
d.balance.total_inventory_error_mL = maxAbs( ...
    (y(:,ix.V) + (tube_mL - y(:,ix.drained)) + y(:,ix.collected)) - (P.V0*1e6 + tube_mL));
d.balance.accounting_identity_error_mL = maxAbs( ...
    y(:,ix.collected) - (y(:,ix.discharged)+y(:,ix.drained)));
end

% -------------------------------------------------------------------------
function idx = sideIndex(t,mask)
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
s.t_start_s = tt(1); s.t_end_s = tt(end); s.duration_s = tt(end)-tt(1);
s.min_Q_mL_s = min(q); s.max_Q_mL_s = max(q);
s.rmse_time_weighted_mL_s = sqrt(trapz(tt,e.^2)/s.duration_s);
s.rmse_unweighted_mL_s = sqrt(mean(e.^2));
s.max_abs_error_mL_s = max(abs(e));
s.max_abs_rel_error_pct = 100*max(abs(q/S.q_plateau-1));
end

function v = maxAbs(x), if isempty(x), v = NaN; else, v = max(abs(x)); end, end
function v = trapzOn(t,x,i), if numel(i)<2, v=0; else, v=trapz(t(i),x(i)); end, end
function v = ternary(c,a,b), if c, v=a; else, v=b; end, end
function s = modeName(m)
switch round(m)
    case 1, s = 'fixed step rate';
    case 2, s = 'pure PI';
    otherwise, s = 'geometric FF + PI';
end
end
