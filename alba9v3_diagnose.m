function d=alba9v3_diagnose(t,y,S,stopTime)
% v3: dimensional margins, explicit completion, no state clipping.
ix=alba9b_log_index();P=alba9_const();B=alba9b_const();
vc=y(:,10)>.5;io=side(t,~vc);ic=side(t,vc);
d=struct('target_mL',S.target_mL,'mode',S.mode,'mode_name',modeName(S.mode), ...
 'case_definition','Upper-feasible / near-limit case under nominal model assumptions', ...
 'controller_type','continuous-time PI; no discrete sampling period', ...
 'resistance_name','bench effective resistance','Rh',P.Rh,'hardware_validated',false);
d.valve_open_to_closed=sum(diff(double(vc))>0);d.valve_closed_to_open=sum(diff(double(vc))<0);
d.valve_single_transition=d.valve_open_to_closed==1 && d.valve_closed_to_open==0;
names=alba9b_guard_names();units={'rad','m','m^3','m','rad/s','N','m^3/s'};
d.guards=struct('name',{},'unit',{},'minimum_margin',{},'time_s',{},'violated',{});
for j=1:7
 idx=(1:numel(t))';if j==7,idx=io;end
 [val,k]=min(y(idx,ix.guards(j)));
 d.guards(j)=struct('name',names{j},'unit',units{j},'minimum_margin',val,'time_s',t(idx(k)),'violated',val<=0);
end
d.guard_violated=any([d.guards.violated]);
d.boundary_events=struct('name',{},'unit',{},'margin',{},'time_s',{},'Q_actual_mL_s',{},'force_N',{},'phase_deg',{},'actual_omega_rad_s',{},'previous_sample_s',{});
for j=1:7
 idx=find(y(:,ix.guards(j))<=0 & (j~=7 | ~vc),1);
 if ~isempty(idx)
 d.boundary_events(end+1)=struct('name',names{j},'unit',units{j},'margin',y(idx,ix.guards(j)), ...
 'time_s',t(idx),'Q_actual_mL_s',y(idx,ix.Q),'force_N',y(idx,ix.Fn),'phase_deg',y(idx,ix.phi)*180/pi, ...
 'actual_omega_rad_s',y(idx,ix.omega),'previous_sample_s',t(max(1,idx-1))); %#ok<AGROW>
 end
end
d.completed=~d.guard_violated && y(end,ix.done)>.5 && d.valve_single_transition;
d.status=v3_classify_status(d.completed,d.guard_violated,t(end)>=stopTime-1e-6);
d.stop_time_s=t(end);d.t_close_s=first(t,vc);d.q_close_left_mL_s=NaN;
if ~isempty(ic),d.q_close_left_mL_s=y(io(end),ix.Q);end
d.q0_captured_mL_s=y(end,ix.q0);d.tail_duration_s=NaN;d.tail_end_s=NaN;
if any(vc)
 d.tail_duration_s=2*S.tail_vol_mL/d.q0_captured_mL_s;te=y(end,ix.tclose)+d.tail_duration_s;
 if te<=t(end)+1e-7,d.tail_end_s=te;end
end
active=find(y(:,ix.n)>S.stop_n_steps_s,1);
d.command_start_s=first(t,y(:,ix.n)>S.stop_n_steps_s);
d.command_stop_s=NaN;d.mechanical_stop_s=NaN;d.mechanical_stop_confirmed_s=NaN;
if ~isempty(active)
 k=find((1:numel(t))'>active & y(:,ix.n)<=S.stop_n_steps_s,1);
 if ~isempty(k),d.command_stop_s=t(k);end
 % Report mechanical settling independently of n, using the last violation
 % of the ACTUAL velocity/acceleration thresholds, not the command timer.
 moving=abs(y(:,ix.omega))>S.stop_omega_rad_s | abs(y(:,ix.acc))>S.stop_acc_m_s2;
 k=find(moving,1,'last');
 if ~isempty(k) && k<numel(t) && t(end)-t(k+1)>=S.stop_dwell_s
   d.mechanical_stop_s=t(k+1);d.mechanical_stop_confirmed_s=t(k+1)+S.stop_dwell_s;
   d.mechanical_event_bracket_s=[t(k) t(k+1)];
 end
end
d.stop_criterion=struct('command_steps_s',S.stop_n_steps_s,'actual_omega_rad_s',S.stop_omega_rad_s, ...
 'actual_acc_m_s2',S.stop_acc_m_s2,'continuous_dwell_s',S.stop_dwell_s,'zero_command_is_not_power_off',true);
d.max_abs_phase_deg=max(abs(y(:,ix.phi)))*180/pi;
d.electrical_phase_boundary_margin_deg=90-d.max_abs_phase_deg;
d.phase_margin_definition='Distance to electrical-phase model boundary; NOT closed-loop stability phase margin';
[d.max_abs_torque_Nm,k]=max(abs(y(:,ix.torque)));d.peak_torque_time_s=t(k);
d.torque_envelope_Nm=B.torque_cap_Nm;d.torque_utilisation_pct=100*d.max_abs_torque_Nm/B.torque_cap_Nm;
d.torque_margin_Nm=B.torque_cap_Nm-d.max_abs_torque_Nm;
d.max_rotor_rpm=max(abs(y(:,ix.omega)))*60/(2*pi);d.max_n_steps_s=max(y(:,ix.n));
d.max_contact_force_N=max(y(:,ix.Fn));d.min_bladder_mL=min(y(:,ix.V));
d.min_signed_overlap_m=min(y(:,ix.overlap_raw));d.maximum_radial_gap_um=max(0,-d.min_signed_overlap_m)/(2*pi)*1e6;
d.recontact_s=first(t,t>1e-5 & y(:,ix.Fn)>1e-3);d.min_cable_tension_N=min(y(:,ix.tension));
d.max_cable_travel_m=max(y(:,ix.x));d.travel_margin_m=B.xmax-d.max_cable_travel_m;
d.postclose_peak_torque_Nm=mx(y(ic,ix.torque));d.postclose_peak_force_N=mx(y(ic,ix.Fn));
d.postclose_travel_m=NaN;d.postclose_pressure_increase_Pa=NaN;
if ~isempty(ic)
 d.postclose_travel_m=max(y(ic,ix.x))-y(ic(1),ix.x);
 d.postclose_pressure_increase_Pa=max(y(ic,ix.pressure))-y(ic(1),ix.pressure);
end
d.collected_mL=y(end,ix.collected);d.volume_error_mL=d.collected_mL-S.target_mL;
d.meets_volume_2mL=abs(d.volume_error_mL)<=2 && d.completed;d.volume_error_is_bookkeeping_only=true;
d.startup=stats(t,y,io(t(io)<=S.rise_s),ix,S);
d.plateau=stats(t,y,io(t(io)>S.rise_s & y(io,ix.Qref)>=S.q_plateau-1e-6),ix,S);
d.deceleration=stats(t,y,io(t(io)>S.rise_s & y(io,ix.Qref)<S.q_plateau-1e-6),ix,S);
d.plateau_rmse_mL_s=d.plateau.rmse_time_weighted_mL_s;d.plateau_err_pct=d.plateau.max_abs_rel_error_pct;
d.meets_plateau_5pct=d.plateau_err_pct<=5;
d.max_accel_open_steps_s2=mx(y(io,ix.n_dot));d.max_accel_closed_steps_s2=mx(y(ic,ix.n_dot));
d.accel_open_ok=d.max_accel_open_steps_s2<=S.a_max+1e-6;
d.accel_closed_ok=isempty(ic)||d.max_accel_closed_steps_s2<=S.a_brake+1e-6;
d.upper_saturation_s=duration(t,io,y(:,ix.sat_flag)>0);d.lower_saturation_s=duration(t,io,y(:,ix.sat_flag)<0);
d.slew_limited_s=duration(t,io,y(:,ix.slew_flag)>.5);d.final_integral_command=y(end,ix.ni);
tube=P.pipe_area*P.pipe_length*1e6;
d.balance=struct('state_integral_vs_bladder_mL',mx(y(:,ix.intQ)-y(:,ix.discharged)), ...
 'flow_integral_vs_bladder_mL',integ(t,y(:,ix.Q),io)-(y(io(1),ix.V)-y(io(end),ix.V)), ...
 'tail_integral_vs_drained_mL',integ(t,y(:,ix.q_tail),ic)-y(end,ix.drained), ...
 'total_inventory_identity_error_mL',mx(y(:,ix.V)+tube-y(:,ix.drained)+y(:,ix.collected)-(P.V0*1e6+tube)), ...
 'V_est_vs_state_mL',mx(y(:,ix.V_est)-y(:,ix.V)),'minimum_pipe_inventory_mL',tube-max(y(:,ix.drained)));
end
function i=side(t,m),i=find(m);if ~isempty(i),i=i([diff(t(i))>0;true]);end,end
function v=first(t,m),k=find(m,1);if isempty(k),v=NaN;else,v=t(k);end,end
function v=mx(x),if isempty(x),v=NaN;else,v=max(abs(x));end,end
function v=integ(t,x,i),if numel(i)<2,v=0;else,v=trapz(t(i),x(i));end,end
function v=duration(t,i,m),if numel(i)<2,v=0;else,dt=[diff(t(i));0];v=sum(dt(m(i)));end,end
function s=stats(t,y,i,ix,S)
s=struct('duration_s',0,'min_Q_mL_s',NaN,'max_Q_mL_s',NaN,'rmse_time_weighted_mL_s',NaN, ...
 'max_abs_rel_error_pct',NaN,'overshoot_pct',NaN,'settled_in_5pct_band_s',NaN);
if numel(i)<2,return;end
tt=t(i);q=y(i,ix.Q);e=y(i,ix.e_track);s.duration_s=tt(end)-tt(1);
s.min_Q_mL_s=min(q);s.max_Q_mL_s=max(q);s.rmse_time_weighted_mL_s=sqrt(trapz(tt,e.^2)/s.duration_s);
s.max_abs_rel_error_pct=100*max(abs(q/S.q_plateau-1));s.overshoot_pct=100*(max(q)/S.q_plateau-1);
bad=find(abs(q/S.q_plateau-1)>.05,1,'last');
if isempty(bad),s.settled_in_5pct_band_s=tt(1);elseif bad<numel(i),s.settled_in_5pct_band_s=tt(bad+1);end
end
function n=modeName(m),nn={'fixed step rate','pure PI','geometric FF + PI'};n=nn{round(m)};end
