function R = run_alba9b(targets,modes,varargin)
%RUN_ALBA9B Stage-9 B runner: the SAME controller on the FULL reference plant.
%   No retuning in this stage. One independently named result file per case.
%
%   Failures are kept. If a run stops because a model-validity guard crossed
%   zero, that is reported as an incomplete case with the guard named; the
%   controller is not adjusted and no envelope is relaxed to make it finish.
if nargin < 1 || isempty(targets), targets = 300; end
if nargin < 2 || isempty(modes),   modes   = 3; end

opt = struct('tag','B_full_reference','stoptime',120);
scenArgs = {};
k = 1;
while k <= numel(varargin)
    if isfield(opt,varargin{k}), opt.(varargin{k}) = varargin{k+1};
    else, scenArgs(end+1:end+2) = varargin(k:k+1); %#ok<AGROW>
    end
    k = k + 2;
end

here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fileparts(here));
outdir = fullfile(fileparts(here),'results',opt.tag);
if ~isfolder(outdir), mkdir(outdir); end

mdl = 'alba9b_full_loop';
if ~isfile(fullfile(here,[mdl '.slx'])), build_alba9b_model(mdl); end
load_system(fullfile(here,[mdl '.slx']));

P = alba9_const();
assignin('base','V0_m3',P.V0);

R = struct('name',{},'scenario',{},'diag',{},'t',{},'y',{});
for target = targets(:).'
    for mode = modes(:).'
        [S,S_vec] = alba9_scenario(target,mode,scenArgs{:});
        name = sprintf('%dmL_mode%d_%s',target,mode,slug(modeName(mode)));
        assignin('base','S_vec',S_vec);
        set_param(mdl,'StopTime',num2str(opt.stoptime));

        out = sim(mdl);
        L = out.alba9b_log;
        t = L.time;
        y = squeeze(L.signals.values);
        d = alba9b_diagnose(t,y,S,opt.stoptime);

        rec = struct('name',name,'scenario',S,'diag',d,'t',t,'y',y);
        R(end+1) = rec; %#ok<AGROW>
        save(fullfile(outdir,['case_' name '.mat']),'-struct','rec');
        fid = fopen(fullfile(outdir,['case_' name '.log']),'w','n','UTF-8');
        fprintf(fid,'%s\n',jsonencode(d,'PrettyPrint',true));
        fclose(fid);
        report(d);
    end
end
fprintf('\n%d cases written to %s\n',numel(R),outdir);
end

% -------------------------------------------------------------------------
function report(d)
fprintf('\n=== FULL MODEL | %d mL | %s ===\n',d.target_mL,d.mode_name);
if d.completed
    fprintf('  RUN     completed, no guard violated. Closest guard: %s at %.3e (t = %.4f s)\n', ...
        d.min_guard_name,d.min_guard_value,d.min_guard_time_s);
else
    fprintf(2,'  RUN     *** INCOMPLETE *** stopped at %.4f s; guard "%s" reached %.3e\n', ...
        d.stopped_early_at_s,d.min_guard_name,d.min_guard_value);
end
fprintf('  MOTOR   peak torque %.6f N*m = %.1f%% of the %.6f N*m ASSUMED envelope (margin %.6f)\n', ...
    d.max_abs_torque_Nm,d.torque_utilisation_pct,d.torque_envelope_Nm,d.torque_margin_Nm);
fprintf('          peak phase %.2f deg of 90 (margin %.2f deg); peak rotor %.2f rpm; peak cable %.5f m/s\n', ...
    d.max_abs_phase_deg,d.phase_margin_deg,d.max_rotor_rpm,d.max_cable_speed_m_s);
fprintf('          cable travel %.6f m of %.6f m (margin %.6f); min tension %+.3e N\n', ...
    d.max_cable_travel_m,d.travel_limit_m,d.travel_margin_m,d.min_cable_tension_N);
fprintf('          max STEP/s %.2f; accel open %.2f/6000 %s, braking %.2f/2000 %s\n', ...
    d.max_n_steps_s,d.max_accel_open_steps_s2,pass(d.accel_open_ok), ...
    d.max_accel_closed_steps_s2,pass(d.accel_closed_ok));
fprintf('  PLANT   max contact force %.3f N; min bladder %.3f mL\n',d.max_contact_force_N,d.min_bladder_mL);
fprintf('  FLOW    plateau %.4f..%.4f mL/s, worst %+.2f%%, time-weighted RMSE %.4f  %s\n', ...
    d.plateau.min_Q_mL_s,d.plateau.max_Q_mL_s,d.plateau_err_pct,d.plateau_rmse_mL_s, ...
    pass(d.meets_plateau_5pct));
fprintf('  VALVE   %d open->closed, %d closed->open %s; closes %.4f s at %.4f mL/s\n', ...
    d.valve_open_to_closed,d.valve_closed_to_open,pass(d.valve_single_transition), ...
    d.t_close_s,d.q_close_left_mL_s);
fprintf('          motor stops %.4f s, tail ends %.4f s (T = %.4f s)\n', ...
    d.motor_stop_s,d.tail_end_s,d.tail_duration_s);
fprintf('  VOLUME  collected %.5f mL, error %+.5f mL (bookkeeping only)\n', ...
    d.collected_mL,d.volume_error_mL);
b = d.balance;
fprintf('  BALANCE independent flow-integral state vs bladder  %+.3e mL\n',b.state_integral_vs_bladder_mL);
fprintf('          log trapezoid vs bladder                    %+.3e mL\n',b.trapz_valve_flow_vs_bladder_mL);
fprintf('          tail integral vs drained state              %+.3e mL\n',b.tail_vs_drained_mL);
end

function s = pass(ok), if ok, s = 'ok'; else, s = '<-- FLAG'; end, end
function s = slug(x), s = regexprep(lower(x),'[^a-z0-9]+','_'); end
function s = modeName(m)
switch round(m)
    case 1, s = 'fixed step rate';
    case 2, s = 'pure PI';
    otherwise, s = 'geometric FF + PI';
end
end
