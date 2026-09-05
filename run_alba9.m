function R = run_alba9(targets,modes,varargin)
%RUN_ALBA9 Stage-9 v2 runner. One independently named result file per case.
%
%   R = run_alba9(300,3)
%   R = run_alba9([200 300 400],3,'tag','A_program_fixes')
%   R = run_alba9(300,[1 2 3],'tight',true,'Kp',80)
%
%   v2 CHANGES
%     * every case is saved to its own file, results/<tag>/case_<name>.mat, and a
%       matching .log. v1 wrote one alba9_results.mat that consecutive calls
%       silently overwrote.
%     * each case is cross-checked against the independent ode15s realisation on
%       the WHOLE trajectory (flow, volume, force, command), not only the final
%       volume, and optionally re-run at tighter tolerance and smaller max step.
%
%   This is a simulation of a REDUCED surrogate of the conditional reference
%   model. It is not experimental validation, and the controller is not yet
%   qualified: the same controller still has to be replayed on the full
%   five-state reference with its motor torque, phase, travel and stall guards.
if nargin < 1 || isempty(targets), targets = 300; end
if nargin < 2 || isempty(modes),   modes   = [1 2 3]; end

opt = struct('tag','A_program_fixes','tight',false);
scenArgs = {};
k = 1;
while k <= numel(varargin)
    if isfield(opt,varargin{k})
        opt.(varargin{k}) = varargin{k+1};
    else
        scenArgs(end+1:end+2) = varargin(k:k+1); %#ok<AGROW>
    end
    k = k + 2;
end

here = fileparts(mfilename('fullpath'));
addpath(here);
outdir = fullfile(here,'results',opt.tag);
if ~isfolder(outdir), mkdir(outdir); end

mdl = 'alba9_closed_loop';
if ~isfile(fullfile(here,[mdl '.slx'])), build_alba9_model(mdl); end
load_system(fullfile(here,[mdl '.slx']));

P = alba9_const();
assignin('base','V0mL',P.V0*1e6);

R = struct('name',{},'scenario',{},'diag',{},'t',{},'y',{},'ode',{},'cmp',{},'cmp_tight',{});
for target = targets(:).'
    for mode = modes(:).'
        [S,S_vec] = alba9_scenario(target,mode,scenArgs{:});
        name = sprintf('%dmL_mode%d_%s',target,mode,slug(modeName(mode)));
        assignin('base','S_vec',S_vec);

        set_param(mdl,'RelTol','2e-8','AbsTol','1e-7','MaxStep','0.025');
        out = sim(mdl);
        L = out.alba9_log;
        t = L.time;
        y = squeeze(L.signals.values);
        d = alba9_diagnose(t,y,S);

        ode = alba9_check_ode(S);
        d.q0_exact_left_limit_mL_s = ode.q0_exact_left_limit;
        d.q0_capture_error_vs_exact_mL_s = d.q0_captured_mL_s - ode.q0_exact_left_limit;
        cmp = alba9_compare(t,y,ode,'independent ode15s, event-based closure');

        cmpT = [];
        if opt.tight
            set_param(mdl,'RelTol','1e-10','AbsTol','1e-9','MaxStep','0.0125');
            outT = sim(mdl);
            LT = outT.alba9_log;
            odeT = alba9_check_ode(S,struct('tight',true));
            cmpT = alba9_compare(LT.time,squeeze(LT.signals.values),odeT, ...
                                 'tightened solver, both implementations');
            set_param(mdl,'RelTol','2e-8','AbsTol','1e-7','MaxStep','0.025');
        end

        rec = struct('name',name,'scenario',S,'diag',d,'t',t,'y',y, ...
                     'ode',rmfield(ode,{'sol_open','sol_closed'}), ...
                     'cmp',cmp,'cmp_tight',cmpT);
        R(end+1) = rec; %#ok<AGROW>

        save(fullfile(outdir,['case_' name '.mat']),'-struct','rec');
        writeLog(fullfile(outdir,['case_' name '.log']),d,cmp,cmpT);
        report(d,cmp,cmpT);
    end
end
fprintf('\n%d cases written to %s (one .mat and one .log each)\n',numel(R),outdir);
end

% -------------------------------------------------------------------------
function report(d,cmp,cmpT)
fprintf('\n=== %d mL | %s ===\n',d.target_mL,d.mode_name);
fprintf('  VALVE   open->closed %d, closed->open %d   %s   (raw comparator: %d / %d)\n', ...
    d.valve_open_to_closed,d.valve_closed_to_open,pass(d.valve_single_transition), ...
    d.raw_comparator_open_to_closed,d.raw_comparator_closed_to_open);
fprintf('          closes %.4f s, latch lag %.1f us = %.2e mL of extra flow\n', ...
    d.t_close_s,d.latch_lag_s*1e6,d.latch_lag_volume_mL);
fprintf('          q0 captured %.5f, exact left limit %.5f, difference %+.2e mL/s\n', ...
    d.q0_captured_mL_s,d.q0_exact_left_limit_mL_s,d.q0_capture_error_vs_exact_mL_s);
fprintf('          motor stops %.4f s, tail ends %.4f s (T = %.4f s, extrapolated rule)\n', ...
    d.motor_stop_s,d.tail_end_s,d.tail_duration_s);
fprintf('  VOLUME  collected %.5f mL, error %+.5f mL  (BOOKKEEPING ONLY)\n', ...
    d.collected_mL,d.volume_error_mL);
fprintf('          real tail 1.0/1.5/2.5/3.0 mL would give %+.1f/%+.1f/%+.1f/%+.1f mL\n', ...
    d.tail_sensitivity(1).final_volume_error_mL,d.tail_sensitivity(2).final_volume_error_mL, ...
    d.tail_sensitivity(4).final_volume_error_mL,d.tail_sensitivity(5).final_volume_error_mL);
fprintf('  FLOW    startup      RMSE %8.4f (time-weighted) over %.2f s\n', ...
    d.startup.rmse_time_weighted_mL_s,d.startup.duration_s);
fprintf('          plateau      RMSE %8.4f (time-weighted) vs %8.4f (unweighted, v1 method)\n', ...
    d.plateau.rmse_time_weighted_mL_s,d.plateau.rmse_unweighted_mL_s);
fprintf('          plateau      %.4f..%.4f mL/s, worst %+.2f%% (target 5%%)   %s\n', ...
    d.plateau_min_mL_s,d.plateau_max_mL_s,d.plateau_err_pct,pass(d.meets_plateau_5pct));
fprintf('          deceleration RMSE %8.4f (time-weighted) over %.2f s\n', ...
    d.deceleration.rmse_time_weighted_mL_s,d.deceleration.duration_s);
fprintf('  MOTOR   max %.2f STEP/s (ceiling 3000); accel open %.2f/6000 %s, braking %.2f/2000 %s\n', ...
    d.max_n_steps_s,d.max_accel_open_steps_s2,pass(d.accel_open_ok), ...
    d.max_accel_closed_steps_s2,pass(d.accel_closed_ok));
fprintf('          saturated %.4f s, slew limited %.4f s, final integral state %.3f STEP/s\n', ...
    d.n_saturated_time_s,d.slew_limited_time_s,d.final_ni);
fprintf('  ROM     max contact force %.3f N %s (audited to 118 N), min stiffness %.0f N/m %s\n', ...
    d.max_F_N,pass(d.F_within_audited_118N),d.min_keff,pass(d.min_keff>0));
fprintf('          min bladder volume %.3f mL %s\n',d.min_V_mL,pass(d.V_above_Vmin));
b = d.balance;
fprintf('  BALANCE valve-flow integral vs bladder decrease  %+.3e mL\n',b.valve_flow_vs_bladder_mL);
fprintf('          tail integral vs drained state          %+.3e mL\n',b.tail_vs_drained_mL);
fprintf('          total inventory (incl. %.4f mL tube)    %+.3e mL\n',b.tube_initial_mL,b.total_inventory_error_mL);
fprintf('          accounting identity (cannot fail)       %+.3e mL\n',b.accounting_identity_error_mL);
fprintf('  VS ODE  dt_close %.2e s, excluded window %.2e s; whole-trajectory dQ %.2e mL/s, dV %.2e mL, dF %.2e N, dn %.2e STEP/s\n', ...
    cmp.t_close_difference_s,cmp.excluded_window_s,cmp.worst_flow_mL_s,cmp.worst_volume_mL, ...
    cmp.worst_force_N,cmp.worst_command_sps);
if ~isempty(cmpT)
    fprintf('  TIGHT   dt_close %.2e s; dQ %.2e mL/s, dV %.2e mL, dF %.2e N, dn %.2e STEP/s\n', ...
        cmpT.t_close_difference_s,cmpT.worst_flow_mL_s,cmpT.worst_volume_mL, ...
        cmpT.worst_force_N,cmpT.worst_command_sps);
end
end

function writeLog(file,d,cmp,cmpT)
fid = fopen(file,'w','n','UTF-8');
c = onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(struct('diagnostics',d,'vs_ode',cmp,'vs_ode_tight',cmpT), ...
        'PrettyPrint',true));
end

function s = pass(ok)
if ok, s = 'ok'; else, s = '<-- FLAG'; end
end

function s = modeName(m)
switch round(m)
    case 1, s = 'fixed step rate';
    case 2, s = 'pure PI';
    otherwise, s = 'geometric FF + PI';
end
end

function s = slug(x)
s = regexprep(lower(x),'[^a-z0-9]+','_');
end
