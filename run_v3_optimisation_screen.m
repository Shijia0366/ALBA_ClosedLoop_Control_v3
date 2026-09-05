function run_v3_optimisation_screen()
% Run only after the numerical validation gate has passed. No gain tuning.
here=fileparts(mfilename('fullpath'));cd(here);addpath(here,fullfile(here,'stage9b_port'));
for mode=1:3
 p=jsondecode(fileread(fullfile(here,'validation','python_trajectories','D_diagnostic_baseline',sprintf('case_300mL_mode%d_comparison.json',mode))));
 c=jsondecode(fileread(fullfile(here,'validation',sprintf('convergence_300_mode%d',mode),'comparison.json')));
 assert(p.passed && c.passed,'Numerical validation gate failed');
end
for rise=[2.5 3]
 tag=sprintf('A_startup_%s',strrep(num2str(rise),'.','p'));
 for mode=[2 3],run_v3_case(300,mode,tag,'rise_s',rise);end
end
% Fixed-rate stays open loop. Its separately named terminal variant reduces
% the fixed rate by a remaining-volume schedule, then issues a zero request.
% 0.05/0.25 mL bracket the nominal ~0.07 mL geometric command-braking estimate
% for 300 mL PI closure (n~211 STEP/s, same open acceleration and tau_slew).
for remaining=[.05 .25]
 tag=sprintf('B_terminal_%s',strrep(num2str(remaining),'.','p'));
 for mode=1:3
  run_v3_case(300,mode,tag,'terminal_manager',1,'stop_remaining_mL',remaining,'stoptime',60);
 end
end
v3_refresh_diagnostics('D_diagnostic_baseline');v3_refresh_diagnostics('D_convergence_tight');
bdclose('all');
end
