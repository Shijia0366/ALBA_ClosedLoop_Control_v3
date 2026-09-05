function rebuild_v3(pythonExe)
% One entry point: rebuild and ACTUALLY RUN Simulink, then independent Python.
% Python 3.13 with requirements-v3.txt must already be installed; no auto-install.
% Example: rebuild_v3('D:\conda\python.exe')
% Run in a separate fresh copy to retain the delivered evidence unchanged.
if nargin<1,pythonExe='python';end
here=fileparts(mfilename('fullpath'));cd(here);addpath(here,fullfile(here,'stage9b_port'));
assert(license('test','Simulink')==1,'Simulink license unavailable; nothing validated');
py('v3_verify.py',{'audit'});py('v3_verify.py',{'points'});
v3_prepare();v3_verify_points();
for m=1:3
 run_v3_case(300,m,'D_diagnostic_baseline');
 run_v3_case(300,m,'D_convergence_tight','tight',true);
 py('v3_verify.py',{'trajectory',fullfile('results','D_diagnostic_baseline',sprintf('case_300mL_mode%d.mat',m))});
end
py('v3_post.py',{'checks'});run_v3_optimisation_screen();
for m=[2 3]
 run_v3_case(300,m,'B_guard_localised','terminal_manager',1,'stop_remaining_mL',.25,'stoptime',60);
end
run_v3_final_validation();bdclose('all');
for target=[300 200 400]
 py('v3_verify.py',{'trajectory',fullfile('results','FINAL_verified',sprintf('case_%dmL_mode2.mat',target))});
end
for m=[2 3]
 py('v3_verify.py',{'trajectory',fullfile('results','B_guard_localised',sprintf('case_300mL_mode%d.mat',m))});
 py('v3_post.py',{'compare',fullfile('frozen_v2_ROM_results',sprintf('case_300mL_mode%d.mat',m)), ...
 fullfile('results','D_diagnostic_baseline',sprintf('case_300mL_mode%d.mat',m)),sprintf('ROM_full_300_mode%d',m),'rom'});
end
py('v3_post.py',{'compare',fullfile('results','FINAL_verified','case_300mL_mode2.mat'), ...
 fullfile('results','FINAL_verified_tight','case_300mL_mode2.mat'),'final_verified_convergence_300','tight'});
py('v3_post.py',{'plots','D_diagnostic_baseline','A_startup_2p5','A_startup_3','B_terminal_0p05','B_terminal_0p25','B_guard_localised','C_combined_300','FINAL_verified','FINAL_verified_tight'});
py('v3_post.py',{'poster','FINAL_verified','2'});py('v3_deliver.py',{'report'});
fprintf('Finished actual Simulink runs. Read REPORT_v3.md for pass/fail, not just this completion message.\n');
 function py(script,args)
  quote=@(s)['"' char(s) '"'];cmd=[quote(pythonExe) ' ' quote(fullfile(here,script))];
  for j=1:numel(args),cmd=[cmd ' ' quote(args{j})];end %#ok<AGROW>
  [status,output]=system(cmd,'-echo');
  assert(status==0,'ALBA:PythonVerification','Python step failed: %s',output);
 end
end
