function v3_refresh_diagnostics(tag)
% Recompute diagnostics only; never alter saved simulation states or controls.
here=fileparts(mfilename('fullpath'));addpath(here,fullfile(here,'stage9b_port'));
files=dir(fullfile(here,'results',tag,'case_*.mat'));
for k=1:numel(files)
 if contains(files(k).name,'failure'),continue;end
 file=fullfile(files(k).folder,files(k).name);r=load(file,'t','y','scenario','solver','wall_seconds');
 diag=alba9v3_diagnose(r.t,r.y,r.scenario,120);save(file,'diag','-append');
 v3_export_csv(r.t,r.y,strrep(file,'.mat','.csv'));
 v3_write_json(strrep(file,'.mat','.json'),struct('scenario',r.scenario,'diagnostics',diag,'solver',r.solver,'wall_seconds',r.wall_seconds));
end
end
