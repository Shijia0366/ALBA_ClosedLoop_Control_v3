function rec=run_v3_case(target,mode,tag,varargin)
% Actual Simulink runner; complete MAT is primary, CSV is an additional view.
here=fileparts(mfilename('fullpath'));addpath(here,fullfile(here,'stage9b_port'));
opt=struct('tight',false,'stoptime',120);args={};
for k=1:2:numel(varargin)
 if isfield(opt,varargin{k}),opt.(varargin{k})=varargin{k+1};else,args(end+1:end+2)=varargin(k:k+1);end %#ok<AGROW>
end
[S,S_vec]=alba9_scenario(target,mode,args{:});P=alba9_const();
dirout=fullfile(here,'results',tag);if ~isfolder(dirout),mkdir(dirout);end
name=sprintf('case_%dmL_mode%d',target,mode);file=fullfile(dirout,name);
diary([file '.run.log']);dc=onCleanup(@()diary('off'));
fprintf('SIMULINK RUN %s %s | %s\n',tag,name,datestr(now,30));
mdl='alba9b_full_loop';load_system(fullfile(here,'stage9b_port',[mdl '.slx']));
assignin('base','S_vec',S_vec);assignin('base','V0_m3',P.V0);
rtol=2e-8;maxstep=.025;factor=1;if opt.tight,rtol=2e-10;maxstep=.005;factor=.01;end
ints={'count','x','v','V','Q','intQ','n','ni','latch','tc','q0','tail','quiet'};
at=[1e-6 2e-11 2e-10 2e-13 2e-14 2e-13 1e-4 1e-4 1e-9 1e-7 1e-7 1e-9 1e-8];
for j=1:numel(ints),set_param([mdl '/Int_' ints{j}],'AbsoluteTolerance',num2str(at(j)*factor,17));end
set_param(mdl,'RelTol',num2str(rtol,17),'MaxStep',num2str(maxstep,17),'StopTime',num2str(opt.stoptime));
rec=struct('name',name,'scenario',S,'t',[],'y',[],'diag',struct(),'solver',struct('rtol',rtol,'maxstep',maxstep,'atol',at*factor,'horizon_s',opt.stoptime));
v3_write_json([file '.status.json'],struct('status','running','start',datestr(now,30)));
lastwarn('');tic;
try
 out=sim(mdl,'CaptureErrors','on');rec.wall_seconds=toc;rec.simulation_metadata=out.SimulationMetadata;
 if ~isempty(out.ErrorMessage),error('ALBA:SolverFailure','%s',out.ErrorMessage);end
 L=out.alba9b_log;rec.t=L.time;rec.y=squeeze(L.signals.values);
 if size(rec.y,1)~=numel(rec.t),rec.y=rec.y.';end
 rec.diag=alba9v3_diagnose(rec.t,rec.y,S,opt.stoptime);rec.simulation_output=out;
 [rec.last_warning,rec.last_warning_id]=lastwarn;
 save([file '.mat'],'-struct','rec','-v7');
 v3_write_json([file '.json'],struct('scenario',S,'diagnostics',rec.diag,'solver',rec.solver,'wall_seconds',rec.wall_seconds));
 v3_write_json([file '.status.json'],struct('status',rec.diag.status,'end',datestr(now,30)));
 fprintf('%s | torque %.9f Nm (%.3f%%), Q RMSE %.6f, close %.6f, actual stop %.6f, volume %.8f\n', ...
 rec.diag.status,rec.diag.max_abs_torque_Nm,rec.diag.torque_utilisation_pct,rec.diag.plateau_rmse_mL_s, ...
 rec.diag.t_close_s,rec.diag.mechanical_stop_s,rec.diag.collected_mL);
 v3_export_csv(rec.t,rec.y,[file '.csv']);
catch ex
 status='implementation_failure';
 if strcmp(ex.identifier,'ALBA:SolverFailure'),status='solver_failure';end
 if contains(lower(ex.identifier),'interrupt')||contains(lower(ex.identifier),'terminated'),status='manual_interruption';end
 rec.diag=struct('status',status,'completed',false,'message',ex.message,'identifier',ex.identifier);
 if exist('out','var'),rec.simulation_output=out;end
 save([file '_failure.mat'],'-struct','rec','-v7');v3_write_json([file '.status.json'],rec.diag);
 fprintf(2,'%s\n',getReport(ex,'extended'));close_system(mdl,0);rethrow(ex);
end
end
