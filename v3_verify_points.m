function rep=v3_verify_points()
here=fileparts(mfilename('fullpath'));addpath(here,fullfile(here,'stage9b_port'));
a=load(fullfile(here,'validation','points_v3.mat'));P=alba9_const();B=alba9b_const();
fields={'phi','torque','Fm','Fn','fr','acc','tension','loss','slope','Qdot','pe','pg','A','D','y','dpdV','overlap','pressure'};
rtol=2e-10;rhsat=[1e-11 1e-8 1e-13 1e-13 1e-13 1e-8];
guardat=[1e-10 1e-11 1e-13 1e-11 1e-10 1e-8 1e-14];
stateScale=[.1 .01 1e-4 1e-5 1e-4 10];derivativeScale=[.01 1 1e-5 1e-3 1e-5 1];
jacat=1e-9*(derivativeScale(:)./stateScale);
worst=zeros(1,4);obsWorst=zeros(1,numel(fields));
for k=1:size(a.inputs,1)
 q=a.inputs(k,:);[dz,o]=alba9b_ref_rhs(q(1),q(2),q(3),q(4),q(5),q(6),P,B);
 oo=cellfun(@(f)o.(f),fields);err=abs(oo-a.obs(k,:))./(a.obs_atol+rtol*abs(a.obs(k,:)));
 obsWorst=max(obsWorst,err);worst(1)=max(worst(1),max(err));
 worst(2)=max(worst(2),max(abs(dz(:)'-a.rhs(k,:))./(rhsat+rtol*abs(a.rhs(k,:)))));
 J=alba9b_ref_jac(q(1),q(2),q(3),q(4),q(5),q(6),P,B);ref=squeeze(a.jac(k,:,:));
 worst(3)=max(worst(3),max(abs(J-ref)./(jacat+rtol*abs(ref)),[],'all'));
 g=alba9b_guards(q(1),q(2),q(3),q(4),q(5),q(6),P,B);
 worst(4)=max(worst(4),max(abs(g(:)'-a.guards(k,:))./(guardat+rtol*abs(a.guards(k,:)))));
end
rep=struct('points',size(a.inputs,1),'passed',all(worst<=1),'max_scaled_errors',worst, ...
 'observable_names',{fields},'observable_scaled_errors',obsWorst,'relative_tolerance',rtol, ...
 'obs_absolute_tolerance',a.obs_atol,'rhs_absolute_tolerance',rhsat,'guard_absolute_tolerance',guardat, ...
 'jacobian_absolute_tolerance',jacat,'criterion','abs(error) <= variable-specific atol + rtol*abs(reference)');
v3_write_json(fullfile(here,'validation','point_verification.json'),rep);disp(rep);assert(rep.passed,'Port verification failed');
end
