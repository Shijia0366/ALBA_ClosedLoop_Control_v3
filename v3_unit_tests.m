function report=v3_unit_tests()
% Diagnostics test the unclipped target, including the equilibrium at a rail.
[S,v]=alba9_scenario(300,2); %#ok<ASGLU>
[dn,~,~,~,~,flag]=alba9_controller(100,0,500,3000,0,0,v);
assert(flag==1 && dn==0,'Upper target clipping must be flagged at equilibrium');
[dn,~,~,~,~,flag]=alba9_controller(0,10,500,0,0,0,v);
assert(flag==-1 && dn==0,'Lower target clipping must be flagged at equilibrium');
[~,~,~,~,~,flag]=alba9_controller(10,10,500,100,100,0,v);assert(flag==0);
[dn,~,~,~,~,flag]=alba9_controller(100,0,500,1000,0,1,v);
assert(dn==-2000 && flag==0,'Closure is a brake request, not saturation of PI');
P=alba9_const();B=alba9b_const();V=300e-6;
[~,~,~,~,yy]=alba9b_geometry(V,P);
o=alba9b_ref_obs(yy-1e-4,0,V,0,0,0,P,B);assert(o.Fn==0);
o=alba9b_ref_obs(yy+116/P.contact_stiffness,0,V,0,0,0,P,B);
assert(abs(o.Fn-116)<1e-8,'Force must not be clipped to legacy 15 N');
o=alba9b_ref_obs(yy+120/P.contact_stiffness,0,V,0,0,0,P,B);
assert(abs(o.Fn-120)<1e-8,'118 N must not clip the full model');
assert(P.Rh==1.6440083045256217e9 && B.torque_cap_Nm==0.29658517619749203);
assert(S.Kp==50 && S.Ki==600 && S.Kaw==10 && S.terminal_manager==0);
assert(strcmp(v3_classify_status(true,false,false),'normal_completion'));
assert(strcmp(v3_classify_status(true,true,false),'boundary_triggered'));
assert(strcmp(v3_classify_status(false,false,true),'timeout'));
assert(strcmp(v3_classify_status(false,false,false,'solver_failure'),'solver_failure'));
assert(strcmp(v3_classify_status(false,false,false,'manual_interruption'),'manual_interruption'));
assert(~strcmp(v3_classify_status(false,false,false),'normal_completion'));
report=struct('passed',15,'failed',0,'scope','diagnostics, completion classification, unilateral contact and frozen constraints');
v3_write_json(fullfile('validation','unit_tests.json'),report);
disp(report);
end
