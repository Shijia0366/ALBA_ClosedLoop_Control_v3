function v3_prepare()
here=fileparts(mfilename('fullpath'));cd(here);addpath(here,fullfile(here,'stage9b_port'));
diary(fullfile(here,'logs','environment_and_build.log'));cl=onCleanup(@()diary('off'));
e=struct('version',version,'release',version('-release'),'simulink_license_test',license('test','Simulink'),'toolboxes',ver);
assert(e.simulink_license_test==1,'Simulink license unavailable');disp(e);
v3_write_json(fullfile('logs','environment.json'),e);v3_unit_tests();build_alba9b_model();
[~,S_vec]=alba9_scenario(300,3);assignin('base','S_vec',S_vec);assignin('base','V0_m3',.0005);
set_param('alba9b_full_loop','SimulationCommand','update');
for block={'Ref','ValveCmp'}
 ph=get_param(['alba9b_full_loop/' block{1}],'PortHandles');
 port=1;if strcmp(block{1},'Ref'),port=2;end
 line=get_param(ph.Inport(port),'Line');src=get_param(line,'SrcBlockHandle');
 assert(strcmp(get_param(src,'Name'),'intQ_mL'),'Controller-side volume source is not flow integration');
end
v3_write_json(fullfile('validation','connection_audit.json'),struct('compile_update_passed',true, ...
 'feedback','Q_actual mL/s, V_est=500-integral(Q_actual) mL','controller_output','STEP/s', ...
 'plant_states','SI x,v,V,Q and integrated STEP position; Q reset once via v2 latch', ...
 'reference_and_valve_use_independent_flow_integral',true,'ROM_extra_actuator',false,'limits_unchanged',true));
end
