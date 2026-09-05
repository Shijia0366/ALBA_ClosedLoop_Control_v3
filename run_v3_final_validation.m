function run_v3_final_validation()
% Select using 300 mL only; then freeze and evaluate held-out volumes.
here=fileparts(mfilename('fullpath'));cd(here);addpath(here,fullfile(here,'stage9b_port'));
opts={'rise_s',3,'terminal_manager',1,'stop_remaining_mL',.05};
rec=run_v3_case(300,2,'C_combined_300',opts{:});
assert(rec.diag.completed && rec.diag.meets_plateau_5pct && rec.diag.meets_volume_2mL, ...
 'Combined 300 mL did not pass; do not silently proceed with another design');
selected=struct('selection_basis','300 mL only, before held-out 200/400 mL', ...
 'reason','Pure PI 3 s had lowest startup torque among plateau-passing screened controllers; 0.05 mL early stop reduced post-close compression. 0.25 mL PI variants failed cable-tension boundary.', ...
 'scenario_300',rec.scenario,'gains_retuned',false,'physical_parameters_changed',false);
v3_write_json(fullfile('validation','selected_controller.json'),selected);
for target=[300 200 400]
 run_v3_case(target,2,'FINAL_verified',opts{:});
end
% Numerical convergence on the final 300 mL combined design, not a retune.
run_v3_case(300,2,'FINAL_verified_tight',opts{:},'tight',true);
end
