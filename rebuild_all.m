function rebuild_all(skipPythonDump)
%REBUILD_ALL Regenerate every derived artefact of Stage 9 v2 from scratch.
%
%   Order matters and is the order the request lays out:
%     1  freeze and cross-check the physical parameters
%     2  build model A (ROM3) with the latched valve
%     3  run the five cases of stage A, with the tightened re-check
%     4  verify the MATLAB port of the FULL reference against the frozen Python
%     5  build model B and run the same five cases on the full reference
%     6  figures, CSV views and the machine-readable summary
%
%   Everything it writes can be deleted and rebuilt. The only inputs are
%   ALBA_Stage8_Results.mat and, for step 4, reference_points.json (produced by
%   stage9b_port/dump_reference_points.py from the frozen Stage-8 Python source).
if nargin < 1, skipPythonDump = true; end
here = fileparts(mfilename('fullpath'));
addpath(here); addpath(fullfile(here,'stage9b_port'));
warning('off','MATLAB:print:ContentTypeImageSuggested');

cases = {'200mL_mode3_geometric_ff_pi','300mL_mode3_geometric_ff_pi', ...
         '400mL_mode3_geometric_ff_pi','300mL_mode1_fixed_step_rate', ...
         '300mL_mode2_pure_pi'};

fprintf('\n=== 1. frozen parameters ===\n');
alba9_make_params(fullfile(here,'ALBA_Stage8_Results.mat'));
alba9_param_table();

fprintf('\n=== 2. model A (ROM3) ===\n');
build_alba9_model();

fprintf('\n=== 3. stage A: program fixes, control parameters unchanged ===\n');
RA1 = run_alba9([200 300 400],3,'tight',true,'tag','A_program_fixes');
RA2 = run_alba9(300,[1 2],     'tight',true,'tag','A_program_fixes');

fprintf('\n=== 4. verify the MATLAB port of the full reference ===\n');
if ~skipPythonDump
    fprintf('  (re-dump reference_points.json with dump_reference_points.py first)\n');
end
alba9b_make_const();
alba9b_verify_port();

fprintf('\n=== 5. model B and stage B on the full reference ===\n');
build_alba9b_model();
RB1 = run_alba9b([200 300 400],3,'tag','B_full_reference');
RB2 = run_alba9b(300,[1 2],     'tag','B_full_reference');

fprintf('\n=== 6. figures, CSV views, summary ===\n');
A = fullfile(here,'results','A_program_fixes');
B = fullfile(here,'results','B_full_reference');
alba9_plot([RA1(1) RA1(2) RA1(3)], ...
    'ALBA Stage 9 v2-A: one controller, three target volumes', ...
    'results/A_program_fixes/A_volumes');
alba9_plot([RA2(1) RA2(2) RA1(2)], ...
    'ALBA Stage 9 v2-A: three controller structures at 300 mL', ...
    'results/A_program_fixes/A_controllers');
alba9b_plot([RB1(1) RB1(2) RB1(3)], ...
    'ALBA Stage 9 v2-B: full reference model, one controller, three volumes', ...
    'results/B_full_reference/B_volumes');
alba9b_plot([RB2(1) RB2(2) RB1(2)], ...
    'ALBA Stage 9 v2-B: full reference model, three controller structures at 300 mL', ...
    'results/B_full_reference/B_controllers');

for i = 1:numel(cases)
    alba9_export_csv(fullfile(A,['case_' cases{i} '.mat']), ...
                     fullfile(A,['case_' cases{i} '.csv']),@alba9_log_index);
    alba9_export_csv(fullfile(B,['case_' cases{i} '.mat']), ...
                     fullfile(B,['case_' cases{i} '.csv']),@alba9b_log_index);
end

alba9_summary_json(fullfile(here,'stage9_v2_summary.json'),cases);
fprintf('\nrebuild_all finished.\n');
end
