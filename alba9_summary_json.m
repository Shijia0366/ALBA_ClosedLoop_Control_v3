function alba9_summary_json(outfile,cases)
%ALBA9_SUMMARY_JSON One machine-readable summary of both stages.
%   Readable without MATLAB. Carries the scope statements alongside the numbers
%   so they cannot be separated from them.
here = fileparts(mfilename('fullpath'));
A = fullfile(here,'results','A_program_fixes');
B = fullfile(here,'results','B_full_reference');

out = struct();
out.package = 'ALBA Stage 9 v2';
out.baseline = 'ALBA_Step9_ClosedLoop_v1, preserved unchanged';
out.stage_A = 'program and statistics fixes only; control parameters identical to v1';
out.stage_B = 'the SAME controller on the full five-state reference; no retuning';
out.stage_C = 'not run in this package';
out.model_reidentified = false;
out.basis_terms_added = false;
out.Rh_re_verified = false;
out.Rh_note = 'bench effective resistance, used exactly as documented; never adjusted to improve a curve';
out.experimental_validation = false;
out.measurement_noise_in_loop = false;
out.controller_qualified = false;

pv = fullfile(here,'stage9b_port','port_verification.json');
if isfile(pv), out.port_verification = jsondecode(fileread(pv)); end
pc = fullfile(here,'parameter_consistency.json');
if isfile(pc)
    d = jsondecode(fileread(pc));
    out.parameter_consistency = d.consistency;
end

rows = {};
for i = 1:numel(cases)
    r = struct('case',cases{i});
    fa = fullfile(A,['case_' cases{i} '.mat']);
    fb = fullfile(B,['case_' cases{i} '.mat']);
    if isfile(fa)
        a = load(fa,'diag'); r.rom3 = a.diag;
    end
    if isfile(fb)
        b = load(fb,'diag'); r.full_reference = b.diag;
    end
    rows{end+1} = r; %#ok<AGROW>
end
out.cases = rows;

out.headline_findings = { ...
 'Stage A: the valve latch removes 12 and 21 closed->open transitions that v1 logged; they were zero-duration and had not corrupted v1 results.'; ...
 'Stage A: the plateau error is 2.38 percent, not the 0.96 percent v1 reported. v1 skipped the first 0.2 s of the plateau, which is where the overshoot sits.'; ...
 'Stage B: peak torque reaches 98.0 percent of the ASSUMED flat envelope during the 2 s startup ramp. 15 mL/s is a near-limit design case and is not declared feasible.'; ...
 'Stage B: the fixed-step-rate termination demands 96.3 percent of the envelope AFTER the valve has closed, while no longer delivering any water.'; ...
 'The final volume error is bookkeeping only: an ideal volume event closes the valve and the delivered tail equals the compensated tail.'};

fid = fopen(outfile,'w','n','UTF-8');
c = onCleanup(@() fclose(fid));
fwrite(fid,jsonencode(out,'PrettyPrint',true),'char');
clear c
fprintf('  %s\n',outfile);
end
