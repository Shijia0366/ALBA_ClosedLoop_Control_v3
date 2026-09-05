function alba9_export_json(outfile,varargin)
%ALBA9_EXPORT_JSON Write every run's diagnostics to JSON.
%   So the results can be read and checked without MATLAB. Time series stay in
%   the .mat files; this is the summary table plus the declared scope.
if nargin < 1 || isempty(outfile)
    outfile = fullfile(fileparts(mfilename('fullpath')),'stage9_diagnostics.json');
end

runs = {};
for k = 1:numel(varargin)
    R = varargin{k};
    for i = 1:numel(R)
        d = R(i).diag;
        s = R(i).scenario;
        d.tail_sensitivity = struct( ...
            'actual_tail_mL',[d.tail_sensitivity.actual_tail_mL], ...
            'final_volume_error_mL',[d.tail_sensitivity.final_volume_error_mL]);
        runs{end+1} = struct('scenario',s,'diagnostics',d); %#ok<AGROW>
    end
end

P = alba9_const();
out = struct();
out.what_this_is = ['Stage-9 closed-loop simulation of the FROZEN Stage-8 three-state ' ...
    'ROM3 surrogate, with ideal noise-free flow feedback.'];
out.what_this_is_not = { ...
    'not experimental validation'; ...
    'not a noisy or sensor-in-the-loop result'; ...
    'not a qualified controller: the same controller has not yet been replayed on the full five-state reference'; ...
    'not a check of motor torque, electrical phase, rotor speed, travel or stall margin - ROM3 has no motor states'; ...
    'none of the acceptance numbers is a medical or regulatory criterion'};
out.model_reidentified = false;
out.basis_terms_added = false;
out.Rh_Pa_s_m3 = P.Rh;
out.Rh_provenance = 'Used exactly as documented by the user. Stage 9 does not re-derive or re-check it.';
out.stiffness_coefficients = P.coef;
out.stiffness_powers = P.powers;
out.stiffness_audited_force_range_N = [0 118];
out.volume_error_caveat = ['The final volume error is near zero BY CONSTRUCTION: the valve ' ...
    'closes on an exact ideal volume event and the tail delivered equals the tail compensated. ' ...
    'See tail_sensitivity for the error that actually matters.'];
out.runs = runs;

fid = fopen(outfile,'w');
cleaner = onCleanup(@() fclose(fid));
fwrite(fid,jsonencode(out,'PrettyPrint',true),'char');
clear cleaner
fprintf('Wrote %s\n',outfile);
end
