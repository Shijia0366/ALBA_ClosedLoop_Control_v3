function alba9_make_params(stage8_mat)
%ALBA9_MAKE_PARAMS Generate alba9_const.m from the FROZEN Stage-8 deliverable.
%   Nothing here re-identifies anything. Coefficients, powers, Rh, geometry
%   and hydraulics are copied verbatim out of ALBA_Stage8_Results.mat, which
%   was exported from the frozen results/selected_rom3.json.
if nargin < 1
    stage8_mat = fullfile(fileparts(mfilename('fullpath')),'ALBA_Stage8_Results.mat');
end
S = load(stage8_mat);
p = S.plant_parameters;
h = S.hydraulic_parameters;

P = struct();
P.Vsilicone            = double(p.Vsilicone);
P.Dref_inner           = double(p.Dref_inner);
P.Dref_outer           = double(p.Dref_outer);
P.mu                   = double(p.mu);
P.cage_mass            = double(p.cage_mass);
P.gravity              = double(p.gravity);
P.gravity_height_slope = double(p.gravity_height_slope);
P.rho                  = double(p.rho);
P.water_viscosity      = double(p.water_viscosity);
P.pipe_diameter        = double(p.pipe_diameter);
P.pipe_area            = double(h.pipe_area_m2);
P.pipe_length          = double(h.pipe_length_m);
P.Rh                   = double(p.Rh);
P.Rlam                 = double(h.R_pipe_laminar_Pa_s_m3);
P.Ih                   = double(h.Ih_Pa_s2_m3);
P.Re_lam               = double(h.Re_laminar);
P.Re_turb              = double(h.Re_turbulent);
P.roughness            = double(h.pipe_absolute_roughness_m);
P.head_Pa              = double(p.rho)*double(p.gravity)*double(p.head);
P.ku                   = 2*pi*double(S.motor_radius_m)/double(S.steps_per_rev);
P.D0                   = double(p.D0);
P.Vmin                 = double(p.Vmin);
P.V0                   = double(p.V0);
P.contact_stiffness    = double(p.contact_stiffness);
P.coef                 = reshape(double(S.stiffness_coefficients),1,[]);
P.powers               = reshape(double(S.stiffness_powers),1,[]);

% ---- the frozen numbers must still be what Stage 8 reported ----
assert(abs(P.Rh-1644008304.5256217) < 1e-3, 'Rh does not match the frozen value');
assert(abs(P.head_Pa-981) < 1e-6, 'head must be the +10 cm assisting head');
assert(isequal(P.powers,[0 1 2 3 4 6]), 'stiffness powers changed');
assert(nnz(P.coef) == 4, 'the frozen support must keep exactly 4 terms');
assert(abs(P.ku-4.908738521234052e-06) < 1e-15, 'ku changed');

folder = fileparts(mfilename('fullpath'));
fid = fopen(fullfile(folder,'alba9_const.m'),'w');
cleaner = onCleanup(@() fclose(fid));
fprintf(fid,'function P = alba9_const()\n');
fprintf(fid,'%%ALBA9_CONST Frozen Stage-8 constants. GENERATED FILE, do not hand-edit.\n');
fprintf(fid,'%%   Source: %s\n',strrep(stage8_mat,'\','/'));
fprintf(fid,'%%   These values are copied, never refitted. Rh is used as documented\n');
fprintf(fid,'%%   and is not re-derived in Stage 9.\n');
fprintf(fid,'%%#codegen\n');
f = fieldnames(P);
for i = 1:numel(f)
    v = P.(f{i});
    if isscalar(v)
        fprintf(fid,'P.%s = %.17g;\n',f{i},v);
    else
        parts = arrayfun(@(x)sprintf('%.17g',x),v,'UniformOutput',false);
        fprintf(fid,'P.%s = [%s];\n',f{i},strjoin(parts,' '));
    end
end
fprintf(fid,'end\n');
clear cleaner
save(fullfile(folder,'alba9_params.mat'),'P');
fprintf('alba9_const.m written with %d frozen constants.\n',numel(f));
end
