function alba9b_make_const(jsonfile)
%ALBA9B_MAKE_CONST Generate alba9b_const.m from the frozen Python reference dump.
%   The actuator and mechanical constants are read out of the dump produced by
%   dump_reference_points.py, which constructed stop_study.Plant exactly as
%   Stage 8 does. Nothing is hand-copied.
here = fileparts(mfilename('fullpath'));
if nargin < 1 || isempty(jsonfile)
    jsonfile = fullfile(here,'reference_points.json');
end
d = jsondecode(fileread(jsonfile));
m = d.meta;

% cross-check the values the request states explicitly
assert(abs(m.rotor_inertia_kg_m2-3.2921961010785e-6) < 1e-21,'rotor inertia changed');
assert(abs(m.nonmotor_equiv_mass_kg-0.5) == 0,'non-motor equivalent mass changed');
assert(abs(m.total_equiv_mass_kg-1.02675137617256) < 5e-15,'total equivalent mass changed');
assert(abs(m.total_equiv_mass_kg - (m.nonmotor_equiv_mass_kg + ...
        m.rotor_inertia_kg_m2/m.radius_m^2)) < 1e-15, ...
        'total equivalent mass must be the sum, counted once');
assert(abs(m.torque_cap_Nm-0.29658517619749203) == 0,'torque envelope changed');
assert(m.supply_voltage_V == 12 && abs(m.current_assumed_A-0.67) < 1e-15,'supply changed');
assert(~m.current_limit_measured && ~m.torque_curve_matches_drive_conditions, ...
       'the torque envelope must stay flagged as unmeasured');
assert(m.damping == 20 && m.friction == 0.1 && m.friction_velocity == 1e-5, ...
       'mechanical assumptions changed');
assert(m.contact_stiffness == 2e5,'contact stiffness changed');
assert(m.steps_per_rev == 3200 && m.full_steps_per_rev == 200 && ...
       m.microsteps_per_full_step == 16,'stepping changed');
assert(abs(m.head_Pa-981) < 1e-9,'head changed');

fid = fopen(fullfile(here,'alba9b_const.m'),'w');
c = onCleanup(@() fclose(fid));
fprintf(fid,'function B = alba9b_const()\n');
fprintf(fid,'%%ALBA9B_CONST Frozen FULL five-state reference constants. GENERATED.\n');
fprintf(fid,'%%   Read out of the frozen Python model, never hand-copied.\n');
fprintf(fid,'%%\n');
fprintf(fid,'%%   ALL of the mechanical values below are ASSUMPTIONS, not identified\n');
fprintf(fid,'%%   or measured quantities: equivalent masses, damping, friction, the\n');
fprintf(fid,'%%   friction smoothing velocity and the contact stiffness.\n');
fprintf(fid,'%%   torque_cap is a CONDITIONAL FLAT ENVELOPE, not a measured 12 V\n');
fprintf(fid,'%%   running-torque curve. 0.67 A is an assumed/estimated current, not a\n');
fprintf(fid,'%%   measured limit. legacy_force_max_N is an old fixture limit and must\n');
fprintf(fid,'%%   NEVER be used in place of the torque-phase model.\n');
fprintf(fid,'%%#codegen\n');
f = fieldnames(m);
for i = 1:numel(f)
    v = m.(f{i});
    if isnumeric(v) && isscalar(v)
        fprintf(fid,'B.%s = %.17g;\n',f{i},double(v));
    elseif islogical(v)
        fprintf(fid,'B.%s = %d;\n',f{i},double(v));
    end
end
% phi = (full_steps_per_rev/4) * (theta_cmd - x/r); keep the factor derived
fprintf(fid,'B.phase_factor = %.17g;   %% full_steps_per_rev/4\n',m.full_steps_per_rev/4);
fprintf(fid,'end\n');
clear c
fprintf('alba9b_const.m written.\n');
end
