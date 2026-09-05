function names = alba9b_guard_names()
%ALBA9B_GUARD_NAMES Names of the seven model-validity guards, in order.
%   Kept out of the codegen path so the guard vector itself stays a plain
%   fixed-size numeric signal.
names = {'phase_boundary','50mm_limit','sphere_limit','extension_boundary', ...
         'speed_envelope_boundary','cable_tension_boundary','reverse_flow_boundary'};
end
