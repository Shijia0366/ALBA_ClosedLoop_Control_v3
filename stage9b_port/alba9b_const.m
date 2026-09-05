function B = alba9b_const()
%ALBA9B_CONST Frozen FULL five-state reference constants. GENERATED.
%   Read out of the frozen Python model, never hand-copied.
%
%   ALL of the mechanical values below are ASSUMPTIONS, not identified
%   or measured quantities: equivalent masses, damping, friction, the
%   friction smoothing velocity and the contact stiffness.
%   torque_cap is a CONDITIONAL FLAT ENVELOPE, not a measured 12 V
%   running-torque curve. 0.67 A is an assumed/estimated current, not a
%   measured limit. legacy_force_max_N is an old fixture limit and must
%   NEVER be used in place of the torque-phase model.
%#codegen
B.total_equiv_mass_kg = 1.02675137617256;
B.rotor_inertia_kg_m2 = 3.2921961010785002e-06;
B.nonmotor_equiv_mass_kg = 0.5;
B.torque_cap_Nm = 0.29658517619749203;
B.radius_m = 0.0025000000000000001;
B.steps_per_rev = 3200;
B.full_steps_per_rev = 200;
B.microsteps_per_full_step = 16;
B.supply_voltage_V = 12;
B.current_assumed_A = 0.67000000000000004;
B.current_limit_measured = 0;
B.torque_curve_matches_drive_conditions = 0;
B.damping = 20;
B.friction = 0.10000000000000001;
B.friction_velocity = 1.0000000000000001e-05;
B.contact_stiffness = 200000;
B.head_Pa = 981;
B.Ih = 28647889.756541163;
B.Rh = 1644008304.5256217;
B.xmax = 0.15816846474859467;
B.Vmin = 3.6388020508979353e-05;
B.D0 = 0.10034658601199008;
B.Vsilicone = 2.9061826440808009e-05;
B.legacy_force_max_N_DO_NOT_USE = 15;
B.max_rate = 3000;
B.max_accel = 6000;
B.closed_brake_accel = 2000;
B.phase_factor = 50;   % full_steps_per_rev/4
end
