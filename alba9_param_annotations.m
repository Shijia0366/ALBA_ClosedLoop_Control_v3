function A = alba9_param_annotations()
%ALBA9_PARAM_ANNOTATIONS Unit, source and assumption status for every parameter.
%
%   STATUS vocabulary - deliberately narrow, so nothing can quietly be promoted
%   to "measured" that was never measured:
%
%     given       a design / nameplate / user-supplied value. Taken as given for
%                 this project; NOT independently verified by any measurement in
%                 Stage 7, 8 or 9.
%     derived     computed from other entries by a stated formula. No new information.
%     assumption  chosen by the modeller or inherited as a conditional prior.
%                 Not identified, not measured.
%     learned     fitted by the frozen Stage-8 SINDy run against SYNTHETIC
%                 noise-free reference data. Not a measured material constant.
%     solver      a numerical setting, not physics.
%     legacy      present in the file but NOT used by Stage 9, and in some cases
%                 actively dangerous to reuse. Read the note.
%     metadata    a string / flag describing scope, not a number.
%
%   Nothing in this table is an experimentally validated hardware value. The
%   Stage-8 package states this explicitly and Stage 9 does not change it.
row = @(n,u,s,src,note) struct('name',n,'unit',u,'status',s,'source',src,'note',note);

A = [
% ---------------------------------------------------------------- geometry
row('V0','m^3','given','user / design', ...
    'Initial bladder water 500 mL.')
row('Dref_outer','m','given','user / mould drawing', ...
    'Natural (unpressurised) OUTER diameter 70 mm.')
row('Dref_inner','m','given','user / mould drawing', ...
    'Natural INNER diameter 66 mm. The 2 mm per-side wall implied by 70-66 is the MOULD wall thickness; it is not the wall thickness once the bladder is filled.')
row('Dmin','m','given','user / design', ...
    'Minimum allowed bladder outer diameter 50 mm. A domain guard, not a plant coefficient.')
row('cable_gap','m','given','user / design', ...
    'Radial clearance from the cable to the bladder outer wall, 3 mm. Enters the geometry through D0/L0/xmax, not the ROM3 right-hand side directly.')
row('Vsilicone','m^3','derived','shell volume from Dref_outer/Dref_inner', ...
    'Silicone shell volume, added to the water volume before the equivalent sphere diameter is taken.')
row('Vnatural','m^3','derived','from Dref_inner', ...
    'Natural internal cavity volume. Below it the spherical-shell extrapolation is unvalidated.')
row('D0','m','derived','equivalent sphere at V0', ...
    'Outer diameter at the initial 500 mL, 0.10034658601199008 m. y = pi*(D0-D(V)) is the cable coordinate.')
row('L0','m','derived','cable path at V0','Reference cable length.')
row('xmax','m','derived','cable travel guard','Travel limit used by the full reference model guards.')
row('Vmin','m^3','derived','sphere at Dmin', ...
    'Water volume at the 50 mm minimum outer diameter, 36.388 mL. Domain guard.')
% ---------------------------------------------------------------- material
row('material','-','metadata','user confirmed','Ecoflex 00-20.')
row('material_confirmed_by_user','-','metadata','user confirmed','Material identity confirmed; its constants are still not measured on this build.')
row('stress100','Pa','given','material datasheet', ...
    'Stress at 100% elongation, used to set the Neo-Hookean mu.')
row('mu','Pa','derived','from stress100', ...
    'Neo-Hookean shear modulus. Datasheet-derived, not measured on the actual bladder.')
% ---------------------------------------------------------------- payload
row('cage_mass','kg','given','user / design', ...
    'Hoberman sphere mass 0.270 kg.')
row('support_fraction','-','assumption','modeller', ...
    '1.0: the bladder carries the cage entirely.')
row('gravity_height_slope','-','assumption','modeller', ...
    'Conditional centroid-height law. pg already accounts for the cage weight; do NOT add a separate self-weight term.')
row('selfweight_model','-','metadata','modeller','Name of the above law.')
% ---------------------------------------------------------------- fluid path
row('rho','kg/m^3','given','water at room temperature','')
row('gravity','m/s^2','given','standard','')
row('water_viscosity','Pa*s','given','water at room temperature','')
row('head','m','given','user / rig layout', ...
    'Bladder outlet is +0.10 m ABOVE the final outlet, so rho*g*h = +981 Pa ASSISTS the flow. The old snapshot value -0.1 is overridden at load time.')
row('pipe_diameter','m','given','user / tubing spec','4 mm internal bore.')
row('pipe_segments','m','given','user / rig layout','0.03 + 0.03 + 0.30 m.')
row('pipe_length','m','derived','sum(pipe_segments)','0.36 m of straight tube.')
row('pipe_area','m^2','derived','pi*d^2/4','')
row('tube_volume','m^3','derived','pipe_area*pipe_length','4.5239 mL of primed tube water.')
row('Rh','Pa*s/m^3','given','bench effective resistance', ...
    ['BENCH EFFECTIVE RESISTANCE of the whole waterway at low flow: tube, valve, ' ...
     'flowmeter and fittings lumped together. It is NOT a 4 mm outlet resistance ' ...
     'and not an identified valve curve. Stage 9 uses it exactly as supplied and ' ...
     'does not re-derive or re-check its calibration basis.'])
row('R_tube','Pa*s/m^3','derived','128*mu_w*L/(pi*d^4)', ...
    'Theoretical laminar straight-tube resistance. Already CONTAINED in Rh, which is why alba9_loss subtracts it before adding the Darcy branch.')
row('I_tube','Pa*s^2/m^3','derived','rho*L/A','Hydraulic inertance of the known straight tube.')
row('quadratic_loss','Pa*s^2/m^6','assumption','disabled', ...
    'Zero. Not identified. Zero does not prove the hardware has no local nonlinear loss.')
% ---------------------------------------------------------------- actuator
row('shaft_diameter','m','given','user / motor spec','5 mm bare shaft.')
row('effective_radius','m','given','user / motor spec', ...
    'Bare winding radius 0.0025 m.')
row('mass_equivalent','kg','assumption','modeller', ...
    'Non-motor equivalent translating mass 0.5 kg. Not measured.')
row('damping','N*s/m','assumption','modeller','20 N*s/m. Not identified.')
row('friction','N','assumption','modeller','0.1 N Coulomb level. Not identified.')
row('friction_velocity','m/s','assumption','modeller', ...
    '1e-5 m/s regularisation velocity for the tanh friction law. Numerical, not physical.')
row('contact_stiffness','N/m','assumption','modeller', ...
    '2e5 N/m button-to-bladder contact stiffness. Not identified.')
row('force_max','N','legacy','old fixture limit', ...
    ['15 N. A LEGACY force-fixture limit. It is NOT the motor capability bound ' ...
     'for this stage and must never be substituted for the torque-phase model.'])
row('force_ramp','s','legacy','old fixture ramp','Not used by Stage 9.')
row('initial_contact_fact','-','metadata','user confirmed', ...
    'The top button is always in contact with the bladder.')
row('water_column_primed_assumption','-','assumption','modeller', ...
    'Tubing is primed and bubble-free at t = 0.')
row('pressure_reference','-','metadata','modeller','Lumped bladder-outlet gauge pressure.')
% ---------------------------------------------------------------- solver
row('rtol','-','solver','Stage-8 reference','')
row('atol','-','solver','Stage-8 reference','')
row('max_step','s','solver','Stage-8 reference','')
row('time_horizon','s','solver','Stage-8 reference','')
];
end
