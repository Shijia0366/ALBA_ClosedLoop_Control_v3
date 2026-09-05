function A = alba9_hydraulic_annotations()
%ALBA9_HYDRAULIC_ANNOTATIONS Annotations for the hydraulic configuration block.
%   Same status vocabulary as alba9_param_annotations.
row = @(n,u,s,src,note) struct('name',n,'unit',u,'status',s,'source',src,'note',note);

A = [
row('version','-','metadata','Stage-8','Conditional reference, version tag.')
row('pipe_area_m2','m^2','derived','pi*d^2/4','')
row('pipe_length_m','m','derived','sum of straight segments','0.36 m.')
row('R_pipe_laminar_Pa_s_m3','Pa*s/m^3','derived','128*mu_w*L/(pi*d^4)', ...
    'Contained inside Rh. Subtracted before the Darcy branch is added, so the laminar tube loss is never counted twice.')
row('R_residual_Pa_s_m3','Pa*s/m^3','derived','Rh - R_pipe_laminar', ...
    'Lumped remainder (valve, flowmeter, fittings, and any absorbed bias). A conditional decomposition of a fitted total, NOT an independently measured component resistance.')
row('Ih_Pa_s2_m3','Pa*s^2/m^3','derived','rho*L/A', ...
    'Known straight-tube inertance only. Excludes unknown component inertance.')
row('hydraulic_inertia_scale','-','assumption','modeller','1.0.')
row('linear_hydraulic_time_constant_s','s','derived','Ih/Rh','')
row('nonlinear_pipe','-','assumption','modeller','Darcy branch enabled above Re_laminar.')
row('pipe_absolute_roughness_m','m','assumption','smooth-pipe assumption', ...
    'Zero. A smooth-pipe assumption, not a measured roughness.')
row('Re_laminar','-','assumption','engineering threshold', ...
    '2000. Not a measured transition for this rig.')
row('Re_turbulent','-','assumption','engineering threshold', ...
    '4000. Cubic smooth blending between the two, also an assumption.')
row('additional_component_K_Pa_s2_m6','Pa*s^2/m^6','assumption','disabled', ...
    'Zero and unauthorised. Stage 9 does not add any unidentified valve or flowmeter nonlinear coefficient.')
row('additional_component_K_authorized','-','metadata','Stage-8','false.')
row('material_relaxation_enabled','-','metadata','Stage-8','false; no invented Ecoflex relaxation constants.')
row('material_relaxation_source','-','metadata','Stage-8', ...
    'Not identified; no Ecoflex relaxation constants were invented.')
row('valve_scope','-','metadata','Stage-8','Fully open; stop at a closing command.')
row('fluid_scope','-','metadata','Stage-8','Primed, bubble-free, rigid bore, forward liquid flow only.')
row('inertia_source','-','metadata','Stage-8','')
row('roughness_source','-','metadata','Stage-8','')
row('transition_source','-','metadata','Stage-8','')
row('component_loss_source','-','metadata','Stage-8','')
row('additional_component_K_source','-','metadata','Stage-8','')
row('experimental_hardware_validation','-','metadata','Stage-8','false.')
];
end
