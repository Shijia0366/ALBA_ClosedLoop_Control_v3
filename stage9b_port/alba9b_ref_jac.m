function J = alba9b_ref_jac(x,v,V,Q,count,sealed,P,B)
%ALBA9B_REF_JAC Analytic Jacobian of the full reference model, ported from
%stop_study.Plant.jac. Supplied to the stiff solver for the same reason the
%frozen model supplies it: the regularised friction and the contact branch are
%both sharp, and a finite-difference Jacobian across them is unreliable.
%#codegen
o = alba9b_ref_obs(x,v,V,Q,count,sealed,P,B);

A = o.A;
if o.overlap > 0
    k = P.contact_stiffness;
else
    k = 0;
end
fnV = k/A;
dFm_dx = -B.torque_cap_Nm*cos(o.phi)*B.phase_factor/B.radius_m^2;
dfr = B.friction/B.friction_velocity*(1 - tanh(v/B.friction_velocity)^2);

J = zeros(6,6);
J(1,2) = 1;
J(2,1:3) = [(dFm_dx-k)/B.total_equiv_mass_kg, ...
            (-B.damping-dfr)/B.total_equiv_mass_kg, ...
            -fnV/B.total_equiv_mass_kg];
J(3,4) = -1;
if sealed <= 0.5
    dP = o.dpdV + fnV/A - o.Fn*(2/(pi*o.D))/A^2;
    J(4,[1 3 4]) = [k/A/P.Ih, dP/P.Ih, -o.slope/P.Ih];
end
J(5,4) = 1;
J(6,[1 2 4]) = [v*dFm_dx, ...
                o.Fm - 2*B.damping*v - o.fr - v*dfr, ...
                P.head_Pa - o.loss - Q*o.slope];
end
