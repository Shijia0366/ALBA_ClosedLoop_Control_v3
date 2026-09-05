function loss = alba9_loss(Q,P)
%ALBA9_LOSS Frozen whole-waterway pressure loss. Q in m^3/s, loss in Pa.
%   Rh already contains the laminar straight-tube part, so the Darcy branch
%   only adds the excess above it. Rh is used exactly as documented by the
%   user and is NOT re-checked or re-derived anywhere in Stage 9.
%#codegen
aq = abs(Q);
Re = P.rho*aq*P.pipe_diameter/(P.water_viscosity*P.pipe_area);
loss = P.Rh*Q;
if Re > P.Re_lam
    s = min((Re-P.Re_lam)/(P.Re_turb-P.Re_lam),1);
    w = s*s*(3-2*s);
    arg = 6.9/Re + (P.roughness/(3.7*P.pipe_diameter))^1.11;
    fD = (-1.8*log10(arg))^(-2);
    turb = P.pipe_length/P.pipe_diameter*P.rho/(2*P.pipe_area^2)*fD*aq^2;
    loss = loss + sign(Q)*w*(turb - P.Rlam*aq);
end
end
