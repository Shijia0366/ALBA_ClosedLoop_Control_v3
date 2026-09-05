function [loss,slope] = alba9b_loss(Q,P)
%ALBA9B_LOSS Frozen waterway loss AND its analytic slope, ported from
%stop_study.Plant.loss. The slope is what the reference Jacobian needs.
%   Identical law to alba9_loss; Rh is used exactly as supplied and the laminar
%   straight-tube part it already contains is subtracted before the Darcy branch
%   is added, so nothing is counted twice.
%#codegen
aq = abs(Q);
re = P.rho*P.pipe_diameter/(P.water_viscosity*P.pipe_area)*aq;
Rlam = P.Rlam;
loss = P.Rh*Q;
slope = P.Rh;
if re > P.Re_lam
    span = P.Re_turb - P.Re_lam;
    s = min((re-P.Re_lam)/span,1);
    w  = s*s*(3-2*s);
    dw = 6*s*(1-s)/span;
    arg = 6.9/re + (P.roughness/(3.7*P.pipe_diameter))^1.11;
    den = -1.8*log10(arg);
    fD  = den^-2;
    dfD = -2*den^-3*(1.8*6.9/(log(10)*arg*re^2));
    pipeC = P.pipe_length/P.pipe_diameter*P.rho/(2*P.pipe_area^2);
    re_factor = P.rho*P.pipe_diameter/(P.water_viscosity*P.pipe_area);
    turb  = pipeC*fD*aq*aq;
    dturb = pipeC*(2*fD*aq + dfD*re_factor*aq*aq);
    loss  = loss + sign(Q)*w*(turb - Rlam*aq);
    slope = slope + dw*re_factor*(turb - Rlam*aq) + w*(dturb - Rlam);
end
end
