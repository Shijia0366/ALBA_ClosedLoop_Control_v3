function [df,dVm,dQm,Ax,keff,F] = alba9_plant(f,Vm,Qm,n,valve_closed,P,clamp_contact)
%ALBA9_PLANT Frozen three-state ROM3. States f = F/100 [N], Vm [mL], Qm [mL/s].
%   n is the commanded STEP/s. Flow is an ODE state solved by this object; it
%   is never assigned to the reference. Open-valve, taut-cable, contacting
%   branch only. This is not a stall, reverse-winding or lost-contact model.
%
%   CLAMP_CONTACT selects HOW the unilateral contact condition f >= 0 is
%   realised, not WHETHER it holds:
%     true  (default) - clamp dF/dt inside the equation, as the frozen Stage-8
%                       code does. Correct for a plain ODE integration.
%     false           - leave dF/dt alone and let a SATURATED integrator hold
%                       the state at zero. Required in Simulink: the in-equation
%                       clamp is a step discontinuity in f, and ode15s builds a
%                       numerical Jacobian across it, which stalls the state
%                       near zero instead of releasing it.
%   Both realise the same physics; only the numerics differ.
%#codegen
if nargin < 7
    clamp_contact = true;
end
V = Vm*1e-6;
Q = Qm*1e-6;
F = max(f*100,0);

[~,Ax,pe,pg] = alba9_geometry(V,P);

keff = 1e5*sum(P.coef.*(F/100).^P.powers);
rate = keff*(P.ku*n - Q/Ax);
if clamp_contact && f <= 0 && rate < 0
    rate = 0;                 % unilateral contact: no tensile "pulling" force
end

df  = rate/100;
dVm = -Qm;
if valve_closed > 0.5
    dQm = 0;                  % ideal instantaneous closure; Q is reset to zero
else
    dQm = ((pe + pg + F/Ax + P.head_Pa - alba9_loss(Q,P))/P.Ih)*1e6;
end
end
