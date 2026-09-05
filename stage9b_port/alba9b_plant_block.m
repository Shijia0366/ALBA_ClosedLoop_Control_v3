function [dx,dv,dV,dQ,Fn,phi,torque,omega,tension,acc,violations,gidx,stop,g,overlap_raw,pressure] = ...
        alba9b_plant_block(x,v,V,Q,count,valve_closed)
%ALBA9B_PLANT_BLOCK Simulink face of the verified FULL five-state reference plant.
%
%   Flow is solved from the motor input, the actual mechanical motion and the
%   hydraulic equation together. Nothing here imposes a flow and nothing sets the
%   actual speed equal to the commanded speed: the motor acts only through the
%   torque-phase relation, phi = (Nfull/4)*(theta_cmd - x/r), and the resulting
%   force fights damping, friction and the bladder contact.
%
%   The guards are evaluated but NOT enforced: this function reports where the
%   model leaves its declared domain, it does not clip anything to keep the run
%   looking healthy. The caller stops the simulation on a violation.
%#codegen
P = alba9_const();
B = alba9b_const();

[dz,o] = alba9b_ref_rhs(x,v,V,Q,count,valve_closed,P,B);
dx = dz(1);  dv = dz(2);  dV = dz(3);  dQ = dz(4);

Fn = o.Fn;  phi = o.phi;  torque = o.torque;
omega = o.omega;  tension = o.tension;  acc = o.acc;

g = alba9b_guards(x,v,V,Q,count,valve_closed,P,B);
if valve_closed > 0.5
    g(7) = inf;        % disarm reverse flow once sealed; do NOT resize the vector
end
% Never rank rad, m, m^3, rad/s, N and m^3/s against one another.
violations = 0; gidx = 0;
for k=1:7
    if g(k)<=0
        violations = violations+1;
        if gidx==0, gidx=k; end
    end
end
stop = double(violations > 0);
overlap_raw = x-o.y; % signed gap is diagnostic; unilateral law is unchanged
pressure = o.pressure;
end
