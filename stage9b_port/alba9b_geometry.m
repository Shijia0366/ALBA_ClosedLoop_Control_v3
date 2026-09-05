function [D,A,pe,pg,y,dpdV] = alba9b_geometry(V,P)
%ALBA9B_GEOMETRY Frozen reference geometry, ported from stop_study.Plant.geometry.
%   Same law the ROM uses, plus the two quantities the full model needs and the
%   ROM does not: the cable coordinate y and dp/dV for the Jacobian.
%   V in m^3. P is alba9_const().
%#codegen
ri = (3*V/(4*pi))^(1/3);
D  = (6*(V+P.Vsilicone)/pi)^(1/3);
A  = D*D/2;
ua = P.Dref_inner/(2*ri);
ub = P.Dref_outer/D;
pe = -2*P.mu*(ua + ua^4/4 - ub - ub^4/4);
pg = P.cage_mass*P.gravity*P.gravity_height_slope/(pi*A);
dpdV = (2*P.mu/3)*(ua*(1+ua^3)/V - ub*(1+ub^3)/(V+P.Vsilicone)) ...
       - 2*pg/(3*(V+P.Vsilicone));
y = pi*(P.D0 - D);
end
