function [D,Ax,pe,pg] = alba9_geometry(V,P)
%ALBA9_GEOMETRY Frozen spherical-shell geometry, Neo-Hookean pe and cage pg.
%   V in m^3. Identical to the Stage-8 reference; nothing is refitted here.
%   Ax is |dV/d(cable coordinate)|, not a projected sphere area.
%#codegen
D  = (6*(V+P.Vsilicone)/pi)^(1/3);
Ax = D*D/2;
ua = P.Dref_inner/(2*(3*V/(4*pi))^(1/3));
ub = P.Dref_outer/D;
pe = -2*P.mu*(ua+ua^4/4-ub-ub^4/4);
pg = P.cage_mass*P.gravity*P.gravity_height_slope/(pi*Ax);
end
