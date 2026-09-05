function dz = rom3_rhs(t,z,u,p,h,r,steps_per_rev,coef,powers)
% Physics-informed SINDy ROM. SI states: [contact force N; water m^3; Q m^3/s].
% Open-valve, taut-cable operating branch only; not a stepper-stall model.
% MATLAB source supplied, but was not executed in MATLAB in this delivery.
F=max(z(1),0); V=z(2); Q=z(3);
if V<=0, error('ROM outside positive-volume domain'); end
D=(6*(V+p.Vsilicone)/pi)^(1/3); A=D^2/2;
ua=p.Dref_inner/(2*(3*V/(4*pi))^(1/3)); ub=p.Dref_outer/D;
pe=-2*p.mu*(ua+ua^4/4-ub-ub^4/4);
pg=p.cage_mass*p.gravity*p.gravity_height_slope/(pi*A);
k=1e5*sum(coef(:).*(F/100).^double(powers(:)));
if k<=0, error('Learned stiffness outside positive-stiffness domain'); end
ku=2*pi*r/double(steps_per_rev);
dF=k*(ku*u(t)-Q/A);
if z(1)<=0 && dF<0, dF=0; end
aq=abs(Q); Re=p.rho*aq*p.pipe_diameter/(p.water_viscosity*h.pipe_area_m2);
loss=p.Rh*Q;
if Re>h.Re_laminar
    s=min((Re-h.Re_laminar)/(h.Re_turbulent-h.Re_laminar),1);
    w=s*s*(3-2*s);
    arg=6.9/Re+(h.pipe_absolute_roughness_m/(3.7*p.pipe_diameter))^1.11;
    f=(-1.8*log10(arg))^(-2);
    turb=h.pipe_length_m/p.pipe_diameter*p.rho/(2*h.pipe_area_m2^2)*f*aq^2;
    loss=loss+sign(Q)*w*(turb-h.R_pipe_laminar_Pa_s_m3*aq);
end
head=p.rho*p.gravity*p.head;
dz=[dF;-Q;(pe+pg+F/A+head-loss)/h.Ih_Pa_s2_m3];
end
