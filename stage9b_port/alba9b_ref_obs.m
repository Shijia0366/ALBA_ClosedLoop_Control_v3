function o = alba9b_ref_obs(x,v,V,Q,count,sealed,P,B)
%ALBA9B_REF_OBS Observables of the FULL five-state reference model.
%   Ported line for line from stop_study.Plant.obs. The motor force comes from a
%   torque-PHASE model, never from a prescribed speed: the electrical phase is
%   the difference between the commanded angle and the actual cable position, so
%   the actual rotor motion is solved, not imposed.
%
%     phi    = (full_steps_per_rev/4) * (2*pi*count/steps_per_rev - x/r)
%     torque = torque_cap * sin(phi)          conditional flat envelope
%     acc    = (Fm - damping*v - friction - Fn) / total_equiv_mass
%
%   The legacy 15 N force_max plays no part here and must not be reintroduced.
%#codegen
[D,A,pe,pg,y,dpdV] = alba9b_geometry(V,P);

overlap = max(x - y, 0);
Fn = P.contact_stiffness*overlap;

phi    = B.phase_factor*(2*pi*count/B.steps_per_rev - x/B.radius_m);
torque = B.torque_cap_Nm*sin(phi);
Fm     = torque/B.radius_m;
fr     = B.friction*tanh(v/B.friction_velocity);
acc    = (Fm - B.damping*v - fr - Fn)/B.total_equiv_mass_kg;
tension = B.nonmotor_equiv_mass_kg*acc + B.damping*v + fr + Fn;

drive = pe + pg + Fn/A + P.head_Pa;
[loss,slope] = alba9b_loss(Q,P);
if sealed > 0.5
    qdot = 0;
else
    qdot = (drive - loss)/P.Ih;
end

o = struct('D',D,'A',A,'pe',pe,'pg',pg,'y',y,'dpdV',dpdV,'overlap',overlap, ...
           'Fn',Fn,'phi',phi,'torque',torque,'Fm',Fm,'fr',fr,'acc',acc, ...
           'tension',tension,'drive',drive,'loss',loss,'slope',slope, ...
           'Qdot',qdot,'pressure',drive-P.head_Pa,'omega',v/B.radius_m);
end
