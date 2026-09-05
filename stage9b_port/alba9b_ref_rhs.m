function [dz,o] = alba9b_ref_rhs(x,v,V,Q,count,sealed,P,B)
%ALBA9B_REF_RHS Right-hand side of the FULL reference model, states
%[x; v; V; Q; integral_Q; energy_work], ported from stop_study.Plant.rhs.
%#codegen
o = alba9b_ref_obs(x,v,V,Q,count,sealed,P,B);
power = o.Fm*v + P.head_Pa*Q - B.damping*v*v - o.fr*v - o.loss*Q;
dz = [v; o.acc; -Q; o.Qdot; Q; power];
end
