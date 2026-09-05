function g = alba9b_guards(x,v,V,Q,count,sealed,P,B)
%ALBA9B_GUARDS Model validity guards, ported from stop_study.Plant.guard_values.
%   Each entry is positive while the model is inside its declared domain and
%   crosses zero when it leaves. They are NOT hardware safety limits; they mark
%   where this conditional model stops being licensed.
%
%     phase_boundary   |phi| < pi/2. Beyond it the torque-phase relation turns
%                      over and the step is being pulled out - the frozen model
%                      does not represent stall past this point.
%     50mm_limit       cable travel below xmax
%     sphere_limit     water volume above Vmin
%     extension        x >= 0
%     speed_envelope   |v/r| below 2*pi rad/s
%     cable_tension    tension >= 0; a cable cannot push
%     reverse_flow     Q >= 0 (open valve only)
%#codegen
o = alba9b_ref_obs(x,v,V,Q,count,sealed,P,B);
g = [pi/2 - abs(o.phi);
     B.xmax - x;
     V - P.Vmin;
     x + 1e-6;
     2*pi - abs(v/B.radius_m);
     o.tension + 1e-5;
     Q + 1e-13];
% All seven are always returned and the vector keeps a FIXED size, exactly as
% the frozen guard_values does. Only the first six are ARMED while the valve is
% sealed, because reverse flow is meaningless once the valve has cut the flow to
% zero; the caller disarms the seventh by raising it, never by resizing the
% vector - a variable-size signal would defeat Simulink's dimension inference.
% Guard names live in alba9b_guard_names, kept out of this codegen path.
end
