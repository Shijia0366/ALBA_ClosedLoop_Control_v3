function [q_tail,done,q_outlet] = alba9_tail(t,valve_closed,q0,tclose,drained,n,Qm,S)
%ALBA9_TAIL Downstream residual-water tail and the FINAL-OUTLET flow.
%
%   WHAT IS GIVEN AND WHAT IS ASSUMED
%     given by the user   the tail is about 2 mL, and about 1.33 s of drain was
%                         observed at one operating point.
%     assumed here        the half-cosine SHAPE, and the rule T = 2*V_tail/q0
%                         which makes the tail start at the pre-closure flow and
%                         end smoothly with integral exactly V_tail.
%
%   Those two assumptions plus V_tail = 2 mL only reproduce 1.33 s when
%   q0 is about 3 mL/s. At a different closing flow the rule EXTRAPOLATES:
%   closing at 8.66 mL/s gives T = 0.462 s, which is a consequence of the
%   assumed rule, not an observed drain time. The four conditions "starts at the
%   pre-closure flow", "half cosine", "2 mL" and "1.33 s" cannot all hold at
%   every operating point, and this function does not pretend they do.
%
%   The tail drains post-valve tubing. It never reduces bladder volume a second
%   time and it feeds back into nothing upstream.
%
%   q_outlet is the flow at the FINAL outlet: the valve flow while the valve is
%   open, the tail drain after it closes. It is a different signal from the
%   valve flow and both are reported.
%#codegen
Vt = S(3);

q_tail = 0;
if valve_closed > 0.5 && q0 > 0 && drained < Vt
    T  = 2*Vt/q0;
    tt = min(max(t-tclose,0),T);
    q_tail = 0.5*q0*(1+cos(pi*tt/T));
end

if valve_closed > 0.5
    q_outlet = q_tail;
else
    q_outlet = Qm;
end

done = 0;
if valve_closed > 0.5 && drained >= Vt-1e-6 && n <= 1e-3
    done = 1; % ROM/legacy command-and-tail flag only; NOT actual mechanical stop
end
end
