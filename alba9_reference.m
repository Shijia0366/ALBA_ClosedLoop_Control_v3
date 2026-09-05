function [Qref,r] = alba9_reference(t,discharged,valve_closed,S)
%ALBA9_REFERENCE Volume-scheduled flow reference in mL/s.
%   Startup is a raised cosine in TIME. Deceleration is scheduled on the
%   REMAINING volume before valve closure, not on a clock, so the deceleration
%   displacement is reserved automatically for every target and the profile
%   stays smooth (dQref/dr vanishes at both ends of the ramp).
%
%   198/298/398 mL are the valve-closing thresholds implied by the assumed
%   ~2 mL tail. They are NOT where deceleration starts. Every timing and shape
%   number here is a design value awaiting validation, not a measurement.
%#codegen
target = S(1); tailc = S(2); qpl = S(4); qcl = S(5);
rise   = S(6); R0    = S(7); rh   = S(8);

close_thresh = target - tailc;
r = close_thresh - discharged;          % remaining mL before the valve closes

sr = min(max(t/rise,0),1);
q_rise = 0.5*qpl*(1-cos(pi*sr));

if r >= R0
    q_dec = qpl;
elseif r <= rh
    q_dec = qcl;
else
    u = (r-rh)/(R0-rh);
    q_dec = qcl + 0.5*(qpl-qcl)*(1-cos(pi*u));
end

Qref = min(q_rise,q_dec);
if valve_closed > 0.5
    Qref = 0;
end
end
