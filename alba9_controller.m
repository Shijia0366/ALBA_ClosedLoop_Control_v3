function [n_in,dni,n_raw,e_track,n_ff,sat_flag,slew_flag] = ...
        alba9_controller(Qref,Qm,Vm,n,ni,valve_closed,S)
%ALBA9_CONTROLLER Positional geometric feedforward + PI, back-calculation
%anti-windup, and command shaping. Output is a STEP/s command, never a flow.
%
%   STRUCTURE (unchanged from v1; do not silently swap it for another form):
%     positional PI      n_raw = n_ff + Kp*e + ni,  d(ni)/dt = Ki*e + Kaw*(n - n_raw)
%     command dynamics   n is a STATE. Its input is slew limited to the
%                        conditional acceleration envelope and the state itself
%                        is clamped to the conditional STEP/s ceiling, so the
%                        rate limit and the saturation act on exactly the signal
%                        the anti-windup term measures against.
%
%   UNITS: e is in mL/s and the output in STEP/s. Kp and Ki carry that
%   conversion. Do not additionally convert the error to SI.
%
%   mode 1 = fixed step rate (open loop), 2 = pure PI, 3 = feedforward + PI.
%
%   v2 FIX: mode 1 previously reported e = 0, which wrote a real tracking error
%   away as zero in the log. The true tracking error is now ALWAYS returned as
%   e_track. Mode 1 remains open loop: e_track is reported, never fed back.
%
%   First round: V_est is the ideal volume and Qm is ideal flow feedback. No
%   sensor range, quantisation, noise or delay yet, and the real contact force
%   is never used as a control signal.
%#codegen
Kp   = S(9);  Ki   = S(10); Kaw  = S(11); nmax = S(12);
amax = S(13); abrk = S(14); tau  = S(15); mode = S(16); nfix = S(17);

P = alba9_const();
[~,Ax,~,~] = alba9_geometry(Vm*1e-6,P);

e_track = Qref - Qm;                 % always the true tracking error, logged
n_ff    = (Qref*1e-6)/(P.ku*Ax);     % geometric inverse: Q = ku*n*Ax at steady state

if mode < 1.5
    n_raw = nfix;                    % open loop: e_track is measured, not used
    e_ctrl = 0;
elseif mode < 2.5
    e_ctrl = e_track;
    n_raw  = Kp*e_ctrl + ni;
else
    e_ctrl = e_track;
    n_raw  = n_ff + Kp*e_ctrl + ni;
end

% v3: optional, separately labelled termination manager. Disabled for the
% diagnostic baseline; it never turns fixed-rate mode into flow feedback.
terminal_stop = false;
if numel(S) >= 21 && S(20) > 0.5
    remaining = S(1)-S(2)-(P.V0*1e6-Vm);
    if mode < 1.5
        u = min(max((remaining-S(8))/(S(7)-S(8)),0),1);
        qdec = S(5)+0.5*(S(4)-S(5))*(1-cos(pi*u));
        n_raw = nfix*qdec/S(4);
    end
    terminal_stop = remaining <= S(21);
end
if valve_closed > 0.5
    n_des = 0;      a = abrk;        % closed-valve brake, still a finite deceleration
elseif terminal_stop
    n_des = 0;      a = amax;
else
    n_des = n_raw;  a = amax;
end
% Compare BEFORE and AFTER target clipping, even when n_dot is already zero.
sat_flag = double(n_des > nmax)-double(n_des < 0);
n_des = min(max(n_des,0),nmax);

n_want = (n_des-n)/tau;
n_in   = min(max(n_want,-a),a);
slew_flag = double(abs(n_want) > a + 1e-9);

% The command state is clamped by the integrator; report when it is riding a
% limit AND being pushed further into it, which is the case anti-windup exists for.

if mode < 1.5 || valve_closed > 0.5 || terminal_stop
    dni = -Kaw*ni;                   % no integral action open loop or after sealing
else
    dni = Ki*e_ctrl + Kaw*(n-n_raw); % back-calculation against the APPLIED command
end
end
